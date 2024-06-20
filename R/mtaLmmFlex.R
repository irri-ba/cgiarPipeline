mtaLmmFlex <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    analysisId=NULL,
    analysisIdForGenoModifications=NULL,
    inputFormulation=NULL, 
    envsToInclude=NULL,
    trait= NULL, traitFamily=NULL, useWeights=TRUE,
    heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf,
    modelType="blup", # either "blup", "pblup", "gblup", "rrblup"
    nMarkersRRBLUP=1000,
    deregress=FALSE,
    maxIters=50, batchSizeToPredict=500, tolParInv=1e-4,
    verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING LMM SOLVER
  mtaAnalysisId <- as.numeric(Sys.time())
  '%!in%' <- function(x,y)!('%in%'(x,y)) 
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  if(modelType %in% c("gblup","ssgblup","rrblup") & is.null(phenoDTfile$data$geno)){
    stop("Please include marker information in your data structure to fit this model type.", call. = FALSE)
  }
  # if user will use markers read them and apply the modifications in the table
  if(modelType %in% c("gblup","ssgblup","rrblup") & !is.null(phenoDTfile$data$geno)){
    Markers <- phenoDTfile$data$geno
    if(is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
      if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and you have not provided a modifications table to impute the genotype data. Please go to the 'Markers QA/QC' module prior to run a gBLUP or rrBLUP model.", call. = FALSE)}
    }else{ # user provided a modifications Id
      modificationsMarkers <- phenoDTfile$modifications$geno
      theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
      if(length(theresMatch) > 0){ # there's a modification file after matching the Id
        modificationsMarkers <- modificationsMarkers[theresMatch,]
        Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
      }else{ # there's no match of the modification file
        if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
      }
    }
    Markers <- Markers[,sample(1:min(c(ncol(Markers), nMarkersRRBLUP)))] # don't use all the markers if goes beyond nK
    Markers <- Markers - 1 # center markers now
  }
  if(is.null(phenoDTfile$metadata$weather)){ # avoid an error when there is no weather information
    phenoDTfile$metadata$weather <- cgiarBase::create_getData_object()$metadata$weather
  }
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait))
  heritUB <- rep(heritUB,length(trait))
  meanLB <- rep(meanLB,length(trait))
  meanUB <- rep(meanUB,length(trait))
  traitOrig <- trait
  
  conditionsModel <- unlist( lapply(inputFormulation, function(x){
    if( (!is.null(x$left)) & (x$center == "|" | x$center == "||") & ("designation" %in% x$right) ){
      return(TRUE)
    }else{return(FALSE)}
  }) )
  if(sum(conditionsModel) == 0){modelType <- "blue"} # assign blue method if the user provided designation in the fixed part
  ############################
  # loading the dataset
  # mydata <- phenoDTfile$predictions #
  # if (nrow(mydata) < 2) stop("Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.", call. = FALSE)
  # mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  
  # add the other available columns to the dataset
  ff <- cgiarBase::formLme4(input0=inputFormulation,object=phenoDTfile, analysisId=analysisId)      
  mydata <<- ff$predictions
  
  # if the user provides two ids with same trait and environments kill the job
  allTraitsInMyData <- unique(na.omit(mydata$trait))
  for(iTrait in allTraitsInMyData){ # iTrait = allTraitsInMyData[1]
    traitByIdCheck <- with(mydata[mydata$trait == iTrait, ], table(environment, analysisId))
    traitByIdCheck <- traitByIdCheck/traitByIdCheck; 
    checkOnPreds <- apply(traitByIdCheck,1,sum, na.rm=TRUE)
    badIdSelection <- which( checkOnPreds > 1)
    if(length(badIdSelection) > 0){
      stop(paste( "You have selected multiple analysisId to be analyzed together but trait",iTrait,"has data for the same environments", paste(names(checkOnPreds)[badIdSelection], collapse =", ") ,"in the IDs provided.") , call. = FALSE)
    }
  }
  
  # some checks after filtering
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  # loading the metrics
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId),]
  #############################
  # defining the name of the surrogate depending on the model
  surrogate <- c("TGV","EBV","GEBV","GEBV","ssGEBV","PV"); names(surrogate) <- c("blup","pblup","gblup","rrblup","ssgblup","blue")
  # if user didn't provide a table for which environments should be included, make it! include all environments as default
  if(is.null(envsToInclude)){
    envsToInclude=  as.data.frame( do.call( rbind, list (with(mydata, table(environment,trait)) ) ) )
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if(nrow(bad) > 0){envsToInclude[bad] = 0}
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }
  # check if traits specified are actually available and remove the ones that are not present
  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  traitToRemove <- character()
  if(!is.null(envsToInclude)){
    for(k in 1:length(trait)){
      if( length(which(envsToInclude[,trait[k]] > 0)) == 0 ){
        if(verbose){ cat(paste0("'", trait[k], "' has not envs selected. It will be removed from trait list \n"))}
        traitToRemove <- c(traitToRemove,trait[k])
      }
    }
    trait <- setdiff(trait,traitToRemove)
  }
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  meanLB <- meanLB[which(traitOrig %in% trait)]
  meanUB <- meanUB[which(traitOrig %in% trait)]
  modelTypeTrait <- rep(modelType, length(trait)); names(modelTypeTrait) <- trait
  ##############################
  ##############################
  ## met analysis
  ##############################
  ##############################
  library(lme4breeding)
  allEnvironments <- na.omit(unique(mydata[,"environment"]))
  predictionsList <- list(); counter=counter2=1
  traitToRemove <- character()
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    failedMarkerModel=FALSE
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    mydataSub <- droplevels(mydataSub[which(!is.na(mydataSub$predictedValue)),])
    # 
    if(!is.null(envsToInclude)){
      goodFieldsUser = rownames(envsToInclude)[which(envsToInclude[,iTrait] > 0)]
    }else{goodFieldsUser <- na.omit(unique(mydataSub[,"environment"]))}
    # remove bad environment based on h2 and r2
    pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2","r2")),]
    goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"environment"])
    goodFields <- intersect(goodFields, goodFieldsUser)
    mydataSub <- mydataSub[which(mydataSub$environment %in% goodFields),]
    # remove bad environment based on environment means
    pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("mean")),]
    if(nrow(pipeline_metricsSub) > 0){ # second reduction for good environments or a given threshold
      goodFieldsMean <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > meanLB[counter2]) & (pipeline_metricsSub$value < meanUB[counter2])),"environment"])
      goodFields <- intersect(goodFields, goodFieldsMean)
      mydataSub <- mydataSub[which(mydataSub$environment %in% goodFields),]
    }
    if(verbose){print(paste("Fields included:",paste(goodFields,collapse = ",")))}
    ## remove records without marker data if marker effects
    if(modelTypeTrait[iTrait] == "rrblup"){ # if rrBLUP model we remove records without markers
      mydataSub <- mydataSub[which(mydataSub$designation %in% rownames(Markers)),]
      if(nrow(mydataSub) == 0 & verbose){
        cat("There was no match between markers and phenotypes. Maybe your threshold to discard individuals has lead to remove from marker information all the individuals that are present in the phenotypic dataset.")
      }
      Mtrait <- Markers[which(rownames(Markers) %in% unique(mydataSub$designation)),]
    }
    
    ## next step
    if(nrow(mydataSub) == 0){
      print(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all envs."))
      traitToRemove <- c(traitToRemove,iTrait)
    }else{
      
      mydataSub <- mydataSub[which(mydataSub$designation != ""),] # remove blank designations
      
      ##############
      # do analysis
      if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
          Ve <- var(mydataSub[,"predictedValue"], na.rm = TRUE)
          # make sure the terms to be fitted have more than one level
          if(deregress){
            mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
          }
          
          if( modelTypeTrait[iTrait] == "blup" ){ # if user doesn't have any special model
            # make sure the matrix only uses the leves for individuals with data
            designationFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"])
            A <- diag(length(designationFlevels))
            colnames(A) <- rownames(A) <- designationFlevels
          }else{
            if(modelTypeTrait[iTrait] == "rrblup"){ # if user will use the rrblup
              # groupTrait <- LGrp
              # genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=designationFlevels, withMarkOnly=onlyInA)
            }else if(modelTypeTrait[iTrait] %in% "gblup"){ # if user wants to do a gblup, pblup or ssgblup model
              designationFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"]))
              if(modelTypeTrait[iTrait] %in% c("pblup","ssgblup")){ # we need to calculate NRM
                paramsPed <- phenoDTfile$metadata$pedigree
                if(length(intersect(paramsPed$value, colnames(phenoDTfile$data$pedigree)))  < 3){
                  stop("Metadata for pedigree (mapping) and pedigree frame do not match. Please reupload and map your pedigree information.", call. = FALSE)
                }
                N <- cgiarBase::nrm2(pedData=phenoDTfile$data$pedigree,
                                     indivCol = paramsPed[paramsPed$parameter=="designation","value"],
                                     damCol = paramsPed[paramsPed$parameter=="mother","value"],
                                     sireCol = paramsPed[paramsPed$parameter=="father","value"]
                )
              }
              if(modelTypeTrait[iTrait] %in% c("gblup","ssgblup")){ # we need to calculate GRM
                commonBetweenMandP <- intersect(rownames(Markers),designationFlevels)
                withoutMarkers <- setdiff(designationFlevels, rownames(Markers))
                if(length(commonBetweenMandP) < 2){ 
                  commonBetweenMandPInOriginal <- intersect(rownames(phenoDTfile$data$geno),designationFlevels)
                  if(length(commonBetweenMandPInOriginal) > 2){
                    warning(paste("Markers could not be matched with phenotypes for trait", iTrait," because marker QA filtering removed the markers for those individuals. A classical BLUP model will be fitted for this trait."), call. = FALSE)   
                    failedMarkerModel=TRUE
                  }else{
                    warning(paste("Markers could not be matched with phenotypes for trait", iTrait,". Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes. A classical BLUP model will be fitted for this trait."), call. = FALSE)              
                  }
                }
                # print(failedMarkerModel)
                if(failedMarkerModel){
                  modelTypeTrait[iTrait] <- "blup"
                  toFixFailure <- na.omit(unique(mydataSub[,"designation"]))
                  A <- diag(length(toFixFailure)); colnames(A) <- rownames(A) <- toFixFailure
                }else{
                  M <- Markers[commonBetweenMandP,]
                  if(ncol(M) > nMarkersRRBLUP){ # we remove that many markers if a big snp chip
                    A <- sommer::A.mat(M[,sample(1:ncol(M), nMarkersRRBLUP)])
                  }else{ A <- sommer::A.mat(M) };  M <- NULL
                  if(modelTypeTrait[iTrait] == "ssgblup"){ # only if ssgblup we merge
                    A <- sommer::H.mat(N,A, tau=1,  omega=1, tolparinv=1e-6)
                  }
                  A <- A + diag(1e-4, ncol(A), ncol(A))
                  if(length(withoutMarkers) > 0){
                    newNames <- c(colnames(A), withoutMarkers)
                    Ax <- diag(length(withoutMarkers))*mean(diag(A)); 
                    A <- bdiag(A, Ax)
                    colnames(A) <- rownames(A) <- newNames
                  }
                }
                
              }
              if(modelTypeTrait[iTrait] %in% c("pblup")){ A <- N  }
              if(!failedMarkerModel){
                badGeno <- which(rownames(A) == "") # should have no covariance with anyone
                if(length(badGeno) > 0){A[badGeno,2:ncol(A)]=0; A[2:nrow(A),badGeno]=0} # make zero covariance with this genotype
                badBlankGenotype <- which(colnames(A)=="")
                if(length(badBlankGenotype) > 0){A <- A[-badBlankGenotype,-badBlankGenotype]}
              }else{
                
              }
            }else{ # blue model
              
            }
          }
          
          if(useWeights){
            mydataSub$w  <- 1/(mydataSub$stdError^2) # add weights column
          }else{
            mydataSub$w  <- 1#/(mydataSub$stdError^2) # add weights column
          }
          
          # save.image(file="strangeBug.RData")
          relmat <- list(designation=A)
          form <<- as.formula( paste("predictedValue", "~", ff$form[[iTrait]])  )
          for(iVcov in 1:length(inputFormulation)){ # iVcov=1
            if("designation" %in% inputFormulation[[iVcov]]$right){
              toAdd <- setdiff(inputFormulation[[iVcov]]$right, "designation")
              if(length(toAdd) > 0){
                for(iTerm in toAdd){ # iTerm = toAdd[1]
                  lvs <- unique(mydataSub[, iTerm])
                  e <- which( inputFormulation[[iVcov]]$right == iTerm)
                  d <- which( inputFormulation[[iVcov]]$right == "designation")
                  E <- diag(length(lvs)); colnames(E) <- rownames(E) <- lvs
                  if(e > d){
                    A <- kronecker(A,E, make.dimnames = TRUE)
                  }else{A <- kronecker(E,A, make.dimnames = TRUE)}
                }
              }
              relmat <- list(A)
              names(relmat) <- paste(inputFormulation[[iVcov]]$right, collapse = ":")
            }
          }
          mydataSub <<-mydataSub
          # A <<- A
          relmat <<- relmat
          mix <- try(
            lmebreed(formula=form, 
                     family = NULL,  #REML = TRUE,
                     # addmat=list(),
                     start = NULL, verbose = TRUE,
                     weights = mydataSub$w,
                     # subset, na.action, offset,
                     contrasts = NULL, dateWarning=TRUE, returnParams=FALSE,
                     rotation=FALSE, coefOutRotation=8,
                     relmat=relmat,#list(designation=A),
                     control = lmerControl(
                       check.nobs.vs.nlev = "ignore",
                       check.nobs.vs.rankZ = "ignore",
                       check.nobs.vs.nRE="ignore"
                     ),
                     data = mydataSub),
            silent = TRUE
          )
          # print(mix)
          # summary(mix)
          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
            
            ## save the model formula used
            currentModeling <- data.frame(module="mtaFlex", analysisId=mtaAnalysisId,trait=iTrait, environment="general",
                                          parameter=c("formula","family","designationEffectType"), 
                                          value=c(as.character(form)[3],traitFamily[iTrait], toupper(modelTypeTrait[iTrait]) ))
            phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
            ## save the environments used goodFields
            currentModeling <- data.frame(module="mtaFlex", analysisId=mtaAnalysisId,trait=iTrait, environment=allEnvironments,
                                          parameter="includedInMta", 
                                          value=ifelse(allEnvironments%in%goodFields, TRUE, FALSE))
            phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
            
            blues <- as.data.frame(summary(mix)$coefficients)
            blues$Estimate <- lme4breeding::adjBeta(blues$Estimate)
            colnames(blues) <- c("predictedValue","stdError","tValue")
            blues$designation <- rownames(blues); blues$environment <- "(Intercept)"
            blues$reliability <- NA
            blues$entryType <- "fixedEffect"
            
            effs <- lme4breeding::ranef(object=mix, condVar=TRUE)
            intercept <- lme4::fixef(mix)[1]
            ## get variance components and fix names
            vc <- lme4::VarCorr(mix); #print(vc,comp=c("Variance"))
            vars <- as.data.frame(vc)#unlist(lapply(vc, function(x){diag(x)}))
            for(iName in names(effs)){
              possibleNames <- paste0(iName,".",1:50)
              for(ipn in possibleNames){ # ipn = possibleNames[1]
                vars$grp[which(vars$grp == ipn)] <- iName
              }
            }
            ## get predictions and metrics
            pp <- means <- sds <- cvs <- r2s <- list()
            
            for(iEffect in names(effs)){ # iEffect = names(effs)[1]
              Vg <- vars[which((vars$grp == iEffect) & (is.na(vars$var2) )), "vcov"]
              provEffects <- as.data.frame(effs[[iEffect]]); 
              
              means[[iEffect]] <- apply(provEffects+intercept,2,function(x){mean(x,na.rm=TRUE)})
              sds[[iEffect]] <- apply(provEffects,2,function(x){sd(x,na.rm=TRUE)})
              cvs[[iEffect]] <- (sds[[iEffect]] / means[[iEffect]]) * 100
              
              provEffects$designation <- rownames(provEffects)
              provEffectsLong <- reshape(provEffects, idvar = "designation", varying = list(1:(ncol(provEffects)-1)),
                                         v.names = "predictedValue", times=colnames(effs[[iEffect]]), timevar = "environment", direction = "long")
              provEffectsLong$predictedValue <- provEffectsLong$predictedValue + intercept
              SEs <- attr(effs[[iEffect]], which="postVar")
              if(is.list(SEs)){
                
                provSe <- do.call(cbind,
                                  lapply(SEs, function(x){
                                    se <- list()
                                    for(k in 1:dim(x)[1]){
                                      se[[k]] <- x[k,k,] # get the diagonal value of each column #apply(x[,,],3,function(y){y[k,k]})
                                    }
                                    return( do.call(cbind, se) )
                                  })
                )
                # provSe <- do.call(cbind, lapply(SEs,function(x){x[,,]}))
                r2 <- provSe
                for(iCol in 1:ncol(provSe)){r2[,iCol] <- (Vg[iCol] - provSe[,iCol])/Vg[iCol]}
                provEffectsLong$stdError <- sqrt( as.vector((provSe)) )
                # provEffectsLong$stdError <- unlist( lapply(SEs,function(x){x[,,]}) )
              }else{
                provSe <- as.matrix(apply(SEs,3,function(x){diag(x)}))
                if(ncol(provSe) > 1){ provSe <- t(provSe) } # only if a complex structure was fitted
                r2 <- provSe
                for(iCol in 1:ncol(provSe)){r2[,iCol] <- (Vg[iCol] - provSe[,iCol])/Vg[iCol]}
                provEffectsLong$stdError <- sqrt( as.vector((provSe)) )
                # provEffectsLong$stdError <- as.vector(t(apply(SEs,3,function(x){diag(x)})))
              }
              r2[which(r2 < 0, arr.ind = TRUE)]=0
              r2s[[iEffect]] <- apply(r2,2,function(x){mean(x, na.rm=TRUE)})
              provEffectsLong$reliability <- as.vector((r2))
              # check if there is an across estimate, otherwise create it
              if("(Intercept)" %!in% unique(provEffectsLong$environment) ){
                ppa <- aggregate(cbind(predictedValue,stdError,reliability) ~ designation, FUN=mean, data=provEffectsLong[which(provEffectsLong$predictedValue != 0),])
                ppa$environment <- "(Intercept)"
                provEffectsLong <- rbind(provEffectsLong,ppa[,colnames(provEffectsLong)])
              }
              provEffectsLong$entryType <- iEffect
              # print(head(provEffectsLong))
              pp[[iEffect]] <- provEffectsLong
            }
            pp <- do.call(rbind,pp)
            pp <- rbind(pp, blues[,colnames(pp)])
            pp$trait <- iTrait # add trait
            
            phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait, environment= vars[which( (is.na(vars$var2) )), "var1"],
                                                    parameter=paste0("V.", vars[which( (is.na(vars$var2) )), "grp"]), method=c("REML"),
                                                    value=vars[which( (is.na(vars$var2) )), "vcov"] ,
                                                    stdError=0
                                         ),
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait, environment= vars[which( (is.na(vars$var2) )), "var1"],
                                                    parameter=paste0("mean.", vars[which( (is.na(vars$var2) )), "grp"]), method=c("sum/n"),
                                                    value=c(unlist(means),NA ),
                                                    stdError=0
                                         ),
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait, environment= vars[which( (is.na(vars$var2) )), "var1"],
                                                    parameter=paste0("CV.", vars[which( (is.na(vars$var2) )), "grp"]), method=c("sum/n"),
                                                    value=c(unlist(cvs),NA ),
                                                    stdError=0
                                         ),
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait, environment= vars[which( (is.na(vars$var2) )), "var1"],
                                                    parameter=paste0("r2.", vars[which( (is.na(vars$var2) )), "grp"]), method=c("sum/n"),
                                                    value=c(unlist(r2s),NA ),
                                                    stdError=0
                                         ),
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait, environment= "general",
                                                    parameter="nEnv", method=c("n"), value= length(goodFields),stdError=0
                                         )
            )
            
          }else{ # if model failed
            if(verbose){ cat(paste("Mixed model failed for this combination. Aggregating and assuming h2 = 0 \n"))}
            pp <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
            pp$stdError <- sd(pp$predictedValue)  # 1
            # pp$status <- "Aggregated"
            pp$reliability <- 1e-6
            pp$trait <- iTrait
            pp$entryType <- ""
            cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
            ## save metrics
            phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                         data.frame(module="mtaFlex",analysisId=mtaAnalysisId, trait=iTrait,
                                                    environment="(Intercept)",
                                                    parameter=c("mean","CV", "r2","Vg","Vr","nEnv"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML","REML","n"),
                                                    value=c(mean(pp$predictedValue, na.rm=TRUE), cv, 0, 0, Ve, length(goodFields) ),
                                                    stdError=c(0,0,0,0, 0,0 )
                                         )
            )
            currentModeling <- data.frame(module="mtaFlex", analysisId=mtaAnalysisId,trait=iTrait, environment="general",
                                          parameter=c("formula","family","designationEffectType"), 
                                          value=c( "None","None","mean" ))
            phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
            
          }
          #######################################
          #######################################
          #######################################
          # Trait run is finished, add entryType
          mydataForEntryType <- droplevels(mydata[which(mydata$trait == iTrait),])
          pp$entryType <- paste( pp$entryType,
                                 apply(data.frame(pp$designation),1,function(x){
                                   found <- which(mydataForEntryType$designation %in% x)
                                   if(length(found) > 0){
                                     x2 <- paste(sort(unique(toupper(trimws(mydataForEntryType[found,"entryType"])))), collapse = "#");
                                   }else{x2 <- "unknownType"}
                                   return(x2)
                                 }), sep = "_")
          mydataForEntryType <- NULL
          if(!inherits(mix,"try-error") ){  # if model run OK
            if(modelTypeTrait[iTrait] == "rrblup"){ # rrblup model
              toUseE <- grep("designation",pp$entryType)
              pp$entryType[toUseE] <- paste( pp$entryType[toUseE],
                                     "GEBV",
                                     ifelse(as.character(pp$designation[toUseE]) %in% rownames(Mtrait) ,"tested", "predicted"),
                                     sep="_")
            }else{ # other model
              toUseE <- grep("designation",pp$entryType)
              pp$entryType[toUseE] <- paste(pp$entryType[toUseE],
                                    ifelse(as.character(pp$designation[toUseE]) %in% setdiff( unique(mydataSub$designation), colnames(A) ), "TGV", surrogate[modelTypeTrait[iTrait]]),
                                    ifelse(as.character(pp$designation[toUseE]) %in% setdiff(colnames(A), unique(mydataSub$designation) ), "predicted", "tested"), # 
                                    sep="_")
            }
          }else{ # if we just averaged
            pp$entryType <- paste(pp$entryType,"average","tested", sep="_")
          }
          ### save predictions
          predictionsList[[counter2]] <- pp;
          counter=counter+1
          failedMarkerModel=FALSE # reset if the trait model failed and was set to TRUE
        }
      }
    }
    counter2 = counter2+1
  }
  #
  trait <- setdiff(trait,traitToRemove)
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId
  ##########################################
  ## add stage and designation code
  entries <- unique(mydata[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(mydata[which(mydata$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4)
    return(y)
  }))
  predictionsBind <- merge(predictionsBind,baseOrigin, by="designation", all.x=TRUE)
  predictionsBind$module <- "mtaFlex"
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="mtaFlex", analysisId=mtaAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="mtaFlex",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  rownames(phenoDTfile$metrics) <- NULL
  rownames(phenoDTfile$modeling) <- NULL
  return(phenoDTfile)
}
