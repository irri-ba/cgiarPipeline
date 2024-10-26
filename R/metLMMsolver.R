metLMMsolver <- function(
    phenoDTfile= NULL, analysisId=NULL,
    fixedTerm= list("1"),  randomTerm=NULL, expCovariates=NULL,
    envsToInclude=NULL, trait= NULL, traitFamily=NULL, useWeights=TRUE,
    calculateSE=TRUE, heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf, nPC=NULL,   # subsetVariable=NULL, subsetVariableLevels=NULL,
    maxIters=50,  verbose=TRUE
){
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING LMM SOLVER
  mtaAnalysisId <- as.numeric(Sys.time())
  namesSeq <- function(x){
    nCharX <- nchar(x)
    maxZeros <- max(nCharX) 
    nZeros <- abs(nCharX - maxZeros)
    zeros <- apply(data.frame(nZeros),1,function(x){paste(rep("0",x), collapse = "")})
    res <- paste0(zeros, as.character(x))
    return(res)
  }
  ##########################################
  ##########################################
  ## CONTROLS FOR MISSPECIFICATION (6 lines)
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the STA analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(!is.null(randomTerm)){
    randomTerm <- unique(randomTerm)
    if(is.null(expCovariates)){expCovariates <- randomTerm; expCovariates <- lapply(expCovariates, function(x){rep("none",length(x))})}else{
      if(length(expCovariates) != length(randomTerm)){
        stop("Please ensure that expCovariates and randomTerm arguments have the same length.", call. = FALSE)
      }else{
        if( sum(unlist(mapply('-', lapply(randomTerm,length), lapply(expCovariates,length), SIMPLIFY = FALSE))) != 0){
          stop("Please ensure that expCovariates and randomTerm arguments have the same length.", call. = FALSE)
        }
      } 
    }
  }else{if(length(randomTerm) == 0){randomTerm <- NULL}}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- "1"}else{fixedTerm <- unique(fixedTerm)}
  traitsForExpCovariates <- unique(phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId),"trait"])
  if(!is.null(nPC)){if(length(intersect(names(nPC), c("geno","weather","pedigree", traitsForExpCovariates))) == 0){stop("The nPC argument needs to be a named numeric vector with names 'geno', 'weather' and 'pedigree' or traits available.", call. = FALSE)}}
  ##########################################
  ##########################################
  ## EXTRACT POSSIBLE EXPLANATORY COVARIATES AND FORM KERNELS (30 lines)
  Weather <- cgiarPipeline::summaryWeather(phenoDTfile, wide=TRUE) # in form of covariates
  Weather <- apply(Weather,2,sommer::imputev)
  colnames(Weather) <- gsub(" ","",colnames(Weather))
  covars <- unique(unlist(expCovariates))
  randomTermForCovars <- unique(unlist(randomTerm))
  if(!is.null(randomTermForCovars)){
    if( any( covars %in% c("geno", "weather","pedigree", traitsForExpCovariates ) ) ){
      if(verbose){message("Checking and calculating kernels requested")}
      ## MARKER KERNEL
      Markers <- phenoDTfile$data$geno # in form of covariates  
      if("geno" %in% covars & !is.null(Markers)){
        classify <- randomTermForCovars[which(covars %in% "geno")]
        if(verbose){message(paste("   Marker kernel for",paste(classify,collapse = ","), "requested"))}
        qas <- which( phenoDTfile$status$module == "qaGeno" ); qas <- qas[length(qas)]
        if(length(qas) > 0){
          modificationsMarkers <- phenoDTfile$modifications$geno[which(phenoDTfile$modifications$geno$analysisId %in% qas),]
          Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
          if(length(which(is.na(Markers))) > 0){Markers <- apply(Markers,2,sommer::imputev)}
        }else{
          Markers <- apply(Markers,2,sommer::imputev)
        }
        G <- sommer::A.mat(Markers-1);  G <- G + diag(1e-5, ncol(G), ncol(G))
        Gchol <- t(chol(G))
        if(nPC["geno"] > 0){
          if(verbose){message("   Eigen decomposition of marker kernel requested")}
          decomp <- RSpectra::svds(Gchol, k = min(c(nPC["geno"], ncol(Gchol))), which = "LM")
          rownames(decomp$u) <- rownames(G); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
          Gchol <- decomp$u 
        }
      }
      ## WEATHER KERNEL
      if("weather" %in% covars & !is.null(Weather)){
        classify <- randomTermForCovars[which(covars %in% "weather")]
        if(verbose){message(paste("   Weather kernel for",paste(classify,collapse = ","), "requested"))}
        rownamesWeather <- rownames(Weather)
        Weather <- apply(Weather, 2, scale)
        Weather <- Weather[,which( !is.na(apply(Weather,2,var)) ), drop=FALSE]
        W <- sommer::A.mat(Weather)
        W <- W + diag(1e-5, ncol(W), ncol(W))
        rownames(W) <- colnames(W) <- rownamesWeather
        Wchol <- t(chol(W))
        if(nPC["weather"] > 0){
          if(verbose){message("   Eigen decomposition of weather kernel requested")}
          decomp <- RSpectra::svds(Wchol, k = min(c(nPC["weather"], ncol(Wchol))), which = "LM")
          rownames(decomp$u) <- rownames(W); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
          Wchol <- decomp$u 
        }
      }
      ## PEDIGREE KERNEL
      Pedigree <- phenoDTfile$data$pedigree
      if("pedigree" %in% covars  &  !is.null(Pedigree)){
        classify <- randomTermForCovars[which(covars %in% "pedigree")]
        if(verbose){message(paste("   Pedigree kernel for",paste(classify,collapse = ","), "requested"))}
        paramsPed <- phenoDTfile$metadata$pedigree
        N <- cgiarBase::nrm2(pedData=phenoDTfile$data$pedigree,
                             indivCol = paramsPed[paramsPed$parameter=="designation","value"],
                             damCol = paramsPed[paramsPed$parameter=="mother","value"],
                             sireCol = paramsPed[paramsPed$parameter=="father","value"]
        )
        Nchol <- t(chol(N))
        if(nPC["pedigree"] > 0){
          if(verbose){message("   Eigen decomposition of pedigree kernel requested")}
          decomp <- RSpectra::svds(Nchol, k = min(c(nPC["pedigree"], ncol(Nchol))), which = "LM")
          rownames(decomp$u) <- rownames(N); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
          Nchol <- decomp$u 
        }
      } # now is in the form of covariates
      # TRAIT-BASED KERNEL (ALWAYS ROW-GROUPED BY DESIGNATION)
      if(any(covars %in% traitsForExpCovariates) ){
        TraitKernels <- list()
        covarsTraits <- intersect(covars,traitsForExpCovariates)
        for(iCovar in covarsTraits){ # iCovar = covarsTraits[1] # for each trait specified in covar
          classify <- randomTermForCovars[which(covars == iCovar)] # identify at what levels should the trait be classified
          if(verbose){message(paste("  ",iCovar,"kernel for",classify, "requested"))}
          baseData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId & (phenoDTfile$predictions$trait == iCovar) ),]
          wideTrait <- reshape(baseData[,c(classify,"designation","predictedValue")], direction = "wide", 
                               idvar = "designation", timevar = classify, v.names = "predictedValue", sep= "_")
          wideTrait <- apply(wideTrait[,-1],2,sommer::imputev)
          colnames(wideTrait) <- gsub("predictedValue_","",colnames(wideTrait))
          S <- cov(wideTrait)
          S <- as.matrix(Matrix::nearPD(x = S, corr = FALSE, 
                                        keepDiag = FALSE, base.matrix = FALSE, do2eigen = TRUE, 
                                        doSym = FALSE, doDykstra = TRUE, only.values = FALSE, 
                                        ensureSymmetry = !isSymmetric(Sigma), eig.tol = 1e-06, 
                                        conv.tol = 1e-07, posd.tol = 1e-08, maxit = 100, conv.norm.type = "I", 
                                        trace = FALSE)$mat)
          Schol <- t(chol(S))
          if(nPC[iCovar] > 0){
            if(verbose){message(paste("   Eigen decomposition of",iCovar," classified by",classify, "kernel requested"))}
            decomp <- RSpectra::svds(Schol, k = min(c(nPC[iCovar], ncol(Schol))), which = "LM")
            rownames(decomp$u) <- rownames(S); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
            Schol <- decomp$u 
          }
          TraitKernels[[iCovar]][[classify]] <- Schol
        }
      }
    }
  }
  ##########################################
  ##########################################
  ## COMPLETE THE CLEANING PARAMETERS (7 lines)
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait)); heritLB <- heritLB[1:length(trait)]; names(heritLB) <- trait
  heritUB <- rep(heritUB,length(trait)); heritUB <- heritUB[1:length(trait)]; names(heritUB) <- trait
  meanLB <- rep(meanLB,length(trait)); meanLB <- meanLB[1:length(trait)]; names(meanLB) <- trait
  meanUB <- rep(meanUB,length(trait)); meanUB <- meanUB[1:length(trait)]; names(meanUB) <- trait
  traitOrig <- trait ## ?????????
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- list("1")} # assign the intercept if there's no fixed effects
  ##########################################
  ##########################################
  # LOAD THE DATASET AND EXTEND IT TO INCLUDE METADATA (16 lines)
  if(verbose){message("Loading the dataset and adding metadata.")}
  mydata <- phenoDTfile$predictions #
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if (nrow(mydata) < 2) stop("Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.", call. = FALSE)
  metaPheno <- phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter %in% c("pipeline","stage","environment","year","season","timepoint","country","location","trial","study","management")),]
  otherMetaCols <- unique(phenoDTfile$data$pheno[,metaPheno$value,drop=FALSE])
  colnames(otherMetaCols) <- cgiarBase::replaceValues(Source = colnames(otherMetaCols), Search = metaPheno$value, Replace = metaPheno$parameter )
  otherMetaCols <- otherMetaCols[which(!duplicated(otherMetaCols[,"environment"])),,drop=FALSE] # we do this in case the users didn't define the environment properly
  mydata <- merge(mydata, otherMetaCols, by="environment", all.x = TRUE)
  WeatherRow <- as.data.frame(Weather); WeatherRow$environment <- rownames(WeatherRow)
  mydata <- merge(mydata, WeatherRow, by="environment", all.x = TRUE)
  # if(!is.null(subsetVariable) & !is.null(subsetVariableLevels)){
  #   if(subsetVariable %in% colnames(otherMetaCols)){
  #     forSubset <- which(mydata[,subsetVariable] %in% subsetVariableLevels)
  #     if(length(forSubset) > 0){mydata <- mydata[forSubset,]}
  #   }
  # }
  ##########################################
  ##########################################
  # CHECK THE ENVS TO INCLUDE PER TRAIT (6 lines)
  if(is.null(envsToInclude)){
    envsToInclude=  as.data.frame( do.call( rbind, list (with(mydata, table(environment,trait)) ) ) )
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if(nrow(bad) > 0){envsToInclude[bad] = 0}
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }; allEnvironments <- rownames(envsToInclude)
  ##########################################
  ##########################################
  # BUILD THE DATASETS FOR MODEL FITTING (100 lines)
  if(verbose){message("Building trait datasets.")}
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId),]
  myDataTraits <- fixedTermTrait <- randomTermTrait <- groupingTermTrait <- Mtrait <- envsTrait <- list()
  for(iTrait in trait){ # iTrait = trait[1]
    # filter for records available
    vt <- which(mydata[,"trait"] == iTrait)
    if(length(vt) > 0){ # we have data for the trait
      prov <- mydata[vt,]
      # filter by the environments to include
      vte <- which(mydata[,"environment"] %in% rownames(envsToInclude)[as.logical(envsToInclude[,iTrait])])
      prov <- prov[vte,]
      # remove bad environment based on h2 and r2
      pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2","r2", apply(expand.grid( c("plotH2","H2","meanR2","r2"), c("designation","mother","father")),1,function(f){paste(f,collapse = "_")}) )),]
      goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value >= heritLB[iTrait]) & (pipeline_metricsSub$value <= heritUB[iTrait])),"environment"])
      prov <- prov[which(prov$environment %in% goodFields),]
      # remove bad environment based on environment means
      pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("mean")),]
      goodFieldsMean <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > meanLB[iTrait]) & (pipeline_metricsSub$value < meanUB[iTrait])),"environment"])
      prov <- prov[which(prov$environment %in% goodFields),]
      if(nrow(prov) > 0){ # if after filters there's still data for this trait we can continue and save the data
        if( var(prov[,"predictedValue"], na.rm = TRUE) > 0 ){ # check that there is variance
          # make new formula for this specific trait if it passed all the filters
          fixedTermProv <- fixedTerm
          
          for(iFixed in 1:length(fixedTermProv)){ # for each element in the list # iFixed=1
            fixedTermProv2 <- fixedTermProv[[iFixed]]
            for(iFixed2 in  fixedTermProv2){ # for each factor in the interactions # iFixed2 = fixedTermProv2[1]
              if(iFixed2 != "1"){
                if( length( table(prov[,iFixed2]) ) == 1 ){ fixedTermProv[[iFixed]] <- setdiff( fixedTermProv[[iFixed]], iFixed2 )}
              }
            }
          }
          fixedTermTrait[[iTrait]] <- unique(fixedTermProv)
          # random formula per trait
          randomTermProv <- randomTerm
          if(!is.null(randomTermProv)){
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=2
              randomTermProv2 <- randomTermProv[[irandom]]
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 2
                if( length( table(prov[,randomTermProv2[irandom2]]) ) == 1 ){ randomTermProv[[irandom]] <- setdiff( randomTermProv[[irandom]], randomTermProv2[irandom2] )}
              }
            }
          }
          # any term that is modified from what user specified we remove it totally, is better than fitting something undesired
          goodTerms <- which( ( unlist(lapply(randomTerm,length)) - unlist(lapply(randomTermProv,length)) ) == 0 )
          randomTermProv <- randomTerm[goodTerms]
          randomTermProv <- unique(randomTermProv)
          # if reduced models reduce the datasets to the needed explanatory covariates
          if(!is.null(randomTermProv)){
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=1
              randomTermProv2 <- randomTermProv[[irandom]]
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 1
                if( expCovariates[[irandom]][irandom2] == "weather"){
                  M <- Wchol
                }else if(expCovariates[[irandom]][irandom2] == "geno"){
                  M = Gchol
                }else if(expCovariates[[irandom]][irandom2] == "pedigree"){
                  M <- Nchol
                }else if(expCovariates[[irandom]][irandom2] %in% traitsForExpCovariates){ # Trait kernel
                  classify <- randomTermForCovars[which(covars == expCovariates[[irandom]][irandom2])]
                  M <- TraitKernels[[expCovariates[[irandom]][irandom2]]][[classify]] # Schol equivalent
                }else{ # No kernel
                  namesZ <- unique(prov[,randomTermProv2[irandom2]])
                  M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                }
                goodLevels <- intersect(unique(prov[,randomTermProv2[irandom2]]), rownames(M) )
                if(length(goodLevels) > 0){ # only if we make a match we reduce the dataset
                  prov <- prov[which(prov[,randomTermProv2[irandom2]] %in% goodLevels),]
                }
              }
            }
          }
          ## build and add the incidence matrices
          groupingTermProv <- Mprov <- envsProv <- list()
          if(!is.null(randomTermProv)){
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=4
              randomTermProv2 <- randomTermProv[[irandom]]
              xxList <- Mlist <- list()
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 1
                if( expCovariates[[irandom]][irandom2] == "weather"){
                  M <- Wchol # Weather
                }else if(expCovariates[[irandom]][irandom2] == "geno"){
                  M = Gchol # Markers
                }else if(expCovariates[[irandom]][irandom2] == "pedigree"){
                  M <- Nchol # Pedigree
                }else if(expCovariates[[irandom]][irandom2] %in% traitsForExpCovariates){ # Trait kernel
                  classify <- randomTermForCovars[which(covars == expCovariates[[irandom]][irandom2])]
                  M <- TraitKernels[[expCovariates[[irandom]][irandom2]]][[classify]] # Schol equivalent
                }else{ # No kernel
                  if( unlist(lapply(prov, class))[randomTermProv2[irandom2]] %in% c("factor","character") ){
                    namesZ <- unique(prov[,randomTermProv2[irandom2]])
                    M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                  } else{ # numeric or integer
                    M <- Matrix::Diagonal(n=1); rownames(M) <- colnames(M) <- randomTermProv2[irandom2]
                  }
                }
                if( unlist(lapply(prov, class))[randomTermProv2[irandom2]] %in% c("factor","character") ){
                  goodLevels <- intersect(unique(prov[,randomTermProv2[irandom2]]), rownames(M) )
                  if(length(goodLevels) == 0){ # if no match then use the regular model matrix
                    namesZ <- unique(prov[,randomTermProv2[irandom2]])
                    M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                  }
                  xx = lme4breeding::redmm(x=prov[,randomTermProv2[irandom2]], M=M, nPC=0) 
                }else{xx <- sommer::isc(prov[,randomTermProv2[irandom2]])$Z}
                # if(nPCuse > 0){colnames(xx) <- paste(randomTermProv2[irandom2], colnames(xx), sep=":")}
                xxList[[irandom2]] = xx
                Mlist[[irandom2]] = M
                if(irandom2 > 1){
                  if(ncol(xxList[[irandom2-1]]) > 1){
                    m1 <- sommer::dsc(xxList[[irandom2-1]])
                  }else{m1 <- sommer::isc(xxList[[irandom2-1]][,1]) }
                  m2 <- sommer::isc(xx)
                  m3 <- sommer::vsc( m1  , m2  )
                  environmentCol <- list()
                  for(o in 1:length(m3$Z)){environmentCol[[o]] <- rep(colnames(m3$theta)[o],nrow(M))}
                  ff <- do.call( "cbind", m3$Z )
                  # envs <- unlist(environmentCol)
                  M <- kronecker(Mlist[[irandom2-1]] , M, make.dimnames = TRUE)
                }else{
                  ff <- xxList[[irandom2]]
                  # envs <- rep("(Intercept)",nrow(M))
                }
              }
              namesForEnvs <- lapply(Mlist,function(x){rownames(x)})
              namesForEnvs=do.call(expand.grid, rev(namesForEnvs))
              if(ncol(namesForEnvs)==1){ # if there's no interactions
                envs <- rep("(Intercept)",nrow(M))
              }else{ # if there's interactions
                nLevsInEnvs <- apply(namesForEnvs,2, function(x){length(unique(x))})
                # remove the one with the biggest number of levels
                namesForEnvs <- namesForEnvs[,-c(which(nLevsInEnvs == max(nLevsInEnvs))),drop=FALSE]
                envs <- apply(namesForEnvs,1,function(x){paste(x,collapse = ":")})
              }
              xxList=NULL;Mlist=NULL
              groupingTermProv[[irandom]] <- c( (ncol(prov)+1) : ( ncol(prov)+ncol(ff) ) ) # build grouping term
              prov <- cbind(prov, as.matrix(ff)) # bind matrix to dataset
              Mprov[[irandom]] <- M # save M matrix
              envsProv[[irandom]] <- envs # save levels for environment
            }
          }
          randomTermTrait[[iTrait]] <- unique(randomTermProv)
          myDataTraits[[iTrait]] <- prov # dataset
          names(groupingTermProv) <- names(envsProv) <- names(randomTermTrait[[iTrait]]) <- names(Mprov) <- unlist(lapply(randomTermProv, function(x){paste(x,collapse = "_")}))
          groupingTermTrait[[iTrait]] <- groupingTermProv
          Mtrait[[iTrait]] <- Mprov
          envsTrait[[iTrait]] <- envsProv
          # finish formula
        }
      }
    }
  }
  ##########################################
  ##########################################
  ## MODEL FITTING
  if(verbose){message("Fitting a model.")}
  predictionsList <- list(); 
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){message(paste("Analyzing trait", iTrait))}
    ## get data subset
    mydataSub <- myDataTraits[[iTrait]] 
    groupingSub <- groupingTermTrait[[iTrait]]
    Msub <- Mtrait[[iTrait]]
    envsSub <- envsTrait[[iTrait]]
    fixedTermSub <- fixedTermTrait[[iTrait]]
    randomTermSub <- randomTermTrait[[iTrait]]
    VarFull <- var(mydataSub[,"predictedValue"], na.rm = TRUE)
    ## deregress if needed
    effectTypeTrait <- phenoDTfile$modeling[which(phenoDTfile$modeling$analysisId == analysisId & phenoDTfile$modeling$trait == iTrait & phenoDTfile$modeling$parameter == "designationEffectType"),"value"]
    if(names(sort(table(effectTypeTrait), decreasing = TRUE))[1] == "BLUP"){ # if STA was BLUPs deregress
      mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
    }
    ## calculate weights
    mydataSub=mydataSub[with(mydataSub, order(environment)), ] # sort by environments
    mydataSub$w <- 1/(mydataSub$stdError^2) # add weights column
    ## get formula
    fix <- paste( unlist(lapply(fixedTermSub, function(x){paste(x, collapse = ":")})), collapse = " + ")
    fix <- paste("predictedValue ~", fix)
    if(length(groupingSub) > 0){
      ranran <- paste("~", paste(paste0("grp(",names(groupingSub),")"), collapse=" + "))
    }else{ranran <- character()}
    if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
    # warnin messages in weights use
    if(useWeights){
      weightsFormulation="w"
      if(verbose){message("   Using weights in the analysis. Residual variance will be fixed to 1.")  }
    }else{
      weightsFormulation=NULL
      if(verbose){message("   Ignoring weights in the analysis. Residual variance will be estimated.")  }
    }
    if(is.null(randomTermSub)){groupingSub=NULL}
    ## model fit
    mix <- try(
      LMMsolver::LMMsolve(fixed =as.formula(fix),
                          random = ranFormulation,
                          # residual=ranres,
                          weights = weightsFormulation,
                          # ginverse = myGinverse,
                          group = groupingSub,
                          family = eval(parse(text = traitFamily[iTrait])),
                          data = mydataSub, maxit = maxIters),
      silent = TRUE
    )
    # print(mix$VarDf)
    if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
      ## save the modeling used
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                    parameter=c("fixedFormula","randomFormula","family"), 
                                    value=c(fix,ifelse(length(ranran)>0,ranran,NA),traitFamily[iTrait] ))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      ## save the environments used goodFields
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment=allEnvironments,
                                    parameter="includedInMta", 
                                    value=ifelse(allEnvironments%in%goodFields, TRUE, FALSE))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      # get variance components
      ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
      Ve <- ss["residual","Variance"] 
      pp <- list()
      mu <- mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
      if(length(mu) > 0){
        pp[["(Intercept)"]] <- data.frame(designation="(Intercept)", predictedValue=mu, stdError=sqrt(as.matrix(solve(mix$C))[1,1]), reliability=NA,
                                        trait=iTrait, entryType="(Intercept)", environment="(Intercept)" )
      }
      fixedEffects <- setdiff(mix$EDdf$Term, mix$VarDf$VarComp)
      fixedEffects <- setdiff(fixedEffects, "(Intercept)")
      for(iGroupFixed in fixedEffects){ # iGroupFixed = fixedEffects[1]
        pick <- mix$ndxCoefficients[[iGroupFixed]]
        shouldBeOne <- which(pick == 0)
        if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
        blue <- mix$coefMME[pick] + mu; names(blue) <- names(pick); blue[1] <- blue[1]-mu
        start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroupFixed) - 1),"Model"]) # we don't add a one because we need the intercept
        nEffects <- length(blue)
        pev <- as.matrix(solve(mix$C))[start:(start+nEffects-1),start:(start+nEffects-1)]
        stdError <- (sqrt(Matrix::diag(pev)))
        pp[[iGroupFixed]] <- data.frame(designation=names(blue), predictedValue=blue, stdError=stdError, reliability=NA,
                                        trait=iTrait, entryType=iGroupFixed, environment="(Intercept)" )
      }; 
      if(!is.null(randomTermSub)){
        for( iGroup in names(groupingSub)){ # iGroup=names(groupingSub)[1]
          pick <- mix$ndxCoefficients[[iGroup]]
          shouldBeOne <- which(pick == 0)
          if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
          blup <- (Msub[[iGroup]] %*% mix$coefMME[pick]) + mu; blup <- as.vector(blup)
          names(blup) <- rownames(Msub[[iGroup]]) # Msub will always be a matrix wither a diagonal or different than but do it across for consistency
          start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
          nEffects <- ncol(Msub[[iGroup]])
          Vg <- ss[iGroup,"Variance"] 
          if(calculateSE){
            if(verbose){message(paste("   Calculating standar errors for",iTrait, iGroup,"predictions"))}
            Cinv <- solve(mix$C) 
            Cinv <- Cinv[start:(start+nEffects-1),start:(start+nEffects-1)]
            Cinv <- as(Cinv, Class = "dgCMatrix")
            startPev <- seq(1, length(blup), 500)
            endPev  <- c(startPev - 1, length(blup)); endPev <- endPev[-1]
            stdError <- list()
            for(s in 1:length(startPev)){
              use <- (startPev[s]:endPev[s])
              stdError[[s]] <-  sqrt(Matrix::diag( Msub[[iGroup]][use,] %*% Matrix::tcrossprod( Cinv, Msub[[iGroup]][use,]) ) ) 
            }
            stdError <- unlist(stdError)
            reliability <- abs((Vg - (stdError^2)) /Vg) # reliability <- abs((Vg - Matrix::diag(pev))/Vg)
          }else{stdError <- reliability <- rep(NA,length(blup))}
          badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
          badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
          pp[[iGroup]] <- data.frame(designation=names(blup), predictedValue=blup, stdError=stdError, reliability=reliability,
                                     trait=iTrait, entryType=iGroup, environment=envsSub[[iGroup]] )
          cv <- (sd(blup,na.rm=TRUE)/mean(blup,na.rm=TRUE))*100
          phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                       data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait= iTrait, 
                                                  environment=paste(unique(envsSub[[iGroup]]), collapse = "_"),
                                                  parameter=paste(c("mean","CV", "r2","Var"),iGroup,sep="_"), 
                                                  method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"),
                                                  value=c(mean(blup, na.rm=TRUE), cv, median(reliability), var(blup, na.rm=TRUE) ),
                                                  stdError=c(NA,NA,sd(reliability, na.rm = TRUE)/sqrt(length(reliability)),NA )
                                       )
          )
        }
      }
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                   data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait= iTrait, environment="across",
                                              parameter=c("Var_residual","nEnv"), method=c("REML","n"), value=c( Ve, length(goodFields) ),stdError=c(NA) )
      )
      pp <- do.call(rbind,pp)
    }else{ # if model failed
      if(verbose){ cat(paste("Mixed model failed for this combination. Aggregating and assuming h2 = 0 \n"))}
      pp <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
      pp$environment <- "(Intercept)"
      pp$stdError <- sd(pp$predictedValue)  
      pp$reliability <- 1e-6
      pp$trait <- iTrait
      cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
      ## save metrics
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                   data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait=iTrait,
                                              environment="across",
                                              parameter=c("mean","CV", "r2","Var_designation","Var_residual","nEnv"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML","REML","n"),
                                              value=c(mean(pp$predictedValue, na.rm=TRUE), cv, NA, NA, NA, length(goodFields) ),
                                              stdError=NA
                                   )
      )
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                    parameter=c("fixedFormula","randomFormula","family","designationEffectType"), 
                                    value=c("None","None","None","mean"))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
    }
    #######################################
    #######################################
    # add additional entry type labels
    if(!inherits(mix,"try-error") ){ 
      mydataForEntryType <- droplevels(mydata[which(mydata$trait == iTrait),])
      entryType <- apply(data.frame(pp$designation),1,function(x){
        found <- which(mydataForEntryType$designation %in% x)
        if(length(found) > 0){
          x2 <- paste(sort(unique(toupper(trimws(mydataForEntryType[found,"entryType"])))), collapse = "#");
        }else{x2 <- ""}
        return(x2)
      })
      pp$entryType <- ifelse(entryType != "", paste(pp$entryType, entryType, sep = "_"), pp$entryType)
      predictionsList[[iTrait]] <- pp;
    }
  }
  ## enf of model fitting
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId
  ##########################################
  ## add timePoint of origin, stage and designation code
  if(verbose){message("Wrapping the results.")}
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
  predictionsBind$module <- "mtaLmms"; rownames(predictionsBind) <- NULL
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="mtaLmms", analysisId=mtaAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="mtaLmms",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
