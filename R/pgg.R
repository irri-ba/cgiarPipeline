pgg <- function(
    phenoDTfile= NULL,
    analysisId=NULL, # sta ID
    trait=NULL, # per trait
    by=NULL,
    percentage=10,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES THE PREDICTED GENETIC GAIN FOR TRAITS
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  pggAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the ID of the analysis to use as input", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(by)){by <- "environment"}
  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions 
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  paramsPheno <- phenoDTfile$metadata$pheno
  paramsPheno <- paramsPheno[which(paramsPheno$parameter != "trait"),]
  colnames(mydata) <- cgiarBase::replaceValues(colnames(mydata), Search = paramsPheno$value, Replace = paramsPheno$parameter )
  
  # add missing columns
  '%!in%' <- function(x,y)!('%in%'(x,y))
  keep <- which( paramsPheno$parameter %!in% c("trait","designation","environment","rep","row","col","iBlock") )
  if(length(keep) > 0){
    toExtractFromData <- paramsPheno[keep, "value"]
    tpe <- unique(phenoDTfile$data$pheno[,c("environment",toExtractFromData)])
    colnames(tpe) <- cgiarBase::replaceValues(colnames(tpe), Search = paramsPheno$value, Replace = paramsPheno$parameter )
    mydata <- merge(mydata, tpe, by="environment", all.x = TRUE)
  }
  
  myPed <- phenoDTfile$data$pedigree
  paramsPed <- phenoDTfile$metadata$pedigree
  colnames(myPed) <- cgiarBase::replaceValues(colnames(myPed), Search = paramsPed$value, Replace = paramsPed$parameter )
  
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  # if(is.null(environment)){environment <- na.omit(unique(mydata$environment))}
  if(is.null(myPed) || (nrow(myPed) == 0 ) ){stop("yearOfOrigin column was not matched in your original file. Please correct.", call. = FALSE)}
  yearsToUse <- as.character(unique(myPed$yearOfOrigin))
  mydata <- merge(mydata, myPed[,c("designation","yearOfOrigin")], by="designation", all.x=TRUE )
  mydata <- mydata[which(!is.na(mydata$yearOfOrigin)),]
  
  # myPheno <- phenoDTfile$data$pheno
  # colnames(myPheno) <- cgiarBase::replaceValues(colnames(myPheno), Search = paramsPheno$value, Replace = paramsPheno$parameter )
  
  if(length(which(paramsPheno$parameter == "year")) == 0){
    mydata$year <- mydata$yearOfOrigin 
    cat("Year column was not mapped, assuming year of origin and first year of testing is the same")
  }
  ############################
  ## gg analysis
  p <- percentage/100
  i <- dnorm(qnorm(1 - p))/p
  # counter=1
  for(iTrait in trait){ # iTrait = trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    uEnvironments <- unique(mydata[,by])
    for(uE in uEnvironments){ # uE <- uEnvironments[1]
      # subset data
      mydataSub <- droplevels(mydata[which((mydata$trait == iTrait) & (mydata[,by] %in% uE)  ),])
      # calculate parameters
      rels <- mydataSub$rel;
      badrels <- which(rels < 0); if(length(badrels) > 0){ rels[badrels] <- 1e-6}
      r <- ifelse(length(na.omit(rels)) > 0, mean(sqrt(na.omit(rels)), na.rm=TRUE), 1e-6)
      sigma<- sd(mydataSub$predictedValue, na.rm = TRUE)
      mu<- mean(mydataSub$predictedValue, na.rm = TRUE)
      min.x<- min(mydataSub$predictedValue, na.rm = TRUE)
      max.x<- max(mydataSub$predictedValue, na.rm = TRUE)
      mydataSubSorted <- mydataSub[with(mydataSub, order(-predictedValue)), ]
      mydataSubSortedSel <- mydataSubSorted[1:round(nrow(mydataSubSorted) * p),]
      age <- mean(mydataSubSortedSel$year, na.rm=TRUE) - mean(mydataSubSortedSel$yearOfOrigin, na.rm=TRUE)
      R <- r * sigma * i
      ggAge =  R/age
      nTrials = length(unique(mydataSub[,"environment"]))
      ##
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                   data.frame(module="pgg",analysisId=pggAnalysisId, trait= iTrait, environment=uE, 
                                              parameter=c("r","sigmaG","meanG","min.G","max.G", "cycleLength","i","R","PGG","nEnvs"), 
                                              method=c("sqrt(r2)","sd(BLUP)","sum(x)/n","min(x)","max(x)","yearTest-yearOrigin","dnorm(qnorm(1 - p))/p","r*sigma*i","R/cycleLength","sum"), 
                                              value=c(r,sigma, mu, min.x, max.x, age, i, R, ggAge,nTrials), 
                                              stdError=0
                                   )
      )
      currentModeling <- data.frame(module="pgg", analysisId=pggAnalysisId,trait=iTrait, environment=uE, 
                                    parameter=c("percentage(%)","verbose"), value=c(percentage, verbose))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
    }
    
  }
  #########################################
  # update databases
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="pgg", analysisId=pggAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="pgg",  analysisId=pggAnalysisId, trait=c("inputObject"), environment="general", 
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)#paste("pgg done:",id))
}
