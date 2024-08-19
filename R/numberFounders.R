numberFounders <- function(
    object= NULL,
    analysisIdForGenoModifications=NULL,
    neExplore=NULL, 
    maxMarker=1000, 
    nSamples=5,
    verbose=FALSE
){
  
  neMarker <- function(M, neExplore=NULL, maxMarker=1000, nSamples=5, verbose=FALSE){
    # maxMarker argument: only used a limited number of markers to avoid this to be too time consuming
    v <- sample(1:ncol(M), min(c(maxMarker, ncol(M))))
    M <- M[,v]
    # calculate the total number of alleles in the population
    nAllelesPop <- apply(M,2,function(x){ifelse(length(table(x)) > 1, 2, 1)})
    nAllelesPopTotal <- sum(nAllelesPop)
    # maxNe argument: define the range to test
    if(is.null(neExplore)){neExplore <- seq(10,100,10)}
    # do the sampling algorithm
    counter <- 1
    allelesCovered <- allelesCoveredSe <- vector(mode="numeric", length = length(neExplore) )
    for(i in neExplore){ # for a possible Ne
      if(verbose){print(paste("Exploring allele coverage (%) at Ne:",i))}
      allelesCoveredSample <- vector(mode="numeric", length = nSamples)
      # nSamples argument: take a couple of samples 
      for(j in 1:nSamples){
        ii <- sample(1:nrow(M),i) # sample i individuals
        nAllelesPopI <- apply(M[ii,],2,function(x){ifelse(length(table(x)) > 1, 2, 1)}) # how many alleles we collect in the sample
        allelesCoveredSample[j] <- sum(nAllelesPopI) # sum them up
      }
      allelesCovered[counter] <- mean(allelesCoveredSample)/nAllelesPopTotal # mean across samples
      allelesCoveredSe[counter] <- ( sd(allelesCoveredSample/nAllelesPopTotal) ) # SE across samples 
      counter <- counter+1
    }
    # save results
    result <- data.frame(allelesCovered=allelesCovered, allelesCoveredSe=allelesCoveredSe, Ne=neExplore)
    return(result)
  }
  ## FUNCTION TO CALCULATE EFFECTIVE POPULATION SIZE USING MARKER INFORMATION
  neAnalysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input marker data file specified.")
  if(is.null(neExplore)){neExplore <- seq(10,100,10)}
  ############################
  # calculate the relationship matrix
  M <- object$data$geno 
  modificationsMarkers <- object$modifications$geno
  theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
  if(length(theresMatch) > 0){ # there's a modification file after matching the Id
    modificationsMarkers <- modificationsMarkers[theresMatch,]
    Markers <- cgiarBase::applyGenoModifications(M=M, modifications=modificationsMarkers)
  }else{ # there's no match of the modification file
    if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
  }
  ne <- neMarker(Markers, neExplore = neExplore, maxMarker=maxMarker, nSamples=nSamples, verbose=verbose)
  ## add metrics
  object$metrics <- rbind(object$metrics,
                               data.frame(module="neMarker",analysisId=neAnalysisId, trait= "none", environment="none",
                                          parameter=c( rep("allelesCovered", length(ne$allelesCovered)), rep("Ne", length(ne$Ne))  ),
                                          method="bootstrap",   value=c( ne$allelesCovered, ne$Ne ),
                                          stdError=c( ne$allelesCoveredSe, rep(NA, length(ne$Ne) ) )
                               )
  )
  ## add modeling
  currentModeling <- data.frame(module="neMarker", analysisId=neAnalysisId,trait="Ne", environment="across",
                                parameter=c("minNe","maxNe","maxMarker","nSamples","analysisIdForGenoModifications"), 
                                value=c(min(neExplore),max(neExplore),maxMarker,nSamples,analysisIdForGenoModifications))
  object$modeling <- rbind(object$modeling,currentModeling[,colnames(object$modeling)] )
  ## add status
  object$status <- rbind( object$status, data.frame(module="neMarker", analysisId=neAnalysisId))
  
  if(verbose){
    cat(paste("Your analysis id is:",neAnalysisId,"\n"))
  }
  return(object)
}
