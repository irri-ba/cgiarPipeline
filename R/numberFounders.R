numberFounders <- function(
    object= NULL,
    analysisIdForGenoModifications=NULL,
    maxNe=100, 
    maxMarker=1000, 
    nSamples=5,
    verbose=FALSE
){
  ## FUNCTION TO CALCULATE EFFECTIVE POPULATION SIZE USING MARKER INFORMATION
  neAnalysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input marker data file specified.")
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
  ne <- sommer::neMarker(Markers, maxNe=maxNe, maxMarker=maxMarker, nSamples=nSamples)
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
                                parameter=c("maxNe","maxMarker","nSamples","analysisIdForGenoModifications"), 
                                value=c(maxNe,maxMarker,nSamples,analysisIdForGenoModifications))
  object$modeling <- rbind(object$modeling,currentModeling[,colnames(object$modeling)] )
  ## add status
  object$status <- rbind( object$status, data.frame(module="neMarker", analysisId=neAnalysisId))
  
  if(verbose){
    cat(paste("Your analysis id is:",neAnalysisId,"\n"))
  }
  return(object)
}
