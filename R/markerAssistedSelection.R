
markerAssistedSelection <- function(
    object= NULL,
    analysisIdForGenoModifications= NULL,
    markersToBeUsed=NULL,
    positiveAlleles=NULL,
    desire=NULL, ploidy=2
){

  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")
  if (is.null(analysisIdForGenoModifications)) stop("No geno clean file specified.")
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if("effectType" %!in% colnames(object$predictions) ){
    object$predictions$effectType <- NA
  }
  # get markers
  Markers <- object$data$geno
  if(is.null(Markers)){stop("This function requires your object to have marker information.", call. = FALSE)}
  # apply modifications
  if(!is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
    modificationsMarkers <- object$modifications$geno
    theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
    if(length(theresMatch) > 0){ # there's a modification file after matching the Id
      modificationsMarkers <- modificationsMarkers[theresMatch,]
      Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
    }
  }else{
    mynames <- rownames(Markers)
    Markers <- apply(Markers,2,sommer::imputev)
    rownames(Markers) <- mynames
  }
  if (is.null(markersToBeUsed)){markersToBeUsed <- 1:ncol(Markers)}else{markersToBeUsed <- intersect(colnames(Markers),markersToBeUsed)}
  Markers <- Markers[,markersToBeUsed]
  ## extract marker matrices and reference alleles
  names(positiveAlleles) <- markersToBeUsed
  X <- result$metadata$geno[markersToBeUsed,]
  positiveAllelesN <- vector(mode = "numeric", length = length(positiveAlleles))
  for(j in 1:nrow(X)){
    positiveAllelesN[j] <- ifelse(X[j,c("refAllele")] == positiveAlleles[X$marker[j]], 2, 0 )
  }
  # if user doesn't express a desire we weight by 1 minus the frequency
  if(is.null(desire) ){
    p <- vector(mode="numeric",length = length(positiveAllelesN))
    for(k in 1:ncol(Markers)){
      ttb <- table(0:ploidy)-1 # table of zeros for dosages
      tto <- table(Markers[,k])
      ttb[names(tto)] <- ttb[names(tto)] + tto
      n <- sum(ttb)
      p[k] <- ( ifelse(positiveAllelesN[k]==0, ttb[1],ttb[length(ttb)]) + (ttb[2:(length(ttb)-1)] / factorial(2:(length(ttb)-1)) ) ) / n
    }
    desire <- (1 - p) # distance to fixation
  }else{ # we apply both, frequencies and weights from user
    # if(var(desire)==0){
    p <- vector(mode="numeric",length = length(positiveAllelesN))
    for(k in 1:ncol(Markers)){
      ttb <- table(0:ploidy)-1 # table of zeros for dosages
      tto <- table(Markers[,k])
      ttb[names(tto)] <- ttb[names(tto)] + tto
      n <- sum(ttb)
      p[k] <- ( ifelse(positiveAllelesN[k]==0, ttb[1],ttb[length(ttb)]) + (ttb[2:(length(ttb)-1)] / factorial(2:(length(ttb)-1)) ) ) / n
    }
    desire <- (1 - p)*desire # distance to fixation
    # }
  }

  G <- cov(Markers)
  Gi <- solve(G + diag(1e-6,ncol(G),ncol(G)))

  w <- Gi %*% desire

  merit <- Markers %*% w

  ###############
  # other tables
  object$status <- rbind( object$status, data.frame(module="mas", analysisId=analysisId))
  ## modeling
  modeling <- data.frame(module="mas",  analysisId=analysisId, trait=c(markersToBeUsed,markersToBeUsed,"inputObject"), environment="across",
                         parameter= c(rep("markerUsed",length(markersToBeUsed)),rep("desire", length(markersToBeUsed)), "analysisId"  ),
                         value= c(markersToBeUsed,round(desire,6),ifelse(is.null(analysisIdForGenoModifications),NA,analysisIdForGenoModifications))
  )
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }
  ##############
  ## predictions
  preds <- data.frame(module="mas",  analysisId=analysisId, pipeline= "unknown",
                      trait="masMerit", gid=1:nrow(Markers), designation=rownames(Markers),
                      mother=NA,father=NA, entryType="test_entry", effectType="designation",
                      environment="across", predictedValue=as.numeric(merit), stdError=NA,
                      reliability=NA
  )
  if(is.null(object$predictions)){
    object$predictions <-  preds
  }else{
    object$predictions <- rbind(object$predictions, preds[, colnames(object$predictions)])
  }

  return(object)

}




