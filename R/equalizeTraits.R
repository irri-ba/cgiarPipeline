equalizeTraits <- function(object, traits, newName=NULL){
  if(is.null(object)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(newName == ""){stop("Please provide a name for equalized trait", call. = FALSE)}
  if(is.null(newName)){stop("Please provide a name for equalized trait", call. = FALSE)}
  if(newName %in% colnames(object$data$pheno)){
    stop(paste0("Column '", newName,"' already exist. Please provide a different name for equalized trait."), call. = FALSE)
  }
  mydata <- object$data$pheno
  newValue <- apply(mydata[,traits, drop=FALSE],1, function(x){mean(x, na.rm=TRUE)})
  transAnalysisId <- as.numeric(Sys.time())
  # if(is.null(newName) ){
  #   newName <- paste(traits,collapse = "_")
  # }
  # if( newName == "" ){
  #   newName <- paste(traits,collapse = "_")
  # }
  newName <- gsub(" ","",newName)
  mydata[,newName] <- newValue
  object$data$pheno <- mydata
  object$metadata$pheno <- rbind( object$metadata$pheno, data.frame(parameter="trait", value=newName))
  ##########################################
  ## update databases
  ## status
  newStatus <- data.frame(module="transE", analysisId=transAnalysisId, analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  # modeling
  currentModeling <- data.frame(module="transE", analysisId=transAnalysisId, trait="inputObject", environment="general",
                                parameter=newName,
                                value=paste0("(",paste(traits, collapse = ' + '),")/",length(traits)))
  object$modeling <- rbind(object$modeling,currentModeling )
  return(object)
}
