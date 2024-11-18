traitTransformation <- function(
    object= NULL,
    trait=NULL, # per trait
    transformation = NULL,
    verbose=FALSE
){
  if(is.null(object)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be converted", call. = FALSE)}
  if(length(trait)==0){stop("Please select at least one trait conversion to apply this function.", call. = FALSE)}
  if(length(trait) != length(transformation)){stop("The vector of transformations required need to be as many as the number of traits and viceversa.", call. = FALSE)}
  cbrt <- function(x){return(x^(1/3))}
  traitOrig <- trait
  transAnalysisId <- as.numeric(Sys.time())
  ###################################
  # loading the dataset
  mydata <- object$data$pheno
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  transformation <- transformation[which(traitOrig %in% trait)]
  #####################################
  # transformation
  counter <- 1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,transformation[counter], sep = "_")] <- do.call(transformation[counter],list(mydata[,iTrait]))
    counter <- counter + 1
  }
  object$data$pheno <- mydata
  ##########################################
  ## update databases
  ## status
  newStatus <- data.frame(module="transP", analysisId=transAnalysisId, analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  # modeling
  currentModeling <- data.frame(module="transP", analysisId=transAnalysisId,trait=trait, environment=NA,
                                parameter=transformation,
                                value=NA)
  object$modeling <- rbind(object$modeling,currentModeling )
  return(object)
}

