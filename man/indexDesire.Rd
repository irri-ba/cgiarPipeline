\name{indexDesire}
\alias{indexDesire}

\title{
Desire selection index pipeline 
}
\description{
Calculates the desire selection index given a table of traits and desire changes in traits (see Details section).
}
\usage{
indexDesire(
    phenoDTfile= NULL,
    analysisId=NULL, 
    environmentToUse=NULL,
    entryTypeToUse=NULL,
    trait= NULL,
    desirev = NULL, 
    scaled=TRUE, 
    verbose=TRUE 
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the selection index.
}
\item{environmentToUse}{
A character vector specifying the levels of environments to use in the index.
}
\item{entryTypeToUse}{
A character vector specifying the levels of entry types to use in the index.
}
\item{trait}{
A character vector with the traits that should be included in the selection index. Once the predictions are subset to the requested analysisId an additional subset for the traits of interest is performed.  
}
\item{desirev}{
A numeric vector with desired trait changes. This should be provided in original-scale units if the next argument (scaled) is set to FALSE, or provided as standard deviations if the scaled argument is set to TRUE. 
}
\item{scaled}{
A TRUE or FALSE statement to specify if the trait predictions from the data object should be scaled or not. 
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

The index follow the theory of Pesek and Baker (1964), where a matrix of traits Y of dimensions n rows (individuals) and t columns (traits) is multiplied by weights to compute the vector m of total merit of individuals:

m = Y w

At the same time the vector of weights 'w' is the cross-product of the inverse of the genetic variance covariance matrix (G) and the vector of desired change in traits (d)

w = G- d

Then the function row binds a new table of predictions to the predictions table of the data object together with some additional information to the modeling, metrics and status tables. 

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to selection indices. Canadian journal of plant science, 49(6), 803-804.

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear selection indices in modern plant breeding (p. 256). Springer Nature.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno

final <- indexDesire(
    phenoDTfile= result, # input data structure
    analysisId=result$status$analysisId[3], # analysis to be picked from predictions database
    trait= c("Plant_Height_cm", "Ear_Height_cm", "Yield_Mg_ha"), # traits to include in the index
    desirev = c(0, 0 , 1), # vector of desired values
    scaled=TRUE, # whether predicted values should be scaled or not
    verbose=TRUE # should we print logs or not
)
tail(final$predictions)

}

