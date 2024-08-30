\name{markerAssistedSelection}
\alias{markerAssistedSelection}

\title{
Apply Marker Assisted Selection
}
\description{
This function takes a matrix of genotypes for individuals, weights for each marker, and builds an index to.
}
\usage{
markerAssistedSelection( 
    object= NULL,
    analysisIdForGenoModifications= NULL,
    markersToBeUsed=NULL,
    positiveAlleles=NULL,
    desire=NULL, ploidy=2
)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisIdForGenoModifications}{
A numeric id to link to the modifications table to be used to modify the marker information in case the user fit an rrBLUP or GBLUP model.
}
  \item{markersToBeUsed}{
Name of markers to be subset.
}
  \item{positiveAlleles}{
A character vector indicating what is the positive allele for the trait of interest.
}
\item{desire}{
A numeric vector indicating the desire value for the desire index. Good initial values can be the coefficients of variation of your traits times their distance to fixation.
}
\item{ploidy}{
A numeric value indicating the ploidy of the individuals.
}

}
\details{

This function takes the markers matrices and performs the cross verification.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{metrics}{a table of individuals' metrics.}
\item{matchMat}{the matrix of probabilities.}
}
}

\references{

None

}

\examples{
####=========================================####

data(DT_example)

marks <- colnames(result$data$geno)[1:10]
positiveAlleles = result$metadata$geno[marks,"refAllele"]

res <- markerAssistedSelection( 
    object= result,
    analysisIdForGenoModifications= result$status$analysisId[8],
    markersToBeUsed=marks,
    positiveAlleles=positiveAlleles,
    desire=NULL, ploidy=2
)

tail(res$predictions)  

}

