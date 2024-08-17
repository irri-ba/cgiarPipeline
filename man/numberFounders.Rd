\name{numberFounders}
\alias{numberFounders}

\title{
Calculates the number of founders based on allele coverage.
}
\description{
Uses a marker matrix to approximate the effective population size (Ne) by discovering how many individuals are needed to sample all possible alleles in a population..
}
\usage{
numberFounders(
    object= NULL,
    analysisIdForGenoModifications=NULL,
    maxNe=100, 
    maxMarker=1000, 
    nSamples=5,
    verbose=FALSE
)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisIdForGenoModifications}{
A numeric id to link to the modifications to be applied to the marker information.
}
\item{maxNe}{
The maximum number of potential founders to try.
}
\item{maxMarker}{
The maximum number of markers to sample. 
}
\item{nSamples}{
number of individuals to sample for the Ne calculation.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

Uses a marker matrix to approximate the effective population size (Ne) by discovering how many individuals are needed to sample all possible alleles in a population.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Lush, J. L. (2013). Animal breeding plans. Read Books Ltd.

Mrode, R. A. (2014). Linear models for the prediction of animal breeding values. Cabi.

}

\examples{
####=========================================####

data("DT_example")

result$status

result2 <- numberFounders(
  object= result,
  analysisIdForGenoModifications=result$status$analysisId[8],
  maxNe=100, 
  maxMarker=1000, 
  nSamples=10,
  verbose=FALSE
)

metrics <- result2$metrics
status <- result2$status

metrics <- metrics[which(metrics$analysisId == status$analysisId[9]),]
metrics <- metrics[which(metrics$parameter=="allelesCovered"),]
plot(metrics$value)
max(metrics$value)

}

