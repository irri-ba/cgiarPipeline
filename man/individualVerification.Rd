\name{individualVerification}
\alias{individualVerification}

\title{
Verification of a cross using markers
}
\description{
This function takes a matrix of genotypes for males, another for females, and another for progeny and calculates the probability of matching the expected genotype. Additional metrics like hetrozygosity and deviation from expectations are returned.
}
\usage{
individualVerification( 
  object= NULL,
  analysisIdForGenoModifications= NULL,
  markersToBeUsed=NULL,
  colsForExpecGeno=NULL, ploidy=2
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
  \item{colsForExpecGeno}{
Columns that define the average of genotype markers to be used as expected genotype (e.g., F1 cross, Backcross, etc.).
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

result$data$pedigree$Hybrid_Parent1 <- sample(result$data$pedigree$Hybrid)
result$data$pedigree$Hybrid_Parent2 <- sample(result$data$pedigree$Hybrid)

res <- individualVerification( 
  object= result,
  analysisIdForGenoModifications= result$status$analysisId[8],
  markersToBeUsed=NULL,
  colsForExpecGeno=NULL, ploidy=2
)

tail(res$modeling)


}

