\name{rggPiepho}
\alias{rggPiepho}

\title{
Realized genetic gain pipeline (Laidig's method).
}
\description{
Performs a realized genetic gain estimation using the two-stage analysis (Laidig et al., 2014).
}
\usage{
rggPiepho(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, 
    deregress=FALSE,
    yearsToUse=NULL,
    entryTypeToUse=NULL,
    verbose=TRUE,
    traitFamily=NULL,
    sampleN=50, 
    bootstrappingN=10,
    forceRules=TRUE
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the analysis.
}
\item{trait}{
A character vector with the traits that should be fitted.
}
\item{deregress}{
Should we deregress the estimates by dividing over the reliability before performing the realized genetic gain calculation.
}
\item{yearsToUse}{
Selection of the years of origin associated to the tested material to use in the calculation.
}
\item{entryTypeToUse}{
A selection of entry types to use for the realized genetic gain calculation.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}
\item{traitFamily}{
Family of the response variable.
}
\item{sampleN}{
Number of individuals to sample per environment to perform the analysis. Included for computational efficiency.
}
\item{bootstrappingN}{
Number of bootstrap samples.
}
\item{forceRules}{
A TRUE/FALSE statement declaring if we should enforce the rule of Excellence in Breeding of minimum 6 years of data to perform realized genetic gain analysis.
}

}
\details{

In order to monitor the efficacy of genetic evaluation across cycles of selection, the realized genetic gain is the preferred process. This option aims to calculate the realized genetic gain using the methods from Laidig et al. (2014). The method uses across-environment means from multiple years of data that have been adjusted based on a good connectivity to then fit a regression of the form means~year.of.origin. In case the means used are BLUPs these can be deregressed.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Laidig, F., Piepho, H. P., Drobek, T., & Meyer, U. (2014). Genetic and non-genetic long-term trends of 12 different crops in German official variety performance trials and on-farm yield trends. Theoretical and Applied Genetics, 127, 2599-2617.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno
result$data$pedigree$Year <- sample(2011:2023, length(result$data$pedigree$Year), replace = TRUE)

final <- rggPiepho(
    phenoDTfile= result,
    analysisId=result$status$analysisId[2],
    trait="Ear_Height_cm", # per trait
    deregress=FALSE,
    yearsToUse=NULL,
    entryTypeToUse=NULL,
    verbose=TRUE,
    traitFamily=NULL,
    sampleN=50, # max number of individuals per environment to sample
    bootstrappingN=10,
    forceRules=TRUE
)

tail(final$metrics,10)

}

