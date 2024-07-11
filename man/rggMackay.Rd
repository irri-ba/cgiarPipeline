\name{rggMackay}
\alias{rggMackay}

\title{
Realized genetic gain pipeline (Mackay's method).
}
\description{
Performs a realized genetic gain estimation using the two-stage analysis (Mackay et al., 2011).
}
\usage{
rggMackay(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    deregressWeight=1,
    partition=FALSE,
    yearsToUse=NULL,
    entryTypeToUse=NULL,
    verbose=TRUE,
    forceRules=TRUE,
    propTopIndsPerYear=1
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
\item{deregressWeight}{
Should any weight be applied to the deregressed value (not recommended but available). 
}
\item{partition}{
When very few years of data are present this option will allow the user to calculate the gain for all 2-year combinations and then average these rates.
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
\item{forceRules}{
A TRUE/FALSE statement declaring if we should enforce the rule of Excellence in Breeding of minimum 6 years of data to perform realized genetic gain analysis.
}
\item{propTopIndsPerYear}{
A value between 0 and 1 to indicate the proportion of individuals to pick per year in case we want to focus the analysis on the top performing individuals per years.
}

}
\details{

In order to monitor the efficacy of genetic evaluation across cycles of selection, the realized genetic gain is the preferred process. This option aims to calculate the realized genetic gain using the methods from Mackay et al. (2011). The method uses across-environment means from multiple years of data that have been adjusted based on a good connectivity to then fit a regression of the form means~year.of.origin. In case the means used are BLUPs these can be deregressed.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Mackay, I., Horwell, A., Garner, J., White, J., McKee, J., & Philpott, H. (2011). Reanalyses of the historical series of UK variety trials to quantify the contributions of genetic and environmental factors to trends and variability in yield over time. Theoretical and Applied Genetics, 122, 225-238.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno
result$data$pedigree$Year <- sample(2011:2023, length(result$data$pedigree$Year), replace = TRUE)

final <- rggMackay(
    phenoDTfile= result,
    analysisId=result$status$analysisId[3],
    trait="Yield_Mg_ha", # per trait
    deregressWeight=1,
    deregress=FALSE,
    partition=FALSE,
    yearsToUse=NULL,
    entryTypeToUse=NULL,
    verbose=TRUE,
    forceRules=TRUE,
    propTopIndsPerYear=1
)

tail(final$metrics)

}

