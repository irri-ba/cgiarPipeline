\name{staLMM}
\alias{staLMM}

\title{
Single-trial analysis pipeline using LMM
}
\description{
Perform the single-trial analysis per trait (see Details section).
}
\usage{
staLMM(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, 
    traitFamily=NULL,
    fixedTerm=c("1"),
    genoUnit = c("designation"),
    maxit=50,
    returnFixedGeno=TRUE,
    verbose=TRUE
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
\item{traitFamily}{
A character vector specifying the family distribution for each trait in case normal gaussian is not a good fit. 
}
\item{fixedTerm}{
A character vector specifying the covariates to be considered fixed.
}
\item{genoUnit}{
One or more of the following; designation, mother, father to indicate which column(s) should be considered the unit of genetic evaluation to compute BLUEs or BLUPs in the single trial analysis step.
}
\item{maxit}{
Maximum number of restricted maximum likelihood iterations to be run for each trait.
}
\item{returnFixedGeno}{
Whether BLUEs or BLUPs should be stored for the second stage.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

The genetic evaluation approach we use known as 'two-step' first analyze trait by trait and trial by trial to remove the spatial noise from experiments using experimental factors like blocking and spatial coordinates. Each trial is one level of the environment column (defined when the user matches the expected columns to columns present in the initial phenotypic input file). Genotype is fitted as both, fixed and random. The user defines which should be returned in the predictions table. By default genotype (designation column) predictions and their standard errors are returned.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Velazco, J. G., Rodriguez-Alvarez, M. X., Boer, M. P., Jordan, D. R., Eilers, P. H., Malosetti, M., & Van Eeuwijk, F. A. (2017). Modelling spatial trends in sorghum breeding field trials using a two-dimensional P-spline mixed model. Theoretical and Applied Genetics, 130, 1375-1392.

Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. (2018). Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial Statistics, 23, 52-71.

Boer M, van Rossum B (2022). _LMMsolver: Linear Mixed Model Solver_. R package version 1.0.4.9000.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno

final <- staLMM(
    phenoDTfile= result,
    analysisId=result$status$analysisId[1:],
    trait=c("Plant_Height_cm", "Ear_Height_cm"), # per trait
    traitFamily=c("gaussian(link = 'identity')", "gaussian(link = 'identity')"),
    fixedTerm=c("1"),
    genoUnit = c("designation"),
    maxit=50,
    returnFixedGeno=TRUE,
    verbose=TRUE
)
tail(final$predictions)

}

