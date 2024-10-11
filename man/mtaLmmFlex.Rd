\name{mtaLmmFlex}
\alias{mtaLmmFlex}

\title{
Multi-environmental trial analysis pipeline using LMM
}
\description{
Perform the multi-trial genetic evaluation using a random regression approach (see Details section).
}
\usage{
mtaLmmFlex(
    phenoDTfile= NULL, analysisId=NULL,
    analysisIdForGenoModifications=NULL,
    inputFormulation=NULL, envsToInclude=NULL,
    trait= NULL, traitFamily=NULL, useWeights=TRUE,
    entryTypeSubset=NULL,
    heritLB= 0.15,  heritUB= 0.95, meanLB=0, meanUB=Inf,
    modelType="blup", nMarkersRRBLUP=1000,
    deregress=FALSE, maxIters=50, batchSizeToPredict=500, 
    tolParInv=1e-4, verbose=TRUE
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the analysis.
}
\item{analysisIdForGenoModifications}{
A numeric id to link to the modifications table to be used to modify the marker information in case the user fit an rrBLUP or GBLUP model.
}
\item{inputFormulation}{
A character vector specifying the terms to be considered in the formulation.
}
\item{envsToInclude}{
A table of fields (rows) by traits (columns) to specify which fields should be included in the analysis for each trait.
}
\item{trait}{
A character vector with the traits that should be fitted.
}
\item{traitFamily}{
A character vector specifying the family distribution for each trait in case normal gaussian is not a good fit. 
}
\item{useWeights}{
A TRUE or FALSE statement to specify if the weights using the standard error should be applied. The default is TRUE since we assume the two-stage model.
}
\item{entryTypeSubset}{
A character vector indicating which entryTypes should be used in the modeling. For example, with multi-year data the user may want to subset to a specific testing stage. The default is NULL which implies that all entry types will be included in the analysis.
}
\item{heritLB}{
A numeric value or vector for trait heritability and reliablity lower-bound-threshold(s). Any trait by environment combination below this value will be excluded from the analysis.
}
\item{heritUB}{
A numeric value or vector for trait heritability and reliablity upper-bound-threshold(s). Any trait by environment combination above this value will be excluded from the analysis.
}
\item{modelType}{
The geentic evaluation model to used; blup, gblup, rrblup, pblup, ssblup.
}
\item{nMarkersRRBLUP}{
Maximum number of markers that should be used if an rrBLUP genetic evaluation method is selected. Under high density marker SNP panels the user may want to reduce the complexity of the marker-based model, the default is 1000. The markers will be randomly sampled without replacement.
}
\item{deregress}{
A TRUE/FALSE value to decide is the STA predictions should be deregressed. Only to be used if you fitted BLUPs in the STA step.
}
\item{maxIters}{
Maximum number of restricted maximum likelihood iterations to be run for each trait.
}
\item{batchSizeToPredict}{
A numeric value indicating how many individuals should be predicted per batch at the end of the genetic evaluation.
}
\item{tolParInv}{
A numeric value to indicate the tolerance value to add when a matrix is not invertible.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

The core algorithm of the genetic evaluation using the two-step approach is the multi-trial analysis. 

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Finlay, K. W., & Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme. Australian journal of agricultural research, 14(6), 742-754.

Henderson Jr, C. R. (1982). Analysis of covariance in the mixed model: higher-level, nonhomogeneous, and random regressions. Biometrics, 623-640.

Odegard, J., Indahl, U., Stranden, I., & Meuwissen, T. H. (2018). Large-scale genomic prediction using singular value decomposition of the genotype matrix. Genetics Selection Evolution, 50(1), 1-12.

}

\examples{
####=========================================####

data(DT_example)
data(inputFormulation)

result$status
result$metadata$pheno

res <- mtaLmmFlex(
    phenoDTfile= result, # analysis to be picked from predictions database
    analysisId=result$status$analysisId[2],
    analysisIdForGenoModifications=NULL,
    inputFormulation=inputFormulation, 
    envsToInclude=NULL,
    trait= c("Yield_Mg_ha","Plant_Height_cm"), traitFamily=NULL, useWeights=FALSE,
    heritLB= 0,  heritUB= 0.95,
    meanLB=0, meanUB=Inf,
    modelType="blup", # either "blup", "pblup", "gblup", "rrblup"
    nMarkersRRBLUP=1000,
    deregress=FALSE,  
    maxIters=50, batchSizeToPredict=500, tolParInv=1e-4,
    verbose=TRUE
)
tail(res$predictions)



}

