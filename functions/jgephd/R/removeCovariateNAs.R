# removeCovariateNAs removes the samples from an expressionSet with NA values in the 
# given covariate
# removeCovariateNAs takes an expressionSet object and a single covariate
# removeCovariateNAs returns the expressionSet but with those NA samples removed
#' Returns the expression set after removing the samples with NAs in the given phenotype data field
#' 
#' @param exprData an expressionSet object with phenotype data included
#' @param covar the name of the phenotype field in exprData from which to remove NA values
#' @return expressionSet object with  NA values in the covar field having been removed
#' @examples
#' removeCovariateNAs(eSet, "GENDER")
removeCovariateNAs <- function(exprData, covar)
{
	covarInd <- match(covar, varLabels(exprData))
	noneNAinds <- which(!is.na(pData(exprData)[,covarInd]))
	return(exprData[,noneNAinds])

}