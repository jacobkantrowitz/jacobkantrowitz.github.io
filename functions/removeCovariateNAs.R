# removeCovariateNAs removes the samples from an expressionSet with NA values in the 
# given covariate
# removeCovariateNAs takes an expressionSet object and a single covariate
# removeCovariateNAs returns the expressionSet but with those NA samples removed
removeCovariateNAs <- function(exprData, covar)
{
	covarInd <- match(covar, varLabels(exprData))
	noneNAinds <- which(!is.na(pData(exprData)[,covarInd]))
	return(exprData[,noneNAinds])

}