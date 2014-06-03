# covariateClasses returns a list of the class of each of the covariates
# covariateClasses takes an expressionSet object and a list of covariates
# the covariates are checked for presence/absence in the pData of the expressionSet
covariateClasses <- function(exprData, covariates)
{
	# before checking for all classes, check that all covariates are in pData(exprData)
	test <- check_covariates(exprData, covariates)
	
	# continue if all covariates are in pData(exprData)
	if (test==0)
	{
		numberCovariates <- length(covariates)
		classes	<- matrix(nrow=numberCovariates, ncol=2)
		classes[,1] <- as.vector(covariates)
		covsInds <- match(covariates, varLabels(exprData))
		classes[,2] <- as.vector((sapply(pData(exprData)[covsInds], class)))
		
		return(classes)
	}
	else
	{
		cat("\nSome of the covariates provided were not found in the phenotype data.\n",
			"Please edit the indicated covariates and resubmit.\n\n")
	}
	
}