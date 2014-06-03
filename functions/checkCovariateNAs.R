# checkCovariateNAs tests for NAs in a set of given phenotypes within an expressionSet
# checkCovariateNAs takes an expressionSet object and a list of covariates
# checkCovariateNAs returns a vector indicating the number of NAs present in each of 
# the given covariates
checkCovariateNAs <- function(exprData, covariates)
{
	cat("Checking covariates for NA values...............\n")
	test <- check_covariates(exprData, covariates)
	
	if(test==0)
	{
		numberNA <- function(x){number_of(is.na(x))}
		covsInds <- match(covariates, varLabels(exprData))
		covsNumNA <- sapply(pData(exprData)[covsInds], numberNA)
		
		return(covsNumNA)
		#numeric(length(covariates))
		
		# check to see if any of the pData(exprData) fields have NAs
		# report those fields that have NAs
	}
}