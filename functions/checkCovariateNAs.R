# checkCovariateNAs tests for NAs in a set of given phenotypes within an expressionSet
# checkCovariateNAs takes an expressionSet object and a list of covariates
# (cannot handle lists of lists of covariates)
# checkCovariateNAs returns a vector indicating the number of NAs present in each of 
# the given covariates

# this function should now only be called via printNACovariates, which in turn should only be called
# from check_covariates. This order enforcement will prevent checkCovariateNAs from having to deal
# with lists of lists of covariates (LoL), which check_covariates already deals with recursively. 
# this function will now be called as part of the base case in check_covariates, when that function
# has determined it is dealing with only a single list of covariates and nota LoL
# This function also assumes that the list of covariates accepted are all present in pData(exprData)
checkCovariateNAs <- function(exprData, covariates)
{
	cat("Checking covariates for NA values...............\n")
#	test <- check_covariates(exprData, covariates)
	
#	if(test==0)
#	{
	numberNA <- function(x){number_of(is.na(x))}
	covsInds <- match(covariates, varLabels(exprData))
	covsNumNA <- sapply(pData(exprData)[covsInds], numberNA)
	
	return(covsNumNA)
	#numeric(length(covariates))
	
	# check to see if any of the pData(exprData) fields have NAs
	# report those fields that have NAs
#	}
}