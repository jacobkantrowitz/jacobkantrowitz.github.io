# printNACovariates prints the number of NAs present in each covariate of a given list
# it calls 'checkCovariateNAs' to calculate the number of NAs
# printNACovariates takes an expressionSet object and a list of covariates
printNACovariates <- function(exprData, covariates)
{
	covNAs <- checkCovariateNAs(exprData, covariates)
	numsNA <- covNAs[covNAs>0]
	cat(numsNA)
	if(length(numsNA)>0)
	{
		for(i in 1:length(numsNA))
		{
			cat(paste("Trait ", names(numsNA)[i], " has ", numsNA[i], " NA values\n", sep=""))
		}
	}
	return(number_of(covNAs>0))
}

