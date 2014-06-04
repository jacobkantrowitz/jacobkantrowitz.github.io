# printNACovariates prints the number of NAs present in each covariate of a given list
# it calls 'checkCovariateNAs' to calculate the number of NAs
# printNACovariates takes an expressionSet object and a list of covariates

# this function is being edited (as of 2014-06-04) to now be called by check_covariates and will
# now not itself call check_covariates. This is necessary because of the way in which lists of lists 
# of covariates can be sent for analysis. This change will make it so that lists of lists (LoL) are sent
# to check_covariates, which checks for the phenotype's presence in the exprData object and then checks
# that the phenotype doesn't have any NAs present in the exprData object. 

# this function assumes that the list of covariates accepted are all present in the pData(exprData)
# should write error checking to make sure this is the case...?
printNACovariates <- function(exprData, covariates)
{
	covNAs <- checkCovariateNAs(exprData, covariates)
	numsNA <- covNAs[covNAs>0]
	if(length(numsNA)>0)
	{
		for(i in 1:length(numsNA))
		{
			cat(paste("Trait ", names(numsNA)[i], " has ", numsNA[i], " NA values\n", sep=""))
		}
	}
  else
  {
   cat("There are 0 NA values present in any of the phenotype fields provided.\n\n")
  }
	return(number_of(covNAs>0))
}

