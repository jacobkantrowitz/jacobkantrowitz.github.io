# cleanNAForAnalysis takes a list of covariates and removes all subjects from the 
# given expressionSet who have NAs in any of the covariates given
# cleanNAForAnalysis takes an expressionSet object and a list of covariates
# cleanNAForAnalysis returns the cleaned expressionSet
cleanNAForAnalysis <- function(exprData, covariates)
{
	tmp <- exprData
	for(i in 1:length(covariates))
	{
		tmp <- removeCovariateNAs(tmp, covariates[i])
	}
	return(tmp)
}