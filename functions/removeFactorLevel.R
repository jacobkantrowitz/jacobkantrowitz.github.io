# removeFactorLevel takes an expressionSet, a covariate, and one of that covariates levels
# and removes the level from the factor, including the name of the factor; this will 
# prevent a factor from being used in an analysis when one of that factor's levels has 0
# instances
removeFactorLevel <- function(exprData, covar, level)
{
	covarInd <- match(covar, varLabels(exprData))
#	removeInd <- 
	exprData2 <- exprData[, pData(exprData)[,covarInd]!=level]
	
	newLevels <- setdiff(levels(pData(exprData2)[,covarInd]), level)
	pData(exprData2)[,covarInd] <- factor(pData(exprData2)[,covarInd], levels=newLevels)
	
	return(exprData2)

}

