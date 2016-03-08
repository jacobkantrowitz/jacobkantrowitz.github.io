# removeFactorLevel takes an expressionSet, a covariate, and one of that covariates levels
# and removes the level from the factor, including the name of the factor; this will 
# prevent a factor from being used in an analysis when one of that factor's levels has 0
# instances
#' Removes samples from an expressionSet from a specific level in a specific phenotype field
#' 
#' @param exprData the expressionSet object from which to remove samples matching the given criteria
#' @param covar the phenotype data field to match and on which to base sample removals
#' @param level the level of the covar phenotype data to select for removal from the expressionSet
#' @return an expressionSet object with no samples matching the level in the covar
#' @description This function takes an expressionSet object, a factor-type covariate and a level from that covariates and removes all samples of the level from the expressionSet. Moreover, the phenotype is releveled so that it can be used in modeling easily.
#' @examples
#' removeFactorLevel(eSet, "GENDER", "MALE") # GENDER must be a phenotype field and MALE one of its levels
removeFactorLevel <- function(exprData, covar, level)
{
	covarInd <- match(covar, varLabels(exprData))
#	removeInd <- 
	exprData2 <- exprData[, pData(exprData)[,covarInd]!=level]
	
	newLevels <- setdiff(levels(pData(exprData2)[,covarInd]), level)
	pData(exprData2)[,covarInd] <- factor(pData(exprData2)[,covarInd], levels=newLevels)
	
	return(exprData2)

}

