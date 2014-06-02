# featureNumber returns the number of features (e.g. genes) from an expressionSet object (biobase)
# featureNumber takes as input an expressionSet object 
featureNumber <- function(exprSet)
{
	return(dim(exprSet)[1])
}