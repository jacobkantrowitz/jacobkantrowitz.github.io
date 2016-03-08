# traitNumber returns the number of traits from a pData(expressionSet) object (biobase)
# traitNumber takes as input an expressionSet object 
traitNumber <- function(exprSet)
{
	return(length(varLabels(exprSet)))
}