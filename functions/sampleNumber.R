# sampleNumber returns the number of samples from an expressionSet object (biobase)
# sampleNumber takes as input an expressionSet object 
sampleNumber <- function(exprSet)
{
	return(dim(exprSet)[2])
}
