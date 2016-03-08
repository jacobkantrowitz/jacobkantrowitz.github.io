# sampleNumber returns the number of samples from an expressionSet object (biobase)
# sampleNumber takes as input an expressionSet object 
#' The number of samples in an expressionSet object
#' 
#' @param exprSet expressionSet object
#' @return numeric giving the number of samples in exprSet
#' @examples
#' sampleNumber(eset)
sampleNumber <- function(exprSet)
{
	return(dim(exprSet)[2])
}
