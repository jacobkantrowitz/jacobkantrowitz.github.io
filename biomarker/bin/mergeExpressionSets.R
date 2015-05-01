# Author: Joe Perez-Rogers
# Date: 2014-05-07
# Merge two or more expression set objects into a single expression set object
# Usage: mergeExpressionSets(es1,es2)
# Input: Two expression sets separated by a comma or a list of expression sets in the form c(es1,es2,es3,...)
# Output: A single expression set object

mergeExpressionSets <- function(x,y=NULL){
	require(affy)
	# N > 2 expression sets to merge
	if(length(x)>1){
		exit = FALSE
		tmp <- x[[1]]
		for(i in 2:length(x)){
			if(exit==TRUE){
				tmp <- NULL
			} else {
				# Check dimensions of expression matrices
				if(dim(exprs(tmp))[1]!=dim(exprs(x[[i]]))[1]){
					cat("Your expression set objects have different exprs() dimensions\n")
					exit <- TRUE
				# Check dimensions of phenotypic matrices
				} else if(dim(pData(tmp))[2]!=dim(pData(x[[i]]))[2]){
					cat("Your expression set objects have different pData() dimensions\n")
					exit <- TRUE
				} else {
					exprs(tmp) <- cbind(exprs(tmp),exprs(x[[i]]))
					pData(tmp) <- rbind(pData(tmp),pData(x[[i]]))
					exit <- FALSE
				}
			}
		}
		ret <- tmp
	
	# N = 2 expression sets to merge
	} else if(!is.null(y)){
		tmp <- x
		# Check dimensions of expression matrices
		if(dim(exprs(tmp))[1]!=dim(exprs(y))[1]){
			cat("Your expression set objects have different exprs() dimensions\n")
			tmp <- NULL
		# Check dimensions of phenotypic matrices
		} else if(dim(pData(tmp))[2]!=dim(pData(y))[2]){
			cat("Your expression set objects have different pData() dimensions\n")
			tmp <- NULL
		} else {
			exprs(tmp) <- cbind(exprs(tmp),exprs(y))
			pData(tmp) <- rbind(pData(tmp),pData(y))
		}
		ret <- tmp
	} else {
		ret <- x
	}
	return(ret)
}



