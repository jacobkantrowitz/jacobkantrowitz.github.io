# Author: Joe Perez-Rogers
# Date: 2014-05-12
# Script to run gene filtration algorithms on an expression set
# Usage:
# Input:
# Output:

runGeneFilter <- function(x,y=NULL,method=NULL){
	if(is.null(method)){
		ret <- x
	} else if(method=="lmFit"){
		design <- model.matrix(~CANCER,data=pData(x))
		fit <- lmFit(exprs(x),design)
		fit2 <- ebayes(fit)
		p <- as.data.frame(fit2$p.value)
		ret <- p
	} else if(method=="corr"){
		if(is.null(y)){
			cat("Error: You must provide two expression set objects to compute correlations\n")
			cat("Returning a NULL object\n")
			ret <- NULL
		} else if((dim(x)!=dim(y))[1] | (dim(x)!=dim(y))[2]){
			cat("Error: Your expression sets are of differing dimensions\n")
			cat("Returning a NULL object\n")
			ret <- NULL
		}
	}
	return(ret)
}
