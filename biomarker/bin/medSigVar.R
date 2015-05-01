# Author: Joseph Perez-Rogers
# Date: 2014-05-21
# Purpose: Script to compute the significance of the variance of a gene from the median variance of all genes
# Usage:
# Input:
# Output:

medSigVar <- function(x){
	variance <- unlist(apply(exprs(x),1,function(r){
		r <- r[-which.min(r)]
		r <- r[-which.max(r)]
		var(r)
	}))
	med.index <- sort(variance,decreasing=T,index.return=T)$ix[length(variance)/2]
	p <- unlist(apply(exprs(x),1,function(z){
		var.test(exprs(x)[med.index,],z,alternative="less")$p.value
	}))
	return(p)
}
