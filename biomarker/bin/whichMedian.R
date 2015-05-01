# Author: Joe Perez-Rogers
# Date: 2014-11-24
# Script: whichMedian.R

whichMedian <- function(x){
	w <- which(x==median(x))
	if(length(w)>0){
		med.idx <- w
	} else {
		diffs <- abs(x-median(x))
		med.idx <- which.min(diffs)
	}
	return(med.idx)
}