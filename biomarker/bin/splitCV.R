# Author: Joe Perez-Rogers
# Date: 2014-05-11
# Script to create k-fold cross-validation splits of an expression set object
# Usage: splitCV()
# Input: x=an expression set object, pct=the percentage of samples you want in your test set
# Output: 

splitCV <- function(x,pct=0.20){
	data <- c(1:nrow(pData(x)))
	refs <- sample(data,round(nrow(pData(x))*pct), replace=F)
	return(refs)
}
