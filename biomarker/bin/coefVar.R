# Author: Joseph Perez-Rogers
# Date: 2014-05-21
# Purpose: Script to compute the coefficient of variation of a gene
# Usage:
# Input:
# Output:

coefVar <- function(x){
	sd(x)/mean(x)
}
