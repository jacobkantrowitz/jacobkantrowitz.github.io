# Author: Joe Perez-Rogers
# Date: 2014-05-11
# Script to set a random seed based on the run index or date
# Usage: setRandomSeed(x,seed,print=FALSE)
# Input: x=an integer indicating which index is currently being run, seed=numeric value of the desired seed to be set, print=whether to print out the seed or not
# Output:

setRandomSeed <- function(x=NULL,seed=NULL,print=FALSE){
	if(!is.null(seed)){
		s <- seed
		set.seed(s)
	} else if(is.null(x)){
		s <- as.numeric(format(Sys.time(),"%d%H%M%S"))
		set.seed(s)
	} else if(!is.null(x)){
		s <- as.numeric(format(Sys.time(),"%d%H%M%S"))
		set.seed(s)
		s <- round(runif(x)[x]*100000000)
		set.seed(s)
	}
	if(print==TRUE){
		cat(paste("Seed: ",s,"\n",sep=""))
	}
}

