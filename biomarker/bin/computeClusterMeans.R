# Author: Joseph Perez-Rogers
# Date: 2014-06-12
# Purpose: Script to compute the mean values of clusters from consensus clustering
# Input:
# Output: 
# Usage: 

computeClusterMeans <- function(x,classes){
  #source("../../bin/assert.R")
  nclust <- length(table(classes))
  #assert(nrow(x)==length(classes),"You must have a class label for each nrow(x)")
  cmeans <- c()
  for(c in 1:nclust){
    l <- length(which(classes==c))
    if(l>1){
      cm <- colMeans(x[classes==c,])
    } else {
      cm <- x[classes==c,]
    }
    cmeans <- cbind(cmeans,cm)
  }
  return(cmeans)
}