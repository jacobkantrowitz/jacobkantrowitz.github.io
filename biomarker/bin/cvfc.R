# Author: Joseph Perez-Rogers
# Date: 2014-06-10
# Purpose: My implementation of FCROS method for fold change calculations
# Input:
# Output:
# Usage:

cvfc <- function(x,labels,n=1000,s=5,seed=sample(1:10**5,1),verbose=FALSE){
  set.seed(seed)
  pct <- 0
  labels <- as.factor(labels)
  if(nlevels(labels)!=2){stop("Your labels must have only two levels")}
  ranks <- rep(0,nrow(x))
  c1 <- which(labels==levels(labels)[1])
  c2 <- which(labels==levels(labels)[2])
  for(i in 1:n){
    fc <- rowMeans(x[,sample(c1,s)])-rowMeans(x[,sample(c2,s)])
    idx <- sort(fc,decreasing=T,index.return=T)$ix
    tmp <- cbind(idx,"score"=(1:nrow(x)))
    ranks <- ranks + tmp[order(tmp[,1]),2]
    if(verbose){
      if(i%in%seq(from=0,to=n,by=(n/10))){
        pct <- pct + 10
        cat(paste(pct,"% Complete\n",sep=""))
      }
    }
  }
  ranks <- ranks/n
  return(ranks)
}