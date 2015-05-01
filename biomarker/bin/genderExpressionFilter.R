# Author: Joe Perez-Rogers
# Date: 2014-10-07
# Purpose: Script to filter genes based on the expression of y-linked genes in females

genderExpressionFilter <- function(x,pop.threshold=0.05,sd.threshold=1.5){
  # Removing genes that are not expressed in at least pop.theshold % of samples
  sex.genes <- c("DDX3Y","KDM5D","RPS4Y1","USP9Y","UTY")
  expression <- as.vector(exprs(x)[sex.genes,which(x$GENDER==1)])
  expression <- expression[-which(expression>6)]
  noise.cutoff <- mean(expression)+sd.threshold*sd(expression)
  pct <- round(ncol(exprs(x))*pop.threshold)
  signal <- apply(exprs(x),1,function(x){length(which(x>noise.cutoff))>pct})
  noisy.genes <- row.names(exprs(x))[which(!signal)]
  newx <- x[!(row.names(exprs(x))%in%noisy.genes),]
  return(noisy.genes)
}