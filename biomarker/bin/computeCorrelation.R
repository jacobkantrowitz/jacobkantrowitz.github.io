# Author: Joe Perez-Rogers
# Date: 2014-10-07
# Purpose: Compute the per-gene correlation coefficients between two expression sets

computeCorrelation <- function(x,y,id="rows",return.all=F){
  continue <- TRUE
  # make sure that the expression sets contain paired data in the same order
  if(id=="rows"){
    x <- x[,row.names(pData(x))%in%row.names(pData(y))]
    idx <- na.omit(match(row.names(pData(x)),row.names(pData(y))))
    y <- y[,idx]
  } else {
    if(is.numeric(id)){
      x <- x[,as.character(pData(x)[,id])%in%as.character(pData(y)[,id])]
      idx <- na.omit(match(pData(x)[,id],pData(y)[,id]))
      y <- y[,idx]
    } else {
      continue <- FALSE
    }
  }
  if(continue){
    # finding the top correlated genes by p-value (p<0.005)
    correlations <- sapply(1:nrow(exprs(x)),function(z){
      correlation <- cor.test(exprs(x)[z,],exprs(y)[z,],alternative="greater");
      c(correlation$estimate,correlation$p.value)
      })
    correlations <- as.data.frame(t(correlations))
    colnames(correlations) <- c("R","p.value")
    row.names(correlations) <- row.names(exprs(x))
    if(return.all){
      return(list("correlations"=correlations,"xdim"=dim(x)))
    } else {
      return(correlations)
    }
  } else {
    return(NULL)
  }
}