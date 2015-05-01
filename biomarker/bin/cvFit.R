# Author: Joe Perez-Rogers
# Date: 2014-07-14
# Purpose: A script to run k-fold cross-validation linear models and report pvalues and t-statistics

cvFit <- function(data,class,model,k=NULL,splits=NULL,reps=20,pct=0.8,verbose=FALSE){
  t.ranks <- c()
  p.ranks <- c()
  if(!is.null(k)){
    if(verbose){cat("Running the k-fold cross-validated linear model\n")}
    train.index <- createFolds(class,k=k,list=TRUE,returnTrain=FALSE)
    for(i in 1:k){
      # Running the linear model checking for each gene's association with the cancer label
      mod <- model[-train.index[[i]],]
      fit <- lmFit(data[,-train.index[[i]]],mod,na.action=na.exclude)
      fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
      fit2 <- ebayes(fit)
      p <- fit2$p.value[,ncol(fit2$p.value)]
      t.ranks <- cbind(t.ranks,fit.t[,ncol(fit.t)])
      p.ranks <- cbind(p.ranks,p)
    }
  } else if(!is.null(splits)){
    if(verbose){cat("Running the random sampling cross-validated linear model\n")}
    train.index <- createDataPartition(y=class,times=reps,p=pct)
    for(i in 1:reps){
      # Running the linear model checking for each gene's association with the cancer label
      mod <- model[train.index[[i]],]
      fit <- lmFit(data[,train.index[[i]]],mod,na.action=na.exclude)
      fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
      fit2 <- ebayes(fit)
      p <- fit2$p.value[,ncol(fit2$p.value)]
      t.ranks <- cbind(t.ranks,fit.t[,ncol(fit.t)])
      p.ranks <- cbind(p.ranks,p)
    }
  }
  return(list("t"=rowMeans(t.ranks),"p"=rowMeans(p.ranks)))  
}