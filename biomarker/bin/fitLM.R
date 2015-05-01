# Author: Joe Perez-Rogers
# Date: 10/29/2014
# Purpose: Script to run linear models and return a series of metrics in a script

fitLM <- function(eset,mod,term=ncol(mod)){
  fit <- lmFit(exprs(eset),design=mod)
  fit2 <- ebayes(fit)
  p.cuts <- c(sum(fit2$p.value[,term]<0.05),sum(fit2$p.value[,term]<0.005))
  fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
  fit.t.sorted <- sort(fit.t[,term],decreasing=T)
  p.adj <- p.adjust(fit2$p.value[,term],"fdr")
  fdr.cuts <- c(sum(p.adj<0.25),sum(p.adj<0.05))
  return(list("fit2"=fit2,"fit.t"=fit.t,"fit.t.sorted"=fit.t.sorted,"p.cuts"=p.cuts,"fdr.cuts"=fdr.cuts))
}