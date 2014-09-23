lmFitWrapper <- function(eset, covariates, varOfInterest, adjust.method="fdr", p.value=0.05, name="Analysis"){
  
  cat(name, "\n")
  toReturn <- list()
  
  # check that all covariates are in eset
  if(sum(is.na(match(covariates, varLabels(eset)))) == 0){
    
    # construct model
    model <- c("~1", covariates)
    model <- formula(paste(c("~1", covariates), collapse="+"))
    design <- model.matrix(model, data=eset)
    
    # run lmFit
    fit <- lmFit(eset, design)
    # run eBayes
    fit <- eBayes(fit)
    
    # determine significant results, correcting if called for
    results <- decideTests(fit, adjust.method=adjust.method, p.value=p.value)
    print(summary(results))
    
    # find the indices and names of the genes of interest
    inds <- which(results[, varOfInterest+1] != 0)
    geneSymbols <- fit$genes$Symbol[inds]
    
    cat("\nNumber of significant genes for", covariates[varOfInterest], ": ", length(inds), "\n")
    
    toReturn$name <- name
    toReturn$hMapName <- cat(name, "\nAdjust:", adjust.method, p.value)
    toReturn$model <- model
    toReturn$design <- design
    toReturn$fit <- fit
    toReturn$results <- results
    toReturn$inds <- inds
    toReturn$geneSymbols <- geneSymbols
    toReturn$eset <- eset
    toReturn$adjust.method <- adjust.method
    toReturn$p.value <- p.value
    toReturn$varOfInterest <- covariates[varOfInterest]
    
    return(toReturn)
  }
  else{
    print("ERROR: Not all covariates are in eset")
  }
  
  
  
  
}