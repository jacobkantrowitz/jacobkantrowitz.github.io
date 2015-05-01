lmFitWrapper <- function(eset, covariates, varOfInterest, adjust.method="fdr", p.value=0.05, name="Analysis"){
  
  toReturn <- list()
    
  # check that all covariates are in eset
  test <- check_covariates(eset, covariates)
  if(test == 0){
    
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
    hMapName <- paste(name, "\nAdjust:", adjust.method, p.value)
    cat(hMapName, "\n")
    print(summary(results))
    
    # find the indices and names of the genes of interest
    inds <- which(results[, varOfInterest+1] != 0)
    geneSymbols <- fit$genes$Symbol[inds]
    
    cat("\nNumber of significant genes for", covariates[varOfInterest], ": ", length(inds), "\n")
    if(adjust.method=="none"){
      cat("Number of significant genes expected by chance :", featureNumber(eset)*p.value)
    }
    
    toReturn$name <- name
    toReturn$hMapName <- hMapName
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
    toReturn$covariates <- covariates
    toReturn$nVarOfInterest <- varOfInterest
    
    return(toReturn)
  }
  else{
    print("ERROR: Not all covariates are in eset")
  } 
}

recalcModelResults <- function(model, p.value, adjust.method="fdr"){
  
  results <- decideTests(model$fit, adjust.method=adjust.method, p.value=p.value)
  hMapName <- paste(model$name, "\nAdjust:", adjust.method, p.value)
  cat(hMapName, "\n")
  print(summary(results))
  
  # find the indices and names of the genes of interest
  inds <- which(results[, model$nVarOfInterest+1] != 0)
  geneSymbols <- model$fit$genes$Symbol[inds]
  
  cat("\nNumber of significant genes for", model$varOfInterest, ": ", length(inds), "\n")
  if(adjust.method=="none"){
    cat("Number of significant genes expected by chance :", featureNumber(model$eset)*p.value)
  }
  
  model$hMapName <- hMapName
  model$results <- results
  model$inds <- inds
  model$geneSymbols <- geneSymbols
  model$adjust.method <- adjust.method
  model$p.value <- p.value
  
  return(model)
  
}