tgtGSVA <- function(genesets, target, toPlot=FALSE,
                    mn="GSVA Correlation Plots",
                    variables, correction="none",
                    pValue=0.05, varOfInterest=1){
  require(GSVA)
  target <- target[!is.na(fData(target)$Symbol), ]
  featureNames(target) <- fData(target)$Symbol
  gsvaResults <- gsva(target, genesets)
  
  toReturn <- list()
  # Plot the results
#  if(toPlot==TRUE){
#    tempCorPlot <- gsvaCorPlot(gsvaResults, mn=mn)
#    toReturn$gsvaCorPlot <- tempCorPlot
#  }
  
  # Test the AEGIS I GSVA results
  gsvaResultsLM <- lmFitWrapper(gsvaResults$es.obs,covariates=variables,
                                varOfInterest=varOfInterest,
                                adjust.method=correction, p.value=pValue)
  
  toReturn$gsvaResults <- gsvaResults
  toReturn$gsvaLM <- gsvaResultsLM
  return(toReturn)
  
  
}