# removeBioReps takes an expressionSet object and removes biological replicates from it
# by removing those CEL names that overlap
# removeBioReps returns the 'clean' duplicate-free expressionSet

removeBioReps <- function(exprData)
{
  celNames <- apply(X=as.data.frame(sampleNames(exprData)), 1, strsplit, split="_")
  temp <- character(length(sampleNumber(exprData)))
  for(i in 1:length(celNames)){ temp[i] <- celNames[[i]]$'sampleNames(exprData)'[5]}
  inds <- match(unique(temp), temp)
  exprData <- exprData[,inds]
  
  return(exprData)
}