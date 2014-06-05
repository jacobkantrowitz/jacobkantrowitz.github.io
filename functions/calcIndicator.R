# calcIndicator takes an expressionSet object and two factors each with 2 levels 
# within the phenoytpe data
# calcIndicator returns an expressionSet object with a new factor, indicator, that
# indicates group membership in the 2x2 table that results from the interesection of
# the two factors
calcIndicator <- function(exprData, factor1, factor2="FinalCaDXc")
{
  ind1 <- which(exprData[[factor1]]==0 & exprData[[factor2]]==0)
  ind2 <- which(exprData[[factor1]]==1 & exprData[[factor2]]==0)
  ind3 <- which(exprData[[factor1]]==0 & exprData[[factor2]]==1)
  ind4 <- which(exprData[[factor1]]==1 & exprData[[factor2]]==1)
  
  indicator <- numeric(sampleNumber(exprData))
  indicator[ind1] <- 1
  indicator[ind2] <- 2
  indicator[ind3] <- 3
  indicator[ind4] <- 4
  indicator <- factor(indicator)
  exprData$indicator <- indicator
  return(exprData)
}