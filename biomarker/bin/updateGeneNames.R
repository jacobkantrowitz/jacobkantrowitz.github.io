updateGeneNames <- function(x){
  fData(x)$Symbol[which(is.na(fData(x)$Symbol))] <- paste("NA_Symbol",1:length(which(is.na(fData(x)$Symbol))),sep="")
  row.names(exprs(x)) <- fData(x)$Symbol
  return(x)
}