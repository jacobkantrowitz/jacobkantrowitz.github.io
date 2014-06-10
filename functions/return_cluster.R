return_cluster <- function(inds, exprData, n.clusters=2, type=ROWS, mthd="average", colClus=exprData$indicator,mn="Figure")
{
  # code to define clusters from heatmap dendrogram
  clabels <- cbind(
    "Indicator" = copdca_colors[exprData$indicator],
    "COPD Status" = copd_colors[exprData$COPD],
    "Smoking Status" = smoking_colors[exprData$SMK],
    "Cancer Status" = cancer_colors[exprData$CANCER],
    "Gender" = gender_colors[exprData$GENDER],
    "Batch" = batch_colors[exprData$BATCH]
  )
  
  pdf(NULL)
  # 1. Main - uses ward clustering; semi-supervised
  #res <- heatmap3(data_to_analyze[inds,], col = bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method="ward"), col.clustering = "semisupervised", ColSideColors = clabels, main = "Figure")
  
  # 2. uses average clustering; semi-supervised
  res <- heatmap3(exprs(exprData)[inds,], col = bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
  
  # 3. unsupervised
  #res <- heatmap3(data_to_analyze[inds,], col = bluered, keep.dendro=TRUE, col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")	
  dev.off()
  
  # For columns, i.e. by patients
  if(type==COLUMNS)
  {
    clusters <- cutree(as.hclust(res$Colv), k=n.clusters)
  }
  # For rows, i.e. by genes
  if(type==ROWS)
  {
    clusters <- cutree(as.hclust(res$Rowv), k=n.clusters)
  }
  
  return(clusters)
}