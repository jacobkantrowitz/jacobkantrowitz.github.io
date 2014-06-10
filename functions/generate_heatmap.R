generate_heatmap  <- function(inds, exprData, rowClusters=NULL, mthd="average", colClus=phen$indicator,mn="Figure")
{
  clabels <- cbind(
    "Indicator" = copdca_colors[exprData$indicator],
    "COPD Status" = copd_colors[exprData$COPD],
    "Smoking Status" = smoking_colors[exprData$SMK],
    "Cancer Status" = cancer_colors[exprData$CANCER],
    "Gender" = gender_colors[exprData$GENDER],
    "Batch" = batch_colors[exprData$BATCH]
  )
  

  
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
  heatmap3(exprData[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabelsNEW, main = mn)
  
  if(!is.null(rowClusters))
  {
    heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, RowSideColors = rowClusters, main = mn)
  }
  
  # Uses ward for clustering; semi-supervised
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="ward"), col.clustering = "semisupervised", ColSideColors = clabelsA, main = "Figure")
  
  # MAIN - Uses average for clustering; semi-supervised
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "semisupervised", ColSideColors = clabels, main = "Figure")
  
  # Uses average for clustering; unsupervised
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")
  
  # Uses centroid for clustering; unsupervised
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="centroid"), col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")
  
  # Includes row labels when RowClusters is used; semi-supervised
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "semisupervised", ColSideColors = clabels, RowSideColors = rowClusters, main = "Figure")
}