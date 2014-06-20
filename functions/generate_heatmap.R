generate_heatmap  <- function(inds, exprData, rowClusters=NULL, mthd="average", colClus=exprData$indicator,mn="Figure", tp="CC")
{
  if (tp == "indicator") {
    clabels <- cbind("Indicator" = copdca_colors[exprData$indicator],
                     "Smoking Status" = smoking_colors[exprData$SMKc],
                     "Gender" = gender_colors[exprData$GENDERc],
                     "Batch" = batch_colors[exprData$BATCH])
    colClus <- exprData$indicator
  }
  
  else if (tp == "COPD2_R7") {
    clabels <- cbind("COPD2_R7" = copd_colors[exprData$COPD2_R7],
                     "Smoking Status" = smoking_colors[exprData$SMKc],
                     "Gender" = gender_colors[exprData$GENDERc],
                     "Batch" = batch_colors[exprData$BATCH])
    colClus <- exprData$COPD2_R7
  }
  
  else if (tp == "COPD2LLN") {
    clabels <- cbind("COPD2LLN" = copd_colors[exprData$COPD2LLN],
                     "Smoking Status" = smoking_colors[exprData$SMKc],
                     "Gender" = gender_colors[exprData$GENDERc],
                     "Batch" = batch_colors[exprData$BATCH])
    colClus <- exprData$COPD2LLN
  }
  
  else if (tp == "FinalCaDXc"){
    clabels <- cbind("FinalCaDXc" = cancer_colors[exprData$FinalCaDXc],
                     "Smoking Status" = smoking_colors[exprData$SMKc],
                     "Gender" = gender_colors[exprData$GENDERc],
                     "Batch" = batch_colors[exprData$BATCH])
    colClus <- exprData$FinalCaDXc
  }
 

  
  #heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
  heatmap3(exprData[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
  
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