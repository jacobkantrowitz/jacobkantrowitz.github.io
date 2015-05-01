# Author: Joseph Perez-Rogers
# Date: 2014-05-21
# Purpose: Script to plot heatmaps
# Usage:
# Input:
# Output:

makeHeatmap <- function(x,colors,method="semisupervised",rowClust="unsupervised",es=TRUE,plotit=FALSE,fname=NULL){
	require(heatmap3)
	brewer=rgb(colorRamp(c("blue", "white", "red"), space="rgb", interpolate="linear")(0:255/255), maxColorValue=255)
  if(es==TRUE){
    dat=exprs(x)
  } else {dat=x}
  if(plotit){
    if(is.null(fname)){
      pdf()
    } else
      pdf(fname)
  }
	res2 <- heatmap3(dat, scale = c("row"),
		 na.rm = TRUE, ColSideColors = colors,
		 hclustfun=function(d) hclust(d, method="ward.D"),
		 col.clustering = method,col = brewer,
		 labCol=NA,row.clustering=rowClust, keep.dendro=TRUE)
  if(plotit){
    dev.off()
  }
  return(res2)
}




