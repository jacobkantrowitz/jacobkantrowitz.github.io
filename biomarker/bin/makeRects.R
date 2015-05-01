# Script to put boxes around cells in heatmap

makeRects <- function(tfMat,border){
  cAbove = expand.grid(1:ncol(tfMat),1:nrow(tfMat))[tfMat,]
  xl=cAbove[,1]-0.50
  yb=cAbove[,2]-0.50
  xr=cAbove[,1]+0.50
  yt=cAbove[,2]+0.50
  rect(xl,yb,xr,yt,border=border,lwd=1)
}