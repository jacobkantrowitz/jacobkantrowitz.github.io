# Author: Joe Perez-Rogers
# Date: 2014-09-03
# Purpose: A function to implement and evaluate the BronchoGen algorithm


evalBronchoGen <- function(eset,class,verbose=TRUE,file=NULL,plotit=FALSE,title=NULL){
  
  require(AUC)
  #   # Preprocessing
  #   ## Removing any samples with an NA for PY
  #   remove <- which(is.na(eset$PY))
  #   if(length(remove>0)){
  #     eset <- eset[,-remove]
  #     if(verbose){
  #       cat(paste("Removed ",length(remove)," sample(s) due to NAs in the PY variable\n",sep=""))
  #       }
  #     }
  #   
  #   ## Setting pack-years to a binary variable
  #   eset$PY[eset$PY<=10] <- 0
  #   eset$PY[eset$PY>10] <- 1
  
  ## Gender Correlate
  gender.gene <- "RPS4Y1"
  GG <- rep(0,ncol(exprs(eset)))
  GG[which(exprs(eset)[gender.gene,]>7.5)] <- 1
  
  ## Smoking Correlate
  intercept <- NA
  weights <- NA
  x<- NA
  intercept <- 40.8579
  weights <- as.vector(c(-0.4462,-2.1298,-1.8256))
  smk.genes <- c("SLC7A11","CLDN10","TKT")
  x <- intercept+weights%*%as.matrix(exprs(eset)[smk.genes,])
  GS <- exp(x)/(1+exp(x))
  
  ## PY Correlate
  intercept <- NA
  weights <- NA
  x<- NA
  intercept <- -5.1429
  weights <- as.vector(c(2.1891,-0.9506))
  py.genes <- c("RUNX1T1","AKR1C2")
  x <- intercept+weights%*%as.matrix(exprs(eset)[py.genes,])
  GPY <- exp(x)/(1+exp(x))
  
  ## Age Correlate
  intercept <- NA
  weights <- NA
  x<- NA
  intercept <- 62.4355
  weights <- as.vector(c(2.3201,1.7632,-1.6377,1.8109,-2.6387,2.1398,1.8630,1.2660,-1.3723,-2.6586,-1.2845,-1.3436,2.1489,1.9108,1.3775,1.4582,-1.8209,-1.0857,1.2493))
  age.genes <- c("CD52","SYT8","TNNT3","ALX1","KLRC4-KLRK1","RASA3","CERS3","ASPA","GRP","APOC1","EPHX3","REEP1","FAM198B","PCDHB4","PCDHB16","FOXD1","SPARC","NKAPL","GPR110")
  # Using "KLRC4-KLRK1" instead of "KLRK1" since "KLRK1" doesn't appear to be in my file
  GA <- intercept+weights%*%as.matrix(exprs(eset)[age.genes,])
  
  ## Cluster Means
  C1A.genes <- as.vector(c("BST1","CD177"))
  C1A <- colMeans(exprs(eset)[C1A.genes,])
  
  C1B.genes <- as.vector(c("ATP12A","TSPAN2"))
  C1B <- colMeans(exprs(eset)[C1B.genes,])
  
  C2.genes <- as.vector(c("GABBR1","MCAM","NOVA1","SDC2"))
  C2 <- colMeans(exprs(eset)[C2.genes,])
  
  C3.genes <- as.vector(c("CDR1","CGREF1","CLDN22","NKX3-1"))
  C3 <- colMeans(exprs(eset)[C3.genes,])
  
  C4A.genes <- as.vector(c("EPHX3","LYPD2"))
  C4A <- colMeans(exprs(eset)[C4A.genes,])
  
  C4B.genes <- as.vector(c("RNF150"))
  C4B <- exprs(eset)[C4B.genes,]
  
  cf <- c(3.317,0.062,0.545,0.166,-3.020,-0.441,-0.340,0.173,0.567,-0.316,-0.379)
  mat <- cbind(1,eset$AGE,GG,as.vector(GS),as.vector(GPY),C1A,C1B,C2,C3,C4A,C4B)
  
  result <- mat%*%cf
  pred.score <- exp(result)/(1+exp(result))
  roc.curve <- AUC::roc(predictions=as.numeric(pred.score),labels=as.factor(class))
  auc.value <- AUC::auc(roc.curve)
  
  # Plot it?
  if(plotit){
    plot(roc.curve,main=paste(title,"\nAUC=",round(auc.value,digits=2),sep=""))
  }
  if(!is.null(file)){
    pdf(file)
    plot(roc.curve,main=paste(round(auc.value,digits=2),sep=""))
    dev.off()
  }
    # Report AUC?
  if(verbose){
    cat(paste("AUC: ",auc.value,"\n",sep=""))
  }
  
  return(list("AUC"=auc.value,"ROC"=roc.curve,"scores"=pred.score))
}