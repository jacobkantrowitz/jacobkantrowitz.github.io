computeGenomicCorrelates <- function(eset){

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
  
  return(list(gender=GG,smoking=GS,py=GPY,age=GA))

}