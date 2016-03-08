geneset.cor.test <- function(geneset1, geneset2, eset, type=c("gsva", "assign", "pca"), num.cors=100){
  require(GSVA)
  #require(ASSIGN)
  genesets <- list(set1=geneset1, set2=geneset2)
  # run the initial gsva analysis
  true.gsvaResults <- gsva(eset, genesets)
  #gsvaResults <- gsva(eset, genesets, no.bootstraps=10, method="zscore")
  # calculate the true correlation between geneset GSVA scores
  #true.cor <- as.numeric(cor(t(exprs(true.gsvaResults$es.obs[1, ])), t(exprs(true.gsvaResults$es.obs[2, ]))))
  true.cor.test <- cor.test(t(exprs(true.gsvaResults$es.obs[1, ])), t(exprs(true.gsvaResults$es.obs[2, ])))
  true.cor <- true.cor.test$estimate
  N <- sampleNumber(eset)
  #true.t.cor <- true.cor/(sqrt(1-(true.cor^2))/(N-2))
  true.t.cor <- true.cor.test$statistic
    
  false.cors <- matrix(nrow=num.cors,ncol=2)
  for(i in 1:num.cors){
    # generate the false gensets
    false.gns <- generate.random.genesets(geneset1, geneset2, eset)
    #false.geneset1 <- false.gns[[1]]; false.geneset2 <- false.gns[[2]]
    # calculate the GSVA scores
    false.gsvaResults <- gsva(eset, false.gns)
    tempCor <- cor.test(t(exprs(false.gsvaResults$es.obs[1, ])), t(exprs(false.gsvaResults$es.obs[2, ])))
    #rho <- cor(t(exprs(false.gsvaResults$es.obs[1, ])), t(exprs(false.gsvaResults$es.obs[2, ])))
    #false.cors[i, 1] <- rho
    false.cors[i, 1] <- tempCor$estimate
    #false.t.cor <- rho/(sqrt(1-(rho^2))/(N-2))
    #false.cors[i, 2] <- false.t.cor
    false.cors[i, 2] <- tempCor$statistic
  }
  
  print("out of the loop")
  if(true.cor >= 0){
    print("in the first if")
    print(false.cors)
    print(true.t.cor)
    print(num.cors)
    print("try looking at column2 only")
    print(false.cors[,2])
    print("calc temp1")
    temp1 <- length(which(false.cors[, 2] > true.t.cor))
    print("calc p")
    p <- temp1/num.cors
    #p <- sum(false.cors[, 2] > true.t.cor)/num.cors
    print("end1")
  } else{
    print("in the first else")
    p <- sum(false.cors[, 2] < true.t.cor)/num.cors
    print("end2")
  }

  print("out of the if-else")
  return(list(true.cor=true.cor, p=p, false.cors=false.cors))
  # the generation of false correlations is complete; now calculate
  # the t-distribution following statistic for each of the correlations
  # t = r/sqrt(1-r^2)/(N-2)
  
  # use this distribution to generate an empiric p-value for true.cor
  # what proportion of values (of the default=100) are above or below 
  # the actual true.cor value
}

####
# Generate false genesets
generate.random.genesets <- function(geneset1, geneset2, eset){
  
  # calculate the number of overlapping genes, num.Overlapping
  # this will heavily depend on how the gensets are submitted - now just limit to personal use
  num.overlapping = sum(as.character(geneset1) %in% as.character(geneset2))
  size.geneset1 <- length(geneset1); size.geneset2 <- length(geneset2)
  
  # define the features to select from for the genesets
  features <- featureNames(eset)
  
  if(num.overlapping >= 1){
    # select a random of number of genes from eset equal to the number of overlapping genes
    false.geneset.overlap <- sample(features, size=num.overlapping, replace=FALSE)
    
    # remove the featureNames chosen for the overlapping set from the remainder of names to choose
    # there will be no duplicate genes within one geneset 
    features <- setdiff(features, false.geneset.overlap)
    
    # fill out the remainder of the two genesets
    false.geneset1 <- union(sample(features, size=size.geneset1-num.overlapping, replace=FALSE),
                            false.geneset.overlap)
    
    false.geneset2 <- union(sample(features, size=size.geneset2-num.overlapping, replace=FALSE),
                            false.geneset.overlap)
  }  
  else {
    false.geneset1 <- sample(features, size=size.geneset1, replace=FALSE)
    false.geneset2 <- sample(features, size=size.geneset2, replace=FALSE)
  }
  
  return(list(false.geneset1=false.geneset1, false.geneset2=false.geneset2))
  
  
}