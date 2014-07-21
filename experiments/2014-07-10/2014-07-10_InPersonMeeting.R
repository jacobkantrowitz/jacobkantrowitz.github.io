
## ----GlobalVariables, echo=FALSE-----------------------------------------
cacheOption = FALSE



## ----Setup, cache=cacheOption, echo=FALSE, warning=FALSE, results="hide", message=FALSE----
# set the working directory to the directory of this script
setwd("/restricted/projectnb/pulmarray/LinGA_protected/Allegro/COPD_Cancer/experiments/2014-07-13")

# load the data using a script that removes patients with:
#  Cancer = NA
#  SMK    = 3
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/scripts/AllegroSetup.R")
library(znorm)

# remove the patients who have a DK for Cancer status as this is the a
# phenotype of interest (along with COPD)
eset <- removeFactorLevel(eset, "FinalCaDXc", "DK")

# remove the biological replicates
eset <- removeBioReps(eset)

# create a subset of all those subjects with PFTs and remove any NA/DK values in necessary fields
esetPFT <- cleanNAForAnalysis(eset, "COPD2_R7")
esetPFT <- removeFactorLevel(esetPFT, "COPD2_R7", "DK")
esetPFT <- calcIndicator(esetPFT, "FinalCaDXc", "COPD2_R7")

iterSubs <- function(eset){
  esetPFTFilterSub <- eset
  
  disFreeCurrent <- which(esetPFTFilterSub$SMKc==1 & esetPFTFilterSub$indicator==1)
  disFreeFormer <- which(esetPFTFilterSub$SMKc==2 & esetPFTFilterSub$indicator==1)
  disFreeFormer <- sample(disFreeFormer, size=length(disFreeCurrent), replace=FALSE)
  
  cancerCurrent <- which(esetPFTFilterSub$SMKc==1 & esetPFTFilterSub$indicator==2)
  cancerFormer <- which(esetPFTFilterSub$SMKc==2 & esetPFTFilterSub$indicator==2)
  cancerFormer <- sample(cancerFormer, size=length(cancerCurrent), replace=FALSE)
  
  copdCurrent <- which(esetPFTFilterSub$SMKc==1 & esetPFTFilterSub$indicator==3)
  copdFormer <- which(esetPFTFilterSub$SMKc==2 & esetPFTFilterSub$indicator==3)
  copdCurrent <- sample(copdCurrent, size=length(copdFormer), replace=FALSE)
  
  bothCurrent <- which(esetPFTFilterSub$SMKc==1 & esetPFTFilterSub$indicator==4)
  bothFormer <- which(esetPFTFilterSub$SMKc==2 & esetPFTFilterSub$indicator==4)
  bothFormer <- sample(bothFormer, size=length(bothCurrent), replace=FALSE)
  
  subSample <- c(disFreeCurrent, disFreeFormer, cancerCurrent, cancerFormer,
                 copdCurrent, copdFormer, bothCurrent, bothFormer)
  
  esetPFTFilterSub <- esetPFTFilterSub[, subSample]
  
  # perform signal/noise filtering (as above) based on median
  md <- median(exprs(esetPFTFilterSub))
  passFilter <- logical(featureNumber(esetPFTFilterSub))
  for(i in 1:featureNumber(esetPFTFilterSub)) {
    passFilter[i] <- sum(exprs(esetPFTFilterSub)[i, ] > md) > 0
  }
  
  # remove those genes not passing the median filter
  esetPFTSub <- esetPFTFilterSub
  esetPFTFilterSub <- esetPFTFilterSub[passFilter, ]
  
  ## ----goal3aInteractionModel----------------------------------------------
  designIntrx <- model.matrix(~ 1 + AGEcalc + SMKc + COPD2_R7*FinalCaDXc, data=esetPFTFilterSub)
  
  fitIntrx <- lmFit(esetPFTFilterSub, designIntrx)
  fitIntrx <- eBayes(fitIntrx)
  
  resultsIntrx <- decideTests(fitIntrx, adjust.method="fdr", p.value=0.05)
  resultsIntrx2 <- decideTests(fitIntrx, adjust.method="fdr", p.value=0.25)
  summary(resultsIntrx)
  
  returnResults <- list()
  returnResults[[1]] <- sum(resultsIntrx[, 6] != 0)
  returnResults[[2]] <- esetPFTFilterSub$BARCODE
  returnResults[[3]] <- featureNames(esetPFTFilterSub)[resultsIntrx[, 6] != 0]
  returnResults[[4]] <- sum(resultsIntrx2[, 6] != 0)
  
  return(returnResults)
}

#generate_heatmap(which(resultsIntrx[, 6] != 0), esetPFTFilterSub, tp="indicator")

iterResults <- list()
for(i in 1:1000){
  print(i)
  iterResults[[i]] <- iterSubs(esetPFT)
}

iterNum <- numeric(1000)
iterNum2 <- numeric(1000)
iterBarcode <- matrix(nrow=1000, ncol=308)
iterGenes <- list()

for(i in 1:1000){
  iterNum[i] <- iterResults[[i]][[1]]
  iterNum2[i] <- iterResults[[i]][[4]]
  iterBarcode[i, ] <- iterResults[[i]][[2]]
  iterGenes[[i]] <- iterResults[[i]][[3]]
}

pdf("numSigGenesIter.pdf")
hist(iterNum[iterNum>0])
dev.off()

# produce maps representing the 10 deciles of the non-zero significant lists
o <- sort(iterNum, decreasing=TRUE, index.return=TRUE)
totalG0 <- sum(iterNum>0)
decileInds <- round(seq(from=1, to=totalG0, by=totalG0/10))
for(i in 1:length(decileInds)){
  # find the relevant genes
  genes <- iterGenes[[o$ix[decileInds[i]]]]
  
  # find the subset of barcodes included in the analysis
  barcodes <- unique(iterBarcode[o$ix[decileInds[i]], ])
  tempEset <- esetPFT[, match(barcodes, esetPFT$BARCODE)]
  
  #create smoking residuals
  designSMK <- model.matrix(~1 + SMKc, data=tempEset)
  fitSMK <- lmFit(tempEset, designSMK)
  fitSMK <- eBayes(fitSMK)
  residSMK <- residuals(fitSMK, tempEset)
  exprs(tempEset) <- residSMK
  geneInds <- match(genes, rownames(exprs(esetPFT)))
  
  pm <- paste("Decile", i, "-", length(geneInds), "GenesSigFDR05SMKresid", sep="")
  # generate a heatmap for the genes in the subset
  pdf(paste(pm, ".pdf", sep=""))
  generate_heatmap(geneInds, tempEset, tp="indicator", mn=pm)
  dev.off()
  
}


# find the genes that appear the most in the iterGenes lists
# keep a data frame with names of genes and counts of appearances in iterGenes
# there are 19684 genes possible to appear

geneAppCounts <- matrix(data=numeric(19684), nrow=19684, ncol=1)
rownames(geneAppCounts) <- rownames(esetPFT)
for(i in 1:1000){
  addInd <- match(iterGenes[[i]], rownames(geneAppCounts))  
  geneAppCounts[addInd] <- geneAppCounts[addInd] + 1
}

# find the genes with at least a few appearances
genes1at <- rownames(geneAppCounts)[which(geneAppCounts>1)]
genes1 <- fitIntrx$genes[match(genes1at, rownames(fitIntrx$genes)),1]
genes2at <- rownames(geneAppCounts)[which(geneAppCounts>2)]
genes2 <- fitIntrx$genes[match(genes2at, rownames(fitIntrx$genes)),1]

genes3 <- rownames(geneAppCounts)[which(geneAppCounts>3)]
genes4 <- rownames(geneAppCounts)[which(geneAppCounts>4)]
genes5 <- rownames(geneAppCounts)[which(geneAppCounts>5)]
genes6 <- rownames(geneAppCounts)[which(geneAppCounts>6)]
genes7 <- rownames(geneAppCounts)[which(geneAppCounts>7)]
genes8at <- rownames(geneAppCounts)[which(geneAppCounts>8)]
genes8 <- fitIntrx$genes[match(genes8at, rownames(fitIntrx$genes)),1]

genes20at <- rownames(geneAppCounts)[which(geneAppCounts>20)]
genes20 <- fitIntrx$genes[match(genes20at, rownames(fitIntrx$genes)),1]
hmInds <- match(genes20at, rownames(esetPFTFilterSub))
generate_heatmap(hmInds, esetPFTFilterSub, tp="indicator")

#output the gene symbols
genes20df <- data.frame(genes20)
colnames(genes20df) <- c("GeneSymbol")

write.table(genes20df, file="genes20GeneSymbols", sep=",", row.names=FALSE, quote=FALSE)

genes40at <- rownames(geneAppCounts)[which(geneAppCounts>40)]
genes40 <- fitIntrx$genes[match(genes40at, rownames(fitIntrx$genes)),1]
hmInds2 <- match(genes40at, rownames(esetPFTFilterSub))
generate_heatmap(hmInds2, esetPFTFilterSub, tp="indicator")
genes40df <- data.frame(genes40)
colnames(genes40df) <- c("GeneSymbol")
write.table(genes40df, file="genes40GeneSymbols", sep=",", row.names=FALSE, quote=FALSE)