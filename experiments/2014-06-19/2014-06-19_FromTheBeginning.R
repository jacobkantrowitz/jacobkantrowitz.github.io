
## ----GlobalVariables-----------------------------------------------------
cacheOption = FALSE



## ----LoadData, echo=FALSE, cache=cacheOption-----------------------------
# set the working directory to the directory of this script
setwd("/restricted/projectnb/pulmarray/LinGA_protected/Allegro/COPD_Cancer/experiments/2014-06-19")

# load the data using a script that removes patients with:
#  Cancer = NA
#  SMK    = 3
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/scripts/AllegroSetup.R")

# remove the patients who have a DK for Cancer status as this is the a
# phenotype of interest (along with COPD)
eset <- removeFactorLevel(eset, "FinalCaDXc", "DK")



## ----QCData, fig.height=7, fig.width=8, fig.align='center', cache=cacheOption----

library(znorm)
pca <- prcomp(znorm(x=exprs(eset), margin=1), center=FALSE, scale.=FALSE)
plotTitle <- "PC1 v PC2"
plot(pca$rotation[,1], pca$rotation[,2], main=plotTitle, xlab="PC1", ylab="PC2")


## ----HighlightQCFigure1, fig.height=7, fig.width=8, fig.align='center', cache=cacheOption----

o <- order(pca$rotation[,1], decreasing=FALSE)
cols <- rep(1, sampleNumber(eset))
cols[o[1:4]] <- 2

plotTitle <- "PC1 v PC2: highlighting those samples to remove"
plot(pca$rotation[,1], pca$rotation[,2], col=cols, main=plotTitle,
     xlab="PC1", ylab="PC2")

plotTitle <- "PC1 v PC2: samples removed"
samplesToKeep <- o[5:sampleNumber(eset)]
plot(pca$rotation[samplesToKeep, 1], pca$rotation[samplesToKeep, 2], 
     main=plotTitle, xlab="PC1", ylab="PC2")

esetQC1 <- eset[,samplesToKeep]



## ----RemovingReplicates, fig.height=7, fig.width=8, fig.align='center', cache=cacheOption----
esetRepsRemoved <- removeBioReps(esetQC1)
pca2 <- prcomp(znorm(exprs(esetRepsRemoved), margin=1), center=FALSE, scale.=FALSE)

plotTitle <- "PC1 v PC2: eSet with Replicates Removed"
plot(pca2$rotation[,1], pca2$rotation[,2],
     main=plotTitle, xlab="PC1", ylab="PC2")


## ----HighlightQCFigure2, fig.height=7, fig.width=8, fig.align='center', cache=cacheOption----

# now highlight those for removal
o2 <- order(pca2$rotation[,1], decreasing=FALSE)
cols <- rep(1, sampleNumber(esetRepsRemoved))
cols[o2[c(1,692:695)]] <- 2

plotTitle <- "PC1 v PC2: highlighting those samples to remove"
plot(pca2$rotation[,1], pca2$rotation[,2], col=cols, main=plotTitle,
     xlab="PC1", ylab="PC2")

# now plot the PCA results having removed those subjects
plotTitle <- "PC1 v PC2: samples removed"
samplesToKeep <- o2[-c(1,692:695)]
plot(pca2$rotation[samplesToKeep, 1], pca2$rotation[samplesToKeep, 2], 
     main=plotTitle, xlab="PC1", ylab="PC2")

esetQC2 <- esetRepsRemoved[,samplesToKeep]



## ----ModelingCancerSignal, cache=cacheOption-----------------------------

# define an appropriate data set to work with
esetCa <- esetQC2

# define the model matrix starting with the simplest
# and exanding from there if appropriate or necessary
designCa1 <- model.matrix(~ 0 + FinalCaDXc, data=esetCa)
colnames(designCa1) <- c("noCancer", "Cancer")

# fit the design matrix to the data
fitCa1 <- lmFit(esetCa, designCa1)

# make contrasts to pull out cancer signal (can do this more simply)
contMatrixCa1 <- makeContrasts(HealthyVCancer = noCancer - Cancer,
                              levels=designCa1)

# calculate the contrasts
fitCa1_2 <- contrasts.fit(fitCa1, contMatrixCa1)
fitCa1_2 <- eBayes(fitCa1_2)

resultsCa1 <- decideTests(fitCa1_2, adjust.method="fdr", p.value=0.05)
summary(resultsCa1)

# to get gene symbols from eset # fData(eset)$Symbol
# create a .rnk file for use with GSEA
# pull the t-statistics and gene symbols from fitCa1_2

tsCa1 <- fitCa1_2$t
rownames(tsCa1) <- fitCa1_2$genes[,1]
write.table(tsCa1, file="cancerSignalTs.rnk", sep="\t", row.names=TRUE, quote=FALSE)



## ----ModelingCOPDSignal, cache=cacheOption-------------------------------

# define an appropriate data set to work with
esetCp <- cleanNAForAnalysis(esetQC2, "COPD2_R7")
esetCp <- removeFactorLevel(esetCp, "COPD2_R7", "DK")

# define the model matrix starting with the simplest
# and exanding from there if appropriate or necessary
designCp1 <- model.matrix(~ 0 + COPD2_R7, data=esetCp)
colnames(designCp1) <- c("noCOPD", "COPD")

# fit the design matrix to the data
fitCp1 <- lmFit(esetCp, designCp1)

# make contrasts to pull out cancer signal (can do this more simply)
contMatrixCp1 <- makeContrasts(HealthyVCOPD = noCOPD - COPD,
                              levels=designCp1)

# calculate the contrasts
fitCp1_2 <- contrasts.fit(fitCp1, contMatrixCp1)
fitCp1_2 <- eBayes(fitCp1_2)

resultsCp1 <- decideTests(fitCp1_2, adjust.method="fdr", p.value=0.05)
summary(resultsCp1)

# to get gene symbols from eset # fData(eset)$Symbol
# create a .rnk file for use with GSEA
# pull the t-statistics and gene symbols from fitCa1_2

tsCp1 <- fitCp1_2$t
rownames(tsCp1) <- fitCp1_2$genes[,1]
 write.table(tsCp1, file="COPDSignalTs.rnk", sep="\t", row.names=TRUE, quote=FALSE)




## ----Heatmaps, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----

# print a heatmap based on the COPD2_R7 model
generate_heatmap(which(p.adjust(fitCp1_2$p.value, method="fdr") < 0.0005), esetCp, tp="COPD2_R7")

# print a heatmap based on the FinalCaDXc model
generate_heatmap(which(p.adjust(fitCa1_2$p.value, method="fdr") < 0.02), esetCa, tp="FinalCaDXc")




## ----SMKResidualHeatmap, fig.height=7, fig.width=8, fig.align='center', cache=cacheOption----

# FOR THE COPD DATA
# first, generate a new model with SMKc as the only term
# maybe include RIN as well
designCp_SMKc <- model.matrix(~ 1 + SMKc, data=esetCp)
colnames(designCp_SMKc)[2] <- "Former"

# fit the design matrix to the data
fitCp_SMKc <- lmFit(esetCp, designCp_SMKc)
fitCp_SMKc <- eBayes(fitCp_SMKc)
resCp_SMKc <- residuals(fitCp_SMKc, esetCp)
esetCpRES <- esetCp
exprs(esetCpRES) <- resCp_SMKc

# print a heatmap based on the COPD2_R7 model
generate_heatmap(which(p.adjust(fitCp1_2$p.value, method="fdr") < 0.0005), esetCpRES, tp="COPD2_R7")

# FOR THE CANCER DATA
# first, generate a new model with SMKc as the only term
# maybe include RIN as well
designCa_SMKc <- model.matrix(~ 1 + SMKc, data=esetCa)
colnames(designCa_SMKc)[2] <- "Former"

# fit the design matrix to the data
fitCa_SMKc <- lmFit(esetCa, designCa_SMKc)
fitCa_SMKc <- eBayes(fitCa_SMKc)
resCa_SMKc <- residuals(fitCa_SMKc, esetCa)
esetCaRES <- esetCa
exprs(esetCaRES) <- resCa_SMKc

# print a heatmap based on the FinalCaDXc model
generate_heatmap(which(p.adjust(fitCa1_2$p.value, method="fdr") < 0.02), esetCaRES, tp="FinalCaDXc")



## ----ResidHeatmapsIndicatorCp, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----
# first need to generate the indicator factor
esetCpRESi <- calcIndicator(esetCpRES, "FinalCaDXc", "COPD2_R7")

# print a heatmap based on the COPD2_R7 model
generate_heatmap(which(p.adjust(fitCp1_2$p.value, method="fdr") < 0.0005), esetCpRESi, tp="indicator")


## ----ResidHeatmapsIndicatorCa, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----
# first need to remove the necessary samples
esetCaRESi <- cleanNAForAnalysis(esetCaRES, "COPD2_R7")
esetCaRESi <- removeFactorLevel(esetCaRESi, "COPD2_R7", "DK")

# now need to generate the indicator factor
esetCaRESi <- calcIndicator(esetCaRESi, "FinalCaDXc", "COPD2_R7")

# print a heatmap based on the FinalCaDXc model
generate_heatmap(which(p.adjust(fitCa1_2$p.value, method="fdr") < 0.02), esetCaRESi, tp="indicator")



## ----GenerateRankedLists, cache=cacheOption------------------------------
# preview the number of significant genes at FDR < 0.05
summary(decideTests(fitCa1_2, adjust.method="fdr", p.value=0.05))
summary(decideTests(fitCp1_2, adjust.method="fdr", p.value=0.05))
summary(decideTests(fitCa_SMKc, adjust.method="fdr", p.value=0.05))
summary(decideTests(fitCp_SMKc, adjust.method="fdr", p.value=0.05))






## ----GenerateGeneSets, cache=cacheOption---------------------------------




## ----ModelInteraction, cache=cacheOption---------------------------------




## ----HairballPCA, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----
# pull the indices of the genes with fdr < 0.05 in each model
#cafdr05 <- which(p.adjust(fitCa1_2$p.value, method="fdr") < 0.05)
#cpfdr05 <- which(p.adjust(fitCp1_2$p.value, method="fdr") < 0.05)
#smfdr05 <- which(p.adjust(fitCa_SMKc$p.value[,2], method="fdr") < 0.05)

# take 500 genes from each set
# for cancer
ca500ind <- order(p.adjust(fitCa1_2$p.value, method="fdr"),
                  decreasing=FALSE)[1:500]
esetCa500L <- esetCa[ca500ind, ]
esetCa500s <- esetCp[ca500ind, ]

# for COPD
cp500ind <- order(p.adjust(fitCp1_2$p.value, method="fdr"),
                  decreasing=FALSE)[1:500]
esetCp500L <- esetCa[cp500ind, ]
esetCp500s <- esetCp[cp500ind, ]

# for Smoking
sm500ind <- order(p.adjust(fitCa_SMKc$p.value[, 2], method="fdr"),
                  decreasing=FALSE)[1:500]
esetSmk500L <- esetCa[sm500ind, ]
esetSmk500s <- esetCp[sm500ind, ]

# thought is that these probably all need to be done in the same
# set of patients so that I can make predictions for each set from
# each other set (created L/s sets Large and small to handle this)

# run PCA on each of the 500 gene sets
# Large sets
pcaLCa <- prcomp(znorm(exprs(esetCa500L), margin=1), center=FALSE, scale.=FALSE)
pcaLCp <- prcomp(znorm(exprs(esetCp500L), margin=1), center=FALSE, scale.=FALSE)
pcaLSm <- prcomp(znorm(exprs(esetSmk500L), margin=1), center=FALSE, scale.=FALSE)

# small sets
pcasCa <- prcomp(znorm(exprs(esetCa500s), margin=1), center=FALSE, scale.=FALSE)
pcasCp <- prcomp(znorm(exprs(esetCp500s), margin=1), center=FALSE, scale.=FALSE)
pcasSm <- prcomp(znorm(exprs(esetSmk500s), margin=1), center=FALSE, scale.=FALSE)

# plots of PC1s from each of the sets against the other
plotTitle1L <- "Large Set, PC1 Smoking vs PC1 Cancer"
plot(pcaLSm$rotation[, 1], pcaLCa$rotation[, 1], main=plotTitle1L,
     xlab="PC1 Smoking", ylab="PC1 Cancer")
cor.test(pcaLSm$rotation[, 1], pcaLCa$rotation[, 1])

plotTitle2L <- "Large Set, PC1 Smoking vs PC1 COPD"
plot(pcaLSm$rotation[, 1], pcaLCp$rotation[, 1], main=plotTitle2L,
     xlab="PC1 Smoking", ylab="PC1 COPD")
cor.test(pcaLSm$rotation[, 1], pcaLCp$rotation[, 1])

plotTitle3L <- "Large Set, PC1 COPD vs PC1 Cancer"
plot(pcaLCp$rotation[, 1], pcaLCa$rotation[, 1], main=plotTitle3L,
     xlab="PC1 COPD", ylab="PC1 Cancer")
cor.test(pcaLCp$rotation[, 1], pcaLCa$rotation[, 1])

# the smaller sets
plotTitle1s <- "Small Set, PC1 Smoking vs PC1 Cancer"
plot(pcasSm$rotation[, 1], pcasCa$rotation[, 1], main=plotTitle1s,
     xlab="PC1 Smoking", ylab="PC1 Cancer")
cor.test(pcasSm$rotation[, 1], pcasCa$rotation[, 1])

plotTitle2s <- "Small Set, PC1 Smoking vs PC1 COPD"
plot(pcasSm$rotation[, 1], pcasCp$rotation[, 1], main=plotTitle2s,
     xlab="PC1 Smoking", ylab="PC1 COPD")
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 1])

plotTitle3s <- "Small Set, PC1 COPD vs PC1 Cancer"
plot(pcasCp$rotation[, 1], pcasCa$rotation[, 1], main=plotTitle3s,
     xlab="PC1 COPD", ylab="PC1 Cancer")
cor.test(pcasCp$rotation[, 1], pcasCa$rotation[, 1])


## ----HairballROCAUC, cache=cacheOption-----------------------------------

library(pROC)
# calculate and plot ROC AUCs
# auROC in limma? roc function in pROC package
# use PC1 from each category as the predictor and the true labels 
# as the responses

# predict COPD with COPD, Smoking, Cancer PC1
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasSm$rotation[, 1]), 
         print.auc=TRUE, main="COPD predicted by Smoking PC1")
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasCp$rotation[, 1]),
         print.auc=TRUE, main="COPD predicted by COPD PC1")
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasCa$rotation[, 1]),
         print.auc=TRUE, main="COPD predicted by Cancer PC1")

# predict Cancer with COPD, Smoking, Cancer PC1
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasCa$rotation[, 1]),
         print.auc=TRUE, main="Cancer predicted by Cancer PC1")
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasSm$rotation[, 1]),
         print.auc=TRUE, main="Cancer predicted by Smoking PC1")
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasCp$rotation[, 1]),
         print.auc=TRUE, main="Cancer predicted COPD PC1")

# predict Smoking with COPD, Smoking, Cancer PC1
plot(roc(response=esetSmk500s$SMKc, predictor=pcasCp$rotation[, 1]),
         print.auc=TRUE, main="Smoking predicted by COPD PC1")
plot(roc(response=esetSmk500s$SMKc, predictor=pcasCa$rotation[, 1]),
         print.auc=TRUE, main="Smoking predicted by Cancer PC1")
plot(roc(response=esetSmk500s$SMKc, predictor=pcasSm$rotation[, 1]),
         print.auc=TRUE, main="Smoking predicted by Smoking PC1")




## ----ModelingSmokingPlusDisease, cache=cacheOption-----------------------
esetSmDis <- esetCp

# define the model matrix including disease term and smoking
designSmCa <- model.matrix(~ 1 + SMKc + FinalCaDXc, data=esetSmDis)
colnames(designSmCa)[2:3] <- c("Former", "Cancer")

# fit the design matrix to the data
fitSmCa1 <- lmFit(esetSmDis, designSmCa)
fitSmCa1 <- eBayes(fitSmCa1)

resultsSmCa1 <- decideTests(fitSmCa1, adjust.method="fdr", p.value=0.1)
summary(resultsSmCa1)



## ----CorrsCOPDvSmkAndCavSmk, cache=cacheOption---------------------------
# correlations of smoking PCs with COPD PCs
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 1])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 2])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 3])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 4])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 5])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 6])
cor.test(pcasSm$rotation[, 1], pcasCp$rotation[, 7])

# correlations of smoking PCs with Cancer PCs
cor.test(pcasSm$rotation[, 1], pcasCa$rotation[, 1])
cor.test(pcasSm$rotation[, 1], pcasCa$rotation[, 2])
cor.test(pcasSm$rotation[, 1], pcasCa$rotation[, 3])
cor.test(pcasSm$rotation[, 1], pcasCa$rotation[, 4])


## ----uncorrPCsPredict, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----
# plots of PC1s from each of the sets against the other
# predict COPD with COPD, Smoking, Cancer PC1
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasSm$rotation[, 1]), 
         print.auc=TRUE, main="COPD predicted by Smoking PC1")
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasCp$rotation[, 7]),
         print.auc=TRUE, main="COPD predicted by COPD PC7")
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasCa$rotation[, 4]),
         print.auc=TRUE, main="COPD predicted by Cancer PC4")

# predict Cancer with COPD, Smoking, Cancer PC1
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasCa$rotation[, 4]),
         print.auc=TRUE, main="Cancer predicted by Cancer PC4")
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasSm$rotation[, 1]),
         print.auc=TRUE, main="Cancer predicted by Smoking PC1")
plot(roc(response=esetCa500s$FinalCaDXc, predictor=pcasCp$rotation[, 7]),
         print.auc=TRUE, main="Cancer predicted COPD PC7")

# predict Smoking with COPD, Smoking, Cancer PC1
plot(roc(response=esetSmk500s$SMKc, predictor=pcasCp$rotation[, 7]),
         print.auc=TRUE, main="Smoking predicted by COPD PC7")
plot(roc(response=esetSmk500s$SMKc, predictor=pcasCa$rotation[, 4]),
         print.auc=TRUE, main="Smoking predicted by Cancer PC4")
plot(roc(response=esetSmk500s$SMKc, predictor=pcasSm$rotation[, 1]),
         print.auc=TRUE, main="Smoking predicted by Smoking PC1")



## ----AdamSuggestions, fig.width=8, fig.height=7, fig.align='center', cache=cacheOption----
# run model with COPD*Cancer + Age and find genes associated with interaction
# use esetSmDis defined above as the same set of patients used for COPD analysis, i.e. those with PFT data available
# esetSmDis <- esetCp

# define the model matrix including disease term and smoking
designIntrx <- model.matrix(~ 1 + AGEcalc + COPD2_R7*FinalCaDXc, data=esetSmDis)
colnames(designIntrx)[2:5] <- c("Age", "COPD2_R7", "Cancer", "Interaction")

# fit the design matrix to the data
fitIntrx <- lmFit(esetSmDis, designIntrx)
fitIntrx <- eBayes(fitIntrx)

resultsIntrx <- decideTests(fitIntrx, adjust.method="fdr", p.value=0.05)
summary(resultsIntrx)

resultsIntrx <- decideTests(fitIntrx, adjust.method="fdr", p.value=0.25)
summary(resultsIntrx)

esetSmDisZ <- esetSmDis
# z-normalizing the expressiond data for the heatmap appears to make no change
exprs(esetSmDisZ) <- znorm(exprs(esetSmDisZ), margin=1)
generate_heatmap(which(resultsIntrx[, 5] != 0), esetSmDisZ, tp="indicator")


# adding smoking status to the model, which heavily influences the genes found above
designIntrx2 <- model.matrix(~ 1 + AGEcalc + SMKc + COPD2_R7*FinalCaDXc, data=esetSmDis)
colnames(designIntrx2)[2:6] <- c("Age", "Former", "COPD2_R7", "Cancer", "Interaction")

# fit the design matrix to the data
fitIntrx2 <- lmFit(esetSmDis, designIntrx2)
fitIntrx2 <- eBayes(fitIntrx2)

resultsIntrx2 <- decideTests(fitIntrx2, adjust.method="fdr", p.value=0.05)
summary(resultsIntrx2)

resultsIntrx2 <- decideTests(fitIntrx2, adjust.method="fdr", p.value=0.25)
summary(resultsIntrx2)

########
# running the simplest interaction model with no covariates
designIntrx3 <- model.matrix(~ 1 + COPD2_R7*FinalCaDXc, data=esetSmDis)
colnames(designIntrx3)[2:4] <- c("COPD2_R7", "Cancer", "Interaction")

# fit the design matrix to the data
fitIntrx3 <- lmFit(esetSmDis, designIntrx3)
fitIntrx3 <- eBayes(fitIntrx3)

resultsIntrx3 <- decideTests(fitIntrx3, adjust.method="fdr", p.value=0.05)
summary(resultsIntrx3)

resultsIntrx3 <- decideTests(fitIntrx3, adjust.method="fdr", p.value=0.25)
summary(resultsIntrx3)

# with no covariates, we are unable to find even 1 gene significant for the interacterm term 

# generate a ranked list for the interaction term
tsIntrx1 <- matrix(fitIntrx$t[, 5])
rownames(tsIntrx1) <- fitIntrx$genes[,1]
write.table(tsIntrx1, file="IntrxTs.rnk", sep="\t", row.names=TRUE, quote=FALSE)

# running PCA for genes associated with interaction
# testing them against smoking, cancer, copd, indicator
# for cancer
cc296ind <- order(p.adjust(fitIntrx$p.value, method="fdr"),
                  decreasing=FALSE)[1:500]
esetCC296s <- esetCp[cc296ind, ]

pcasCC <- prcomp(znorm(exprs(esetCC296s), margin=1), center=FALSE, scale.=FALSE)

# the smaller sets
plotTitleCC1 <- "Small Set, PC1 CC vs PC1 Smoking"
plot(pcasCC$rotation[, 1], pcasSm$rotation[, 1], main=plotTitleCC1,
     xlab="PC1 CC", ylab="PC1 Cancer")
cor.test(pcasCC$rotation[, 1], pcasSm$rotation[, 1])

plotTitleCC2 <- "Small Set, PC1 CC vs PC1 COPD"
plot(pcasCC$rotation[, 1], pcasCp$rotation[, 1], main=plotTitleCC2,
     xlab="PC1 CC", ylab="PC1 COPD")
cor.test(pcasCC$rotation[, 1], pcasCp$rotation[, 1])

plotTitleCC3 <- "Small Set, PC1 CC vs PC1 Cancer"
plot(pcasCC$rotation[, 1], pcasCa$rotation[, 1], main=plotTitleCC3,
     xlab="PC1 CC", ylab="PC1 Cancer")
cor.test(pcasCC$rotation[, 1], pcasCa$rotation[, 1])

esetCp500s <-calcIndicator(esetCp500s, "FinalCaDXc", "COPD2_R7")
# predict Smoking, COPD, Cancer, and indicator with interaction PC1
plot(roc(response=esetCp500s$SMKc, predictor=pcasCC$rotation[, 1]), 
         print.auc=TRUE, main="Smoking predicted by interaction PC1")
plot(roc(response=esetCp500s$COPD2_R7, predictor=pcasCC$rotation[, 1]),
         print.auc=TRUE, main="COPD predicted by interaction PC1")
plot(roc(response=esetCp500s$FinalCaDXc, predictor=pcasCC$rotation[, 1]),
         print.auc=TRUE, main="Cancer predicted by interaction PC1")
plot(roc(response=esetCp500s$indicator, predictor=pcasCC$rotation[, 1]),
         print.auc=TRUE, main="Indicator predicted by interaction PC1")

# generate a heatmap (it won't look good) of interaction term genes
esetCp <- calcIndicator(esetCp, "FinalCaDXc", "COPD2_R7")
generate_heatmap(which(resultsIntrx[,5] != 0), esetCp, tp="indicator")




## ----BoxplotsIntrx-------------------------------------------------------
pdf("IntrxBoxplots2")
for(i in 1:296){ 
  # generate the gene name
  ind <- which(resultsIntrx[,5] != 0)[i]
 
  gn <- paste(fitIntrx$genes[ind,1], ":", fitIntrx$genes[ind,2],  sep="")
 
  boxplot(exprs(esetCp[which(resultsIntrx[,5] != 0)[i], ])[1, ] ~ esetCp$indicator,
          names=c("Disease Free", "Cancer", "COPD", "Cancer+COPD"),
          col=c("green", "red", "blue", "yellow"),
          xlab="Group", ylab="Expression",
          main=gn)
  
 # legend("topleft", c("Disease Free", "Cancer", "COPD", "Cancer+COPD"), 
 #        fill=c("Green", "Red", "Blue", "Yellow"), inset=0.05)
  }

#

dev.off()


## ----MedianExpressionFiltering-------------------------------------------
# find those genes that do not have even 1 sample above the overall expressionSet median
md <- median(exprs(esetCp))
passFilter <- logical(featureNumber(esetCp))
for(i in 1:featureNumber(esetCp)) {
  passFilter[i] <- sum(exprs(esetCp)[i, ] > md) > 0
}

# remove those genes not passing the median filter
esetCpFilter <- esetCp[passFilter, ]

# rerun interaction analysis
fitIntrxF <- lmFit(esetCpFilter, designIntrx)
fitIntrxF <- eBayes(fitIntrx)

resultsIntrxF <- decideTests(fitIntrxF, adjust.method="fdr", p.value=0.25)
summary(resultsIntrxF)

# include SMKc just to see what happens
# define the model matrix including disease term and smoking
designIntrxFSmk <- model.matrix(~ 1 + AGEcalc + COPD2_R7*FinalCaDXc + SMKc,
                            data=esetCpFilter)
colnames(designIntrxFSmk)[2:5] <- c("Age", "COPD2_R7", "Cancer",
                                    "Former", "Interaction")

# fit the design matrix to the data
fitIntrxF2 <- lmFit(esetCpFilter, designIntrxFSmk)
fitIntrxF2 <- eBayes(fitIntrxF2)

resultsIntrxF2 <- decideTests(fitIntrxF2, adjust.method="fdr", p.value=0.05)
summary(resultsIntrxF2)

resultsIntrxF2 <- decideTests(fitIntrxF2, adjust.method="fdr", p.value=0.25)
summary(resultsIntrxF2)


# try including PC1 from COPD and Cancer (above) to correct for smoking-related
# disease signal

designIntrxFPC1 <- model.matrix(~1 + AGEcalc + pcasCp$rotation[, 1] + 
                                  pcasCa$rotation[,1] + COPD2_R7 + FinalCaDXc + 
                                  COPD2_R7:FinalCaDXc, data = esetCpFilter)

colnames(designIntrxFPC1)[2:7] <- c("Age", "COPD_PC1", "Cancer_PC1", "COPD2_R7",
                              "Cancer", "Interaction")

# fit the design matrix to the data
fitIntrxFPC1 <- lmFit(esetCpFilter, designIntrxFPC1)
fitIntrxFPC1 <- eBayes(fitIntrxFPC1)

resultsIntrxFPC1 <- decideTests(fitIntrxFPC1, adjust.method="fdr", p.value=0.05)
summary(resultsIntrxFPC1)


# try modeling with PC1s instead of the actual term (i.e. remove COPD2_R7, FinalCaDXc)
designIntrxFPCs <- model.matrix(~1 + AGEcalc + pcasCp$rotation[, 1]*pcasCa$rotation[,1],
                                data = esetCpFilter)

colnames(designIntrxFPCs)[2:5] <- c("Age", "COPD_PC1", "Cancer_PC1", "PCsInteraction")

# fit the design matrix to the data
fitIntrxFPCs <- lmFit(esetCpFilter, designIntrxFPCs)
fitIntrxFPCs <- eBayes(fitIntrxFPCs)

resultsIntrxFPCs <- decideTests(fitIntrxFPCs, adjust.method="fdr", p.value=0.01)
summary(resultsIntrxFPCs)



## ----KatieIdeasFromEmail-------------------------------------------------

# model to include RIN and to use FEV1% or FEV1/FVC instead of COPD binary
esetCpFilterPFT <- cleanNAForAnalysis(esetCpFilter, "RATIOc")
esetCpFilterPFT <- cleanNAForAnalysis(esetCpFilterPFT, "RIN")
designPFT1 <- model.matrix(~1 + SMKc + AGEcalc + FinalCaDXc + RATIOc + FinalCaDXc:RATIOc 
                           + SMKc:FinalCaDXc + SMKc:RATIOc
                           , data = esetCpFilterPFT)

colnames(designPFT1)[2:8] <- c("Former", "AGE", "Cancer", "Ratio", "DiseaseIntrx",
                               "CancerSmk", "COPDSmk")

# fit the design matrix to the data
fitPFT1 <- lmFit(esetCpFilterPFT, designPFT1)
fitPFT1 <- eBayes(fitPFT1)

resultsPFT1 <- decideTests(fitPFT1, adjust.method="fdr", p.value=0.05)
summary(resultsPFT1)

save_entrez(which(resultsPFT1[, 7] != 0), featureNames(esetCpFilterPFT), filename="genes717fdr05smokingcancerInteraction")

tsSMKCaIntrx <- matrix(fitPFT1$t[, 7])
rownames(tsSMKCaIntrx) <- fitPFT1$genes[,1]
write.table(tsSMKCaIntrx, file="SMKCancerTs.rnk", sep="\t", row.names=TRUE, quote=FALSE)

# try using COPD binary instead of continuous Ratio term
designPFT2 <- model.matrix(~1 + SMKc + AGEcalc + FinalCaDXc + COPD2_R7 + 
                             FinalCaDXc:COPD2_R7 + SMKc:FinalCaDXc + SMKc:COPD2_R7
                           , data = esetCpFilterPFT)

colnames(designPFT2)[2:8] <- c("Former", "AGE", "Cancer", "COPD", "DiseaseIntrx",
                               "CancerSmk", "COPDSmk")

# fit the design matrix to the data
fitPFT2 <- lmFit(esetCpFilterPFT, designPFT2)
fitPFT2 <- eBayes(fitPFT2)

resultsPFT2 <- decideTests(fitPFT2, adjust.method="fdr", p.value=0.05)
summary(resultsPFT2)

esetCpFilterPFT2 <- esetCpFilterPFT
esetCpFilterPFT2$SMKc <- as.factor(as.numeric(esetCpFilterPFT2$SMKc)-1)
esetCpFilterPFT2 <- calcIndicator(esetCpFilterPFT2, "SMKc", "FinalCaDXc")
generate_heatmap(which(resultsPFT2[, 7] != 0), esetCpFilterPFT2, tp="indicator")



## ----output296forEnrichr-------------------------------------------------
intrx296ind <- which(resultsIntrx[, 5] != 0)
intrx296 <- data.frame(fitIntrx$genes[intrx296ind,1])
colnames(intrx296) <- c("GeneSymbol")
intrx296$Group <- resultsIntrx[intrx296ind, 5]
intrx296$Group[intrx296$Group == -1] <- 0

write.table(intrx296, file="intrx296GeneSymbols", sep=",", row.names=FALSE, quote=FALSE)


## ----subsampling---------------------------------------------------------
 esetCpFilterPFT <- calcIndicator(esetCpFilterPFT, "FinalCaDXc", "COPD2_R7")
summary(esetCpFilterPFT$indicator[esetCpFilterPFT$SMKc==1])
summary(esetCpFilterPFT$indicator[esetCpFilterPFT$SMKc==2])
# we should sample to match the proportion of currents and formers in the COPD group
# here that ratio is 0.50
# sample the disease free group down to 0.50 - remove 26 samples from the former
# sample the cancer group down to 0.50 - remove 15 samples from the former
# sample the both group down to 0.50 - remove 2 samples from the former
# also check for more-or-less balance among other covariates (check new/old table 1)

esetPFTsub <- esetCpFilterPFT
disFreeCurrent <- which(esetPFTsub$SMKc==1 & esetPFTsub$indicator==1)
disFreeFormer <- which(esetPFTsub$SMKc==2 & esetPFTsub$indicator==1)
disFreeFormer <- sample(disFreeFormer, size=length(disFreeCurrent), replace=FALSE)

cancerCurrent <- which(esetPFTsub$SMKc==1 & esetPFTsub$indicator==2)
cancerFormer <- which(esetPFTsub$SMKc==2 & esetPFTsub$indicator==2)
cancerFormer <- sample(cancerFormer, size=length(cancerCurrent), replace=FALSE)

copdCurrent <- which(esetPFTsub$SMKc==1 & esetPFTsub$indicator==3)
copdFormer <- which(esetPFTsub$SMKc==2 & esetPFTsub$indicator==3)
copdFormer <- sample(copdFormer, size=length(copdCurrent), replace=FALSE)

bothCurrent <- which(esetPFTsub$SMKc==1 & esetPFTsub$indicator==4)
bothFormer <- which(esetPFTsub$SMKc==2 & esetPFTsub$indicator==4)
bothFormer <- sample(bothFormer, size=length(bothCurrent), replace=FALSE)

subSample <- c(disFreeCurrent, disFreeFormer, cancerCurrent, cancerFormer,
               copdCurrent, copdFormer, bothCurrent, bothFormer)

esetPFTsub <- esetPFTsub[, subSample]
summary(esetPFTsub$indicator[esetPFTsub$SMKc==1])
summary(esetPFTsub$indicator[esetPFTsub$SMKc==2])

# run design with the subsampled data and the median filter
designIntrxSub <- model.matrix(~ 1 + AGEcalc + COPD2_R7*FinalCaDXc, data=esetPFTsub)
colnames(designIntrxSub)[2:5] <- c("Age", "COPD2_R7", "Cancer", "Interaction")

fitIntrxSub <- lmFit(esetPFTsub, designIntrxSub)
fitIntrxSub <- eBayes(fitIntrxSub)

resultsIntrxSub <- decideTests(fitIntrxSub, adjust.method="fdr", p.value=0.25)
summary(resultsIntrxSub)

# run design with the median filter (not subsampled)
designIntrxFilt <- model.matrix(~ 1 + AGEcalc + COPD2_R7*FinalCaDXc, data=esetCpFilterPFT)
colnames(designIntrxFilt)[2:5] <- c("Age", "COPD2_R7", "Cancer", "Interaction")

fitIntrxFilt <- lmFit(esetCpFilterPFT, designIntrxFilt)
fitIntrxFilt <- eBayes(fitIntrxFilt)

resultsIntrxFilt <- decideTests(fitIntrxFilt, adjust.method="fdr", p.value=0.25)
summary(resultsIntrxFilt)

intrxInd <- which(resultsIntrxFilt[, 5] != 0)
intrx486 <- data.frame(fitIntrxFilt$genes[intrxInd,1])
colnames(intrx486) <- c("GeneSymbol")
intrx486$Group <- resultsIntrxFilt[intrInd, 5]
intrx486$Group[intrx486$Group == -1] <- 0

write.table(intrx486, file="intrxFilt486GeneSymbols", sep=",", row.names=FALSE, quote=FALSE)



## ----uniqueToBoth--------------------------------------------------------
# try method one; mash everything together into 2 groups

esetMash <- esetCpFilterPFT
esetMash$biGroup <- numeric(345)
esetMash$biGroup[esetMash$indicator==4] <- 1
esetMash$biGroup <- factor(esetMash$biGroup)

designMash <- model.matrix(~ 1 + AGEcalc + biGroup, data=esetMash)
colnames(designMash)[2:3] <- c("Age", "COPDCancer")

fitMash <- lmFit(esetMash, designMash)
fitMash <- eBayes(fitMash)

resultsMash <- decideTests(fitMash, adjust.method="fdr", p.value=0.05)
summary(resultsMash)

qs <- p.adjust(p=fitMash$p.value[, 3], "fdr")

pdf("biGroupBoxplots.pdf")
for(i in 1:length(which(resultsMash[, 3] != 0))){ 
  # generate the gene name
  ind <- which(resultsMash[, 3] != 0)[i]
 
  gn <- paste(fitMash$genes[ind, 1], ":", fitMash$genes[ind, 2], "\nq=", qs[ind], sep="")
 
  boxplot(exprs(esetMash)[ind, ] ~ esetMash$indicator,
          names=c("Disease Free", "Cancer", "COPD", "Cancer+COPD"),
          col=c("green", "red", "blue", "yellow"),
          xlab="Group", ylab="Expression",
          main=gn)
  }

dev.off()



## ----testingCEACAM5------------------------------------------------------

ceacam5ind <- match("CEACAM5", fitIntrx$genes[, 1])
ceacam5 <- exprs(esetSmDis)[790,]
esetSmDis$indicator
t.test(ceacam5[esetSmDis$indicator==1],ceacam5[esetSmDis$indicator==2])


