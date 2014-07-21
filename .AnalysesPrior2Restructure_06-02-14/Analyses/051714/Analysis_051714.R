# GENDER*CANCER interaction analysis in cancer only (no COPD) Allegro Cohort (IDE & CLIA)
# Jake Kantrowitz
# kantro@bu.edu, jacob.kantrowitz@gmail.com
# 04-02-14

# Goal(s) for the day:
# Run gender*cancer analysis models
######################################################################

init.dir <- getwd()
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Functions_v01.r")
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/data_phenotype_updated.RData")
cat("Loaded OK\n")

#cat("\nLoading cluster data for 657 patients; removed 43 w/o cancer or PY values\n")
#load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Analyses/030214/clusters.RData")
#load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Analyses/031114/clusPhen.RData")
#cat("Loaded OK\n")

save.packages()

## (1) find patients with or without cancer and without COPD
# optionally: split them by gender
cancer.ind <- which((data_phenotype$CANCER==1 | data_phenotype$CANCER==0) & data_phenotype$COPD==0)
can.gender0.ind <- which((data_phenotype$CANCER==1 | data_phenotype$CANCER==0) & data_phenotype$COPD==0 & data_phenotype$GENDER==0)
can.gender1.ind <- which((data_phenotype$CANCER==1 | data_phenotype$CANCER==0) & data_phenotype$COPD==0 & data_phenotype$GENDER==1)

cancer.UNF <- keep_patients(cancer.ind, data_expression, data_phenotype)
cancer.CV <- keep_patients(cancer.ind, data_filtered$cv, data_phenotype)
cancer.CVMN <- keep_patients(cancer.ind, data_filtered$cv_mn, data_phenotype)

can.gender0.UNF <- keep_patients(can.gender0.ind, data_expression, data_phenotype)
can.gender0.CV <- keep_patients(can.gender0.ind, data_filtered$cv, data_phenotype)
can.gender0.CVMN <- keep_patients(can.gender0.ind, data_filtered$cv_mn, data_phenotype)

can.gender1.UNF <- keep_patients(can.gender1.ind, data_expression, data_phenotype)
can.gender1.CV <- keep_patients(can.gender1.ind, data_filtered$cv, data_phenotype)
can.gender1.CVMN <- keep_patients(can.gender1.ind, data_filtered$cv_mn, data_phenotype)



covars.cancer <- list(list(list("BATCH", "RIN", "SMK", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "AGE", "GENDER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "AGE", "CANCER", "GENDER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "AGE", "CANCER", "GENDER", "CANCER*GENDER"), TRUE))
					
covars.can.gender <- list(list(list("BATCH", "RIN", "SMK", "AGE"), TRUE),
						list(list("BATCH", "RIN", "SMK", "AGE", "CANCER"), TRUE))
						

# Import and fix all of the new annotation information
setwd("/protected/projects/pulmarray/Allegro/COPD_Cancer/masterAnnotation/jk_csv")

clia <- read.csv("CLIA New 2014-05-02 1555 AZJ JK.csv")
ide <- read.csv("IDE New 2014-02-30 0008 AZJ (version 1) JK.csv")

# the names here don't exactly overlap oddly enough; two columns are not the same, otherwise exactly match
clia2 <- clia[,1:60]
ide2 <- ide[,1:60]

# merge the IDE and CLIA annotation information together
# this creates a data frame with 745 rows and 62 columns (i.e. all of the patients and all columns)
all.new <- merge(clia2,ide2, all.x=TRUE,all.y=TRUE)

# need to fix the barcodes so they are comparable and so the new can be ordered the same as the old
BARCODEfix <- function(f){
	a <- strsplit(as.character(f),"-")
	if (as.numeric(a[[1]][2])<10){
		a[[1]][2] <- paste("0", a[[1]][2], sep="")
	}
	newBARCODE <- paste(a[[1]][1], a[[1]][2], a[[1]][3], sep="-")
	return(newBARCODE)
}


bc <- data_phenotype$BARCODE
bc2 <- sapply(bc, FUN=BARCODEfix)

bc.ind <- match(bc2, all.new$BARCODE)
all.new.700 <- all.new[bc.ind,]

setwd(init.dir)
saveRDS(all.new.700, file="allNew700.RData")

#create a phenotype/annotation table to use for analysis
data_phenotype$BARCODEfixed <- bc2

# need to merge the all.new.700 and data_phenotype tables into new table
# merge by the fixed BARCODES (BARCODEfixed) and the original BARCODE from data_phenotype and all.new.700, respectively
# do not resort the new table (i.e. maintain the original order for comparability with order of expression data)
# creates a new table 
mergeAnnotation <- merge(data_phenotype, all.new.700, by.x="BARCODEfixed", by.y="BARCODE", sort=FALSE)

new.covars <- list(list(list("BATCH", "RIN", "AGEcalc", "GENDERc", "RACE4PFTs", "PYc", "SMKc", "FinalCaDXc", "COPD_EM_Bronchitisc", "FinalCaDXc*COPD_EM_Bronchitisc"), TRUE))

new.covars.indiv <- list(list(list("BATCH"), TRUE),
						list(list("RIN"), TRUE),
						list(list("AGEcalc"), TRUE),
						list(list("GENDERc"), TRUE),
						#list(list("RACE4PFTs"), TRUE), # has NAs
						#list(list("PYc"), TRUE), # has NAs and is a factor
						list(list("SMKc"), TRUE),
						list(list("FinalCaDXc"), TRUE),
						list(list("COPD_EM_Bronchitisc"), TRUE),
						list(list("FinalCaDXc", "COPD_EM_Bronchitisc","FinalCaDXc*COPD_EM_Bronchitisc"),TRUE))

						
remove.ind <- which(mergeAnnotation$FinalCaDXc=="DK" | mergeAnnotation$COPD_EM_Bronchitisc=="DK")
# results in removing 44 patients -> 656 remain
noNA.all <- remove_patients(remove.ind, data_expression, mergeAnnotation)

# remove DKs from Cancer
tmp <- noNA.all[[2]]$FinalCaDXc
tmp <- factor(tmp, levels=c(0,1))
noNA.all[[2]]$FinalCaDXc <- tmp

# remove DKs from COPD
tmp <- noNA.all[[2]]$COPD_EM_Bronchitisc
tmp <- factor(tmp, levels=c(0,1))
noNA.all[[2]]$COPD_EM_Bronchitisc <- tmp

# remove DKs from Gender
tmp <- noNA.all[[2]]$GENDERc
tmp <- factor(tmp, levels=c(0,1))
noNA.all[[2]]$GENDERc <- tmp

# make SMKc as a factor
tmp <- as.factor(noNA.all[[2]]$SMKc)
noNA.all[[2]]$SMKc <- tmp

#results.indiv <- analyze_lists(new.covars.indiv, noNA.all[[1]], noNA.all[[2]])	
#results.indiv.total <- analyze_lists(new.covars.indiv, noNA.all[[1]], noNA.all[[2]])	

# batch -> no results
# RIN -> 17945 at FDR 0.25
# AGEcalc -> 2419 at FDR 0.25
# GENDERc -> 6380 at FDR 0.25
# SMKc -> 12276 at FDR 0.25
# FinalCaDXc1 -> 744 at FDR 0.25
# COPD_EM_Bronchitisc1 -> 7230  at FDR 0.25
# INTERACTION EFFECT MODEL
#	FinalCaDXc1 -> 1001 at FDR 0.25
#	COPD_EM_Bronchitisc1 -> 4595 at FDR 0.25
#	FinalCaDXc1:COPD_EM_Bronchitisc1 -> 0 at FDR 0.25

#	FinalCaDXc1 -> 2471 at p < 0.05
#	COPD_EM_Bronchitisc1 -> 4314 at p < 0.05
#	FinalCaDXc1:COPD_EM_Bronchitisc1 -> 1769 at p < 0.05 (985 expected by chance)

# SMK status appears to play a large role in the genes significant for COPD*CANCER interaction
# try including smoking in the model
# try mapping residuals from smoking only model

# pull smoking residuals from models run above
#resid.smk <- results.indiv.total[[5]]$model_resids

#save.genes(results.indiv.total[[ind]]$cutoffs_p[[2]]$'FinalCaDXc1:COPD_EM_Bronchitisc1', noNA.all[[1]], noNA.all[[2]], filename="copdCancerIntrx.p.01")
#save.genes(results.indiv.total[[ind]]$cutoffs_p[[3]]$'FinalCaDXc1:COPD_EM_Bronchitisc1', noNA.all[[1]], noNA.all[[2]], filename="copdCancerIntrx.p.005")

#save.genes(results.indiv.total[[ind]]$cutoffs_p[[2]]$'FinalCaDXc1:COPD_EM_Bronchitisc1', resid.smk, noNA.all[[2]], filename="copdCancerIntrx.p.01.19684.403genes.resid.smk")
#save.genes(results.indiv.total[[ind]]$cutoffs_p[[3]]$'FinalCaDXc1:COPD_EM_Bronchitisc1', resid.smk, noNA.all[[2]], filename="copdCancerIntrx.p.005.19684.202genes.resid.smk")


#new.covars.models <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc", "COPD_EM_Bronchitisc", "FinalCaDXc*COPD_EM_Bronchitisc"), TRUE))
#results.models1 <- analyze_lists(new.covars.models, noNA.all[[1]], noNA.all[[2]])



# Cancer and SMK are confounded; fisher exact gives p-value 0.008
# COPD and SMK are confounded; fisher exact gives p-value 0.01
# going to add an interaction term in the model
#new.covars.models2 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc", "COPD_EM_Bronchitisc", "FinalCaDXc*COPD_EM_Bronchitisc", "FinalCaDXc*SMKc", "COPD_EM_Bronchitisc*SMKc"), TRUE))
#results.models2 <- analyze_lists(new.covars.models2, noNA.all[[1]], noNA.all[[2]])




# try using a different definition of COPD
# try using Ratio < 0.7 definition
remove.ind2 <- which(is.na(mergeAnnotation$COPD2_R7) | mergeAnnotation$COPD2_R7=="DK" | mergeAnnotation$FinalCaDXc=="DK")
noNA.copdR7 <- remove_patients(remove.ind2, data_expression, mergeAnnotation)

# remove DKs from Cancer
tmp <- noNA.copdR7[[2]]$FinalCaDXc
tmp <- factor(tmp, levels=c(0,1))
noNA.copdR7[[2]]$FinalCaDXc <- tmp

# remove DKs from COPD
tmp <- noNA.copdR7[[2]]$COPD2_R7
tmp <- factor(tmp, levels=c(0,1))
noNA.copdR7[[2]]$COPD2_R7 <- tmp

# remove DKs from Gender
tmp <- noNA.copdR7[[2]]$GENDERc
tmp <- factor(tmp, levels=c(0,1))
noNA.copdR7[[2]]$GENDERc <- tmp

# make SMKc as a factor
tmp <- as.factor(noNA.copdR7[[2]]$SMKc)
noNA.copdR7[[2]]$SMKc <- tmp

#new.covars.models3 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc", "COPD2_R7", "FinalCaDXc*COPD2_R7", "FinalCaDXc*SMKc", "COPD2_R7*SMKc"), TRUE))

#results.models.COPDR7 <- analyze_lists(new.covars.models3, noNA.copdR7[[1]], noNA.copdR7[[2]])

# need to use smoking residuals from heatmap generation
#new.covars.smk <- list(list(list("SMKc"), TRUE))
#results.smk.COPDR7 <- analyze_lists(new.covars.smk, noNA.copdR7[[1]], noNA.copdR7[[2]])

#new.covars.resid <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc"), TRUE))
#results.resid.COPDR7 <- analyze_lists(new.covars.resid, noNA.copdR7[[1]], noNA.copdR7[[2]])
#resid.full <- results.resid.COPDR7[[1]]$model_resids

# pull out the smoking resids for heatmap
#resid.smk <- results.smk.COPDR7[[1]]$model_resids

#c2 <- return_cluster(results.models.COPDR7[[1]]$cutoffs_p[[4]]$'FinalCaDXc1:COPD2_R71', resid.smk, noNA.copdR7[[2]])

# find vaguely around 400 genes
#ccIntFDR.15 <- which(results.models.COPDR7[[1]]$model_fdr[,7] < 0.15) # this is 338 genes

#c2.2 <- return_cluster(ccIntFDR.15, resid.smk, noNA.copdR7[[2]])

ind0 <- which(noNA.copdR7[[2]]$COPD2_R7==0 & noNA.copdR7[[2]]$FinalCaDXc==0)
ind1 <- which(noNA.copdR7[[2]]$COPD2_R7==1 & noNA.copdR7[[2]]$FinalCaDXc==0)
ind2 <- which(noNA.copdR7[[2]]$COPD2_R7==0 & noNA.copdR7[[2]]$FinalCaDXc==1)
ind3 <- which(noNA.copdR7[[2]]$COPD2_R7==1 & noNA.copdR7[[2]]$FinalCaDXc==1)


# need to fix the code for generating the heatmaps; currently uses the indicator function, which needs to be replaced with new info
# will replace cClus variable in save.genes function with new indicator function
# generate new indicator function
indicator.new <- numeric(nrow(noNA.copdR7[[2]]))
indicator.new[ind0] <- 0
indicator.new[ind1] <- 1
indicator.new[ind2] <- 2
indicator.new[ind3] <- 3
noNA.copdR7[[2]]$indicator <- as.factor(indicator.new)

allOrNone <- numeric(nrow(noNA.copdR7[[2]]))
allOrNone[ind0] <- 0
allOrNone[ind1] <- 0
allOrNone[ind2] <- 0
allOrNone[ind3] <- 1
noNA.copdR7[[2]]$allOrNone <- as.factor(allOrNone)

#new.covars.allOrNone <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "allOrNone"), TRUE))
#results.model.allOrNone <- analyze_lists(new.covars.allOrNone, noNA.copdR7[[1]], noNA.copdR7[[2]])

#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1, noNA.copdR7[[1]], noNA.copdR7[[2]], filename="allOrNone.fdr.25.19684.429genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[2]]$allOrNone1, noNA.copdR7[[1]], noNA.copdR7[[2]], filename="allOrNone.fdr.20.19684.429genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[3]]$allOrNone1, noNA.copdR7[[1]], noNA.copdR7[[2]], filename="allOrNone.fdr.10.19684.429genes")

# need to use smoking residuals from heatmap generation
#new.covars.smk <- list(list(list("SMKc"), TRUE))
#results.smk.COPDR7 <- analyze_lists(new.covars.smk, noNA.copdR7[[1]], noNA.copdR7[[2]])

#smk.resid <- results.smk.COPDR7[[1]]$model_resids
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1, smk.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.25.smk.resid.19684.429genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[2]]$allOrNone1, smk.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.20.smk.resid.19684.429genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[3]]$allOrNone1, smk.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.10.smk.resid.19684.429genes")

# want to use all-term model residuals
#new.covars.resid2 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc"), TRUE))
#results.all.resid <- analyze_lists(new.covars.resid2, noNA.copdR7[[1]], noNA.copdR7[[2]])

#all.resid <- results.model.allOrNone[[1]]$model_resids
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1, all.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.25.all.resid.19684.429genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[2]]$allOrNone1, all.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.20.all.resid.19684.280genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[3]]$allOrNone1, all.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.10.all.resid.19684.84genes")

#cc2 <- return_cluster(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1, all.resid, noNA.copdR7[[2]])
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1[cc2==1], all.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.25.all.resid.19684.429genes.cluster1.243genes")
#save.genes(results.model.allOrNone[[1]]$cutoffs_fdr[[1]]$allOrNone1[cc2==2], all.resid, noNA.copdR7[[2]], filename="allOrNone.fdr.25.all.resid.19684.429genes.cluster2.186genes")

noNA.copdR7.CS <- remove_patients(which(noNA.copdR7[[2]]$AllegroCorticosteroidsc=="DK"), noNA.copdR7[[1]], noNA.copdR7[[2]])
noNA.copdR7.CS[[2]]$AllegroCorticosteroidsc <- factor(noNA.copdR7.CS[[2]]$AllegroCorticosteroidsc, levels=c(1,0))

#new.covars.allOrNone2 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "allOrNone", "AllegroCorticosteroidsc"), TRUE))
#results.model.allOrNone2 <- analyze_lists(new.covars.allOrNone2, noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]])

#save.genes(results.model.allOrNone2[[1]]$cutoffs_fdr[[1]]$allOrNone1, noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]], filename="allOrNone2.fdr.25.modelCS.18genes")

# want to use all-term model residuals + CS
#new.covars.resid3 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "AllegroCorticosteroidsc"), TRUE))
#results.all.resid3 <- analyze_lists(new.covars.resid3, noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]])
#all.resid3 <- results.all.resid3[[1]]$model_resids

#save.genes(results.model.allOrNone2[[1]]$cutoffs_fdr[[1]]$allOrNone1, all.resid3, noNA.copdR7.CS[[2]], filename="allOrNone2.fdr.25.allCS.resid.modelCS.18genes")
#save.genes(results.model.allOrNone2[[1]]$cutoffs_p[[2]]$allOrNone1, all.resid3, noNA.copdR7.CS[[2]], filename="allOrNone2.p.01.allCS.resid.modelCS.333genes")

new.covars.models4 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc", "COPD2_R7", "FinalCaDXc*COPD2_R7", "FinalCaDXc*SMKc", "COPD2_R7*SMKc", "AllegroCorticosteroidsc"), TRUE))
results.all.models4 <- analyze_lists(new.covars.models4, noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]])

fdr16ind <- which(results.all.models4[[1]]$model_fdr[,8] < .16)
save.genes(fdr16ind, noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]], filename="CC.fdr.16.101genes.modelCS")

ind0 <- which(noNA.copdR7.CS[[2]]$indicator==0)
ind1 <- which(noNA.copdR7.CS[[2]]$indicator==2)
no02 <- remove_patients(union(ind0,ind1), noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]])

no02[[2]]$indicator <- factor(no02[[2]]$indicator, levels=c(1,3))
new.covars.models5 <- list(list(list("RIN", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc", "FinalCaDXc*SMKc", "AllegroCorticosteroidsc"), TRUE))
results.all.models5 <- analyze_lists(new.covars.models5, no02[[1]], no02[[2]])


#save.genes(results.all.models4[[1]]$cutoffs_fdr[[2]]$'FinalCaDXc1:COPD2_R71', noNA.copdR7.CS[[1]], noNA.copdR7.CS[[2]], filename="CC.fdr.2.362genes.modelCS")
#save.genes(results.all.models4[[1]]$cutoffs_fdr[[2]]$'FinalCaDXc1:COPD2_R71', all.resid3, noNA.copdR7.CS[[2]], filename="CC.fdr.2.362genes.modelCS")


results.all.models4[[1]]$cutoffs_fdr[[1]]$'FinalCaDXc1:COPD2_R71'
ttt <- results.all.models4[[1]]$cutoffs_fdr[[1]]$'FinalCaDXc1:COPD2_R71'
results.all.models4[[1]]$model_fdr[ttt,8]
ttt2 <- which(results.all.models4[[1]]$model_fdr[ttt,8] < .157191)
ttt[ttt2[1]]
boxplot(noNA.copdR7.CS[[1]][ttt[ttt2[tind]],]~noNA.copdR7.CS[[2]]$indicator)
# boxplot website: http://www.statmethods.net/graphs/boxplot.html
pdf("boxplot.significant.for.copdCancer.pdf")
boxplot(noNA.copdR7.CS[[1]][ttt[ttt2[50]],]~noNA.copdR7.CS[[2]]$indicator, names=c("Healthy", "COPD", "Cancer", "COPD+Cancer"), notch=TRUE,
	col=c("green", "yellow", "orange", "red"), xlab="Group", ylab="Gene Expression", main="Gene significant for COPD*Cancer interaction")
dev.off()

# find 548 of the significant ones
#results.all.models4[[1]]$cutoffs_fdr[[1]]$'FinalCaDXc1:COPD2_R71'[548]
# ind 13266
noNA.copdR7.CS[[1]]
noNA.copdR7.CS[[1]][13266,]


# some helpful examples of running covariate list models, saving genes (txt and heatmap) and returning clusters
#results.UNF <- analyze_lists(covars, data.noNA[[1]], data.noNA[[2]])
#save.genes(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], "cvmnFILT.noNA.smk+copd.p.01")
#cs <- return_cluster(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], colClus="unsupervised", type=ROWS)