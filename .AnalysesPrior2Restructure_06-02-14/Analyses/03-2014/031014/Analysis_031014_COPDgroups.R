# Allegro Analysis - examining gene sets from Grant Duclos
# Jake Kantrowitz
# 03-04-14

# GOAL look for genes differing by cluster 1/2
# Find pathways differing between the clusters
# Find clinical variables differing between the clusters

######################################################################
############## HELPER FUNCTIONS
fold.change <- function(temp, groups, g1=1, g2=2)
{
	g1.ind <- which(groups==g1)
	g2.ind <- which(groups==g2)
	num.genes <- nrow(temp)
	fcs <- matrix(nrow=num.genes, ncol=3)
	
	for(gene in 1:num.genes)
	{
		fcs[gene, 1] <- mean(temp[gene, g1.ind])
		fcs[gene, 2] <- mean(temp[gene, g2.ind])
		fcs[gene, 3] <- fcs[gene,1] - fcs[gene,2]
	}
	
	return(fcs)
}

save.genes <- function(inds, exprData, phen, filename, rClusters=NULL, mtd="average", cClus=phen$indicator,man="Figure")
{
	save_entrez(inds, rownames(exprData), paste(filename, ".txt", sep=""))
	pdf(paste(filename,".pdf", sep=""))
	generate_heatmap(inds, exprData, phen, rowClusters=rClusters, mthd=mtd, colClus=cClus, mn=man)
	dev.off()
	cat("Saved entrez IDs and saved heatmap for ", length(inds), " genes.\n", sep="")
}

######################################################################

init.dir <- getwd()
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Functions_v01.r")
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("../../data_phenotype_updated.RData")
cat("Loaded OK\n")

cat("\nLoading cluster data for 657 patients\n")
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Analyses/030214/clusters.RData")
cat("Loaded OK\n")

cat("\nSaving packages\n")
int <- installed.packages()
dt <- as.character(Sys.Date())
write.table(int, file=paste(dt,"_packages.csv", sep=""), sep=",")
cat("Saved OK\n")

# load clia and ide master annotation files
clia <- read.csv("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/CLIA_clinical_data_master_copy.csv")
ide <- read.csv("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/IDE_clinical_data_master_copy.csv")

############## Remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA.UNF <- remove_patients(ind.NA, data_expression, data_phenotype)
noNA.CV <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)
noNA.CVMN <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)


# set covariates to use for analysis
# analyze only cancer patients - only non-COPD patients
covars <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE)
				)

############## Run with unfiltered data
# remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA <- remove_patients(ind.NA, data_expression, data_phenotype)

resultsUNF <- analyze_lists(covars, noNA[[1]], noNA[[2]])


############## Run with CV filtered data
# remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)

resultsCV <- analyze_lists(covars, noNA[[1]], noNA[[2]])

############## Run with CVMN filtered data
# remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)

resultsCVMN <- analyze_lists(covars, noNA[[1]], noNA[[2]])





######################################################################################################
######################################################################################################

# run linear model with patient clusters as only term or as co-variate 
# cluster based on interaction model CVMN filtered data
clusters <- return_cluster(resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', resultsCVMN[[1]]$model_resids, noNA[[2]], colClus="unsupervised", type=COLUMNS)
clusters <- as.factor(clusters)

clusters2 <- return_cluster(resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', noNA[[1]], noNA[[2]], colClus="unsupervised", type=COLUMNS)

clusPhen <- noNA[[2]]
clusPhen$clusters <- clusters

covar.clusters <- list(list(list("clusters"), TRUE),
					list(list("clusters", "BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE))
					
results.clusters <- analyze_lists(covar.clusters, noNA[[1]], clusPhen)
##

fold.change(noNA[[1]][results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2,], clusPhen$clusters, 1,2)

fcs <- fold.change(noNA[[1]][results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2,], clusPhen$clusters, 1,2)

fcs1.ind <- results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2[which(abs(fcs[,3]) > 1)]
save.genes(fcs1.ind, noNA[[1]], clusPhen, "clusters.fdr.05.fc.1")

fcs0.7.ind <- results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2[which(abs(fcs[,3]) > .7)]
save.genes(fcs.7.ind, noNA[[1]], clusPhen, "clusters.fdr.05.fc.7")

fcs.5.up.ind <- results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2[which(fcs[,3] > .5)]
save.genes(fcs.5.up.ind, noNA[[1]], clusPhen, "clusters.fdr.05.fc.5.up")


fcs.5.down.ind <- results.clusters[[1]]$cutoffs_fdr[[4]]$clusters2[which(fcs[,3] < -.5)]
save.genes(fcs.5.down.ind, noNA[[1]], clusPhen, "clusters.fdr.05.fc.5.down")

# Fisher exact tests against the clusters
fe.cancer <- fisher.test(clusPhen$clusters,clusPhen$CANCER) # p=0.6551876
fe.copd <- fisher.test(clusPhen$clusters,clusPhen$COPD) # p=1.382e-05
fe.healthy <- fisher.test(clusPhen$clusters, clusPhen$HEALTHY) # p=0.0004813916
fe.smk <- fisher.test(clusPhen$clusters, clusPhen$SMK) # p=0.1205316
fe.gender <- fisher.test(clusPhen$clusters, clusPhen$GENDER) # p=0.7685393

# T-tests against the clusters
tt.age <- t.test(clusPhen$AGE[clusPhen$clusters==1], clusPhen$AGE[clusPhen$clusters==2]) # p=0.2640838
tt.py <- t.test(clusPhen$PY[clusPhen$clusters==1], clusPhen$PY[clusPhen$clusters==2]) # p=0.8647443



# Generate some heatmaps with clusters supervising
clabels.copd <- cbind(
		"COPD" = copd_colors[clusPhen$COPD],
		"Cluster" = cancer_colors[clusPhen$clusters]
	)
	
clabels.cancer <- cbind(
		"COPD" = copd_colors[clusPhen$CANCER],
		"Cluster" = cancer_colors[clusPhen$clusters]
	)

pdf("clusters.fdr.05.fc.5.up.COPD.pdf")
heatmap3(noNA[[1]][fcs.5.up.ind,], col=bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering=clusPhen$clusters, ColSideColors=clabels.copd)
dev.off()

pdf("clusters.fdr.05.fc.5.down.COPD.pdf")
heatmap3(noNA[[1]][fcs.5.down.ind,], col=bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering=clusPhen$clusters, ColSideColors=clabels.copd)
dev.off()

pdf("clusters.fdr.05.fc.5.up.CANCER.pdf")
heatmap3(noNA[[1]][fcs.5.up.ind,], col=bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering=clusPhen$clusters, ColSideColors=clabels.cancer)
dev.off()

pdf("clusters.fdr.05.fc.5.down.CANCER.pdf")
heatmap3(noNA[[1]][fcs.5.down.ind,], col=bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering=clusPhen$clusters, ColSideColors=clabels.cancer)
dev.off()


######### WANT TO MAKE SURE THE CLUSTERS ARE THE MOST SENSICAL
# non residuals
heatmap3(noNA[[1]][resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1',], col=bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering="unsupervised", ColSideColors=clabels.cancer)
clusters <- return_cluster(resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', resultsCVMN[[1]]$model_resids, noNA[[2]], colClus="unsupervised", type=COLUMNS)

clusters2 <- return_cluster(resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', noNA[[1]], noNA[[2]], colClus="unsupervised", type=COLUMNS)
clusters2 <- as.factor(clusters2)

clabels.clusters <- clusters2
# residuals
heatmap3(resultsCVMN[[1]]$model_resids[resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1',], col=bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method="average"), col.clustering="unsupervised", ColSideColors=clabels.copd, main="FIGURE")

# use clusters2
heatmap3(resultsCVMN[[1]]$model_resids[resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1',], col=bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method="average"), col.clustering="unsupervised", ColSideColors=clabels.clusters, main="FIGURE")



