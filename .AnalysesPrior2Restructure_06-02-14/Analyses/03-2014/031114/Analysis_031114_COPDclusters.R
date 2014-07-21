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
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/data_phenotype_updated.RData")
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
clia <- read.csv("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/masterAnnotation/CLIA_clinical_data_master_copy.csv")
ide <- read.csv("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/masterAnnotation/IDE_clinical_data_master_copy.csv")

############## Remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA.UNF <- remove_patients(ind.NA, data_expression, data_phenotype)
noNA.CV <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)
noNA.CVMN <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)


# set covariates to use for analysis
# analyze only cancer patients - only non-COPD patients
covars <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), FALSE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD"), FALSE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER"), FALSE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD"), FALSE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE)
				)

############## Run with unfiltered data
# remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA.UNF <- remove_patients(ind.NA, data_expression, data_phenotype)

#results.UNF <- analyze_lists(covars, noNA.UNF[[1]], noNA.UNF[[2]])


############## Run with CV filtered data
# remove NAs in cancer or PY status
noNA.CV <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)

#results.CV <- analyze_lists(covars, noNA.CV[[1]], noNA.CV[[2]])

############## Run with CVMN filtered data
# remove NAs in cancer or PY status
noNA.CVMN <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)

results.CVMN <- analyze_lists(covars, noNA.CVMN[[1]], noNA.CVMN[[2]])





######################################################################################################
######################################################################################################

# run linear model with patient clusters as only term or as co-variate 
# cluster based on interaction model CVMN filtered data
clusters <- return_cluster(results.CVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', results.CVMN[[1]]$model_resids, noNA[[2]], colClus="unsupervised", type=COLUMNS)
clusters <- as.factor(clusters)

clusters2 <- return_cluster(results.CVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', noNA[[1]], noNA[[2]], colClus="unsupervised", type=COLUMNS)

clusPhen <- noNA[[2]]
clusPhen$clusters <- clusters

covar.clusters <- list(list(list("clusters"), TRUE),
					list(list("clusters", "BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE))
					
results.clusters <- analyze_lists(covar.clusters, noNA.CVMN[[1]], clusPhen)


# see if there is a relationship between gene mean expression and cluster coefficient
clusters.ind <- results.clusters[[2]]$cutoffs_fdr[[1]]$clusters2
clusters.p <- results.clusters[[2]]$model_ps[clusters.ind,1]
clusters.sort <- sort(clusters.p, index.return=TRUE)
clusters.ix.2000 <- clusters.sort$ix[1:2000]
clusters.ind.p.2000 <- clusters.ind[clusters.ix.2000]

cluster.coefficients <- numeric(2000)
for(i in 1:2000)
{
	cluster.coefficients[i] <- results.clusters[[2]]$models[[clusters.ind.p.2000[i]]]$coefficients[2]
}

cluster.means <- numeric(2000)
for(i in 1:2000)
{
	cluster.means[i] <- mean(noNA.CVMN[[1]][clusters.ind.p.2000[i],])
}

save(cluster.means, file="cluster.means.RData")
save(cluster.coefficients, file="cluster.coefficients.RData")
save(clusters.ind.p.2000, file="clusters.ind.p.2000.RData")
save(clusters.ix.2000, file="clusters.ix.2000.RData")
save(clusters.sort, file="clusters.sort.RData")
save(clusters.p, file="clusters.p.RData")
save(clusters.ind, file="clusters.ind.RData")
save(clusters, file="clusters.RData")
save(clusPhen, file="clusPhen.RData")

pdf("geneMeanVcoeff.pdf")
plot(cluster.means, cluster.coefficients)
dev.off()

# see if there is a relationship between gene mean expression and RIN coefficient
rin.ind <- results.clusters[[2]]$cutoffs_fdr[[1]]$RIN
rin.p <- results.clusters[[2]]$model_ps[rin.ind,1]
rin.sort <- sort(rin.p, index.return=TRUE)
rin.ix.2000 <- rin.sort$ix[1:2000]
rin.ind.p.2000 <- rin.ind[rin.ix.2000]

rin.coefficients <- numeric(2000)
for(i in 1:2000)
{
	rin.coefficients[i] <- results.clusters[[2]]$models[[rin.ind.p.2000[i]]]$coefficients[7]
}

rin.means <- numeric(2000)
for(i in 1:2000)
{
	rin.means[i] <- mean(noNA.CVMN[[1]][rin.ind.p.2000[i],])
}

pdf("geneMeanVcoeffRIN.pdf")
plot(rin.means, rin.coefficients)
dev.off()

# see if there is a relationship between gene mean expression and SMK coefficient
smk.ind <- results.clusters[[2]]$cutoffs_fdr[[1]]$SMK2
smk.p <- results.clusters[[2]]$model_ps[smk.ind,1]
smk.sort <- sort(smk.p, index.return=TRUE)
smk.ix.2000 <- smk.sort$ix[1:2000]
smk.ind.p.2000 <- smk.ind[smk.ix.2000]

smk.coefficients <- numeric(2000)
for(i in 1:2000)
{
	smk.coefficients[i] <- results.clusters[[2]]$models[[smk.ind.p.2000[i]]]$coefficients[7]
}

smk.means <- numeric(2000)
for(i in 1:2000)
{
	smk.means[i] <- mean(noNA.CVMN[[1]][smk.ind.p.2000[i],])
}

pdf("geneMeanVcoeffSMK.pdf")
plot(smk.means, smk.coefficients)
dev.off()

# see if there is a relationship between gene mean expression and GENDER coefficient
gender.ind <- results.clusters[[2]]$cutoffs_fdr[[1]]$GENDER1
gender.p <- results.clusters[[2]]$model_ps[gender.ind,1]
gender.sort <- sort(gender.p, index.return=TRUE)
gender.ix.2000 <- gender.sort$ix[1:2000]
gender.ind.p.2000 <- gender.ind[gender.ix.2000]

gender.coefficients <- numeric(2000)
for(i in 1:2000)
{
	gender.coefficients[i] <- results.clusters[[2]]$models[[gender.ind.p.2000[i]]]$coefficients[10]
}

gender.means <- numeric(2000)
for(i in 1:2000)
{
	gender.means[i] <- mean(noNA.CVMN[[1]][gender.ind.p.2000[i],])
}

pdf("geneMeanVcoeffGENDER.pdf")
plot(gender.means, gender.coefficients)
dev.off()




smk.cor <- cor.test(smk.means, smk.coefficients)
cluster.cor <- cor.test(cluster.means, cluster.coefficients)
rin.cor <- cor.test(rin.means, rin.coefficients)
gender.cor <- cor.test(gender.means, gender.coefficients)




