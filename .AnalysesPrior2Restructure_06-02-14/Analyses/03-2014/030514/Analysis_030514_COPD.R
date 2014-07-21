# Allegro Analysis - examining gene sets from Grant Duclos
# Jake Kantrowitz
# 03-04-14

# From Grant 02-13-14
# "cil.entrez.txt" consists of genes (621 total) that appear to be involved in ciliogenesis,
# the onset of their expression is delayed during differentiation of airway cells derived
# from smokers. My thought was that since there are indications that some
# ciliogenesis-related genes (mir4423 for example) are down-regulated in airway tissue
# taken from cancer patients, maybe you could pick up a large-scale down-regulation of
# the process in your cancer data

# GOAL: look for genes down-regulated in cancer patients compared to healthy patients

# "er.entrez.txt" consists of genes (521 total) that appear to be involved in endoplasmic
# reticulum stress, the unfolded protein response, and maybe overall expansion
# of the ER...these genes are very high in smoker cells prior to commitment to
# differentiation. I've come across a number of papers that emphasize the potential
# importance of ER stress in development of COPD so maybe these genes would be
# particularly high in samples taken from COPD patients.

# GOAL: look for genes up-regulated in COPD patients compared to healthy patients

# As far as lay-out goes, I think it might be helpful to separate cancer from copd when
# making the heatmaps. Maybe arrange it so that it goes from former smokers to current
# smokers to cancer or copd left to right?

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
setwd("../../")
source("Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("data_phenotype_updated.RData")
cat("Loaded OK\n")
setwd(init.dir)
int <- installed.packages()
dt <- as.character(Sys.Date())
write.table(int, file=paste(dt,"_packages.csv", sep=""), sep=",")


############## Remove NAs in cancer or PY status
ind.NA <- which(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
noNA.UNF <- remove_patients(ind.NA, data_expression, data_phenotype)
noNA.CV <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)
noNA.CVMN <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)


############## Load cilia and endoplasmic reticulum related gene sets from Grant Duclos
cilia.entrez <- as.character(as.matrix(read.csv("../../GeneSets/cil.entrez.txt", head=FALSE)))
er.entrez <- as.character(as.matrix(read.csv("../../GeneSets/er.entrez.txt", head=FALSE)))

array.entrez <- return_entrez(rownms=rownames(data_expression))

cilia.ind <- match(cilia.entrez, array.entrez)[!is.na(match(cilia.entrez, array.entrez))]
er.ind <- match(er.entrez, array.entrez)[!is.na(match(er.entrez, array.entrez))]

# Pull cilia/er genes from unfiltered data that has had NAs removed
cilia.genes <- keep_genes(cilia.ind, noNA.UNF[[1]])
er.genes <- keep_genes(er.ind, noNA.UNF[[1]])




############## Set covariates to use for analysis
covars <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE)
				)

covars.cancer <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE)
				)

covars.copd <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE)
				)


############## Setup patients subsets to run cancer/healthy and copd/healthy
ind.copd <- which(noNA.UNF[[2]]$COPD==1)
ind.cancer <- which(noNA.UNF[[2]]$CANCER==1)

set.cancerControl <- remove_patients(ind.copd, cilia.genes, noNA.UNF[[2]])
set.copdControl <- remove_patients(ind.cancer, er.genes, noNA.UNF[[2]])

# Run cilia genes in cancer/healthy patients
# GOAL: look for genes down-regulated in cancer patients compared to healthy patients
# find genes with fold change > 0 (i.e. down in group 2, in this case cancer)
fc.cancer <- fold.change(set.cancerControl[[1]], set.cancerControl[[2]]$CANCER, g1=0, g2=1)
down.cilia.cancer <- which(fc.cancer[,3] > 0)
set.cancerControl.down <- keep_genes(down.cilia.cancer, set.cancerControl[[1]])

# running analysis with 351 genes in 417 patients (cancer n=303)
# @ p < 0.05 expect 18 by chance
results.cilia.cancerControl.down <- analyze_lists(covars.cancer, set.cancerControl.down, set.cancerControl[[2]])

############## Run er genes in copd/healthy patients
# GOAL: look for genes up-regulated in COPD patients compared to healthy patients
# find genes with fold change < 0 (i.e. up in group 2, in this case copd)
fc.copd <- fold.change(set.copdControl[[1]], set.copdControl[[2]]$COPD, g1=0, g2=1)
up.er.copd <- which(fc.copd[,3] < 0)
set.copdControl.up <- keep_genes(up.er.copd, set.copdControl[[1]])

# running analysis with 313 genes in 163 patients (copd n=49)
# @ p <0.05 expect 16 by chance
results.er.copdControl.up <- analyze_lists(covars.copd, set.copdControl.up, set.copdControl[[2]])


# print some heatmaps for cilia/cancer and save the gene lists
pdf("cilia.cancerControl.up.p.05.pdf")
generate_heatmap(results.cilia.cancerControl.down[[2]]$cutoffs_p[[1]]$CANCER1, set.cancerControl.down, set.cancerControl[[2]])     
dev.off()

pdf("cilia.cancerControl.up.p.05.resid.pdf")
generate_heatmap(results.cilia.cancerControl.down[[2]]$cutoffs_p[[1]]$CANCER1, results.cilia.cancerControl.down[[1]]$model_resids, set.cancerControl[[2]])
dev.off()

save_entrez(results.cilia.cancerControl.down[[2]]$cutoffs_p[[1]]$CANCER1, rownames(set.cancerControl.down), "cilia.cancerControl.up.p.05.txt")

# print some heatmaps for er/copd and save the gene lists
pdf("er.copdControl.up.fdr.25.pdf")
generate_heatmap(results.er.copdControl.up[[2]]$cutoffs_fdr[[1]]$COPD1, set.copdControl.up, set.copdControl[[2]])
dev.off()

pdf("er.copdControl.up.fdr.25.resid.pdf")
generate_heatmap(results.er.copdControl.up[[2]]$cutoffs_fdr[[1]]$COPD1, results.er.copdControl.up[[1]]$model_resids, set.copdControl[[2]])
dev.off()

save_entrez(results.er.copdControl.up[[2]]$cutoffs_fdr[[1]]$COPD1, rownames(set.copdControl.up), "er.copdControl.up.fdr.25.txt")



# combine sets and re-analyze
union.erCilia <- union(cilia.ind, er.ind)
ciliaEr.genes <- keep_genes(union.erCilia, noNA.UNF[[1]])
results.ciliaEr <- analyze_lists(covars, ciliaEr.genes, noNA.UNF[[2]])


############## Run with unfiltered data
#results.cilia <- analyze_lists(covars, cilia.genes, noNA.UNF[[2]])
#results.er <- analyze_lists(covars, er.genes, noNA.UNF[[2]])

# 606 cilia genes = 31 by chance at p < 0.05

# 511 endo retic genes = 26 by chance at p < 0.05

, 













######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# calculate fold change
# already working with log data
# thus, fold change calculated by subtraction (dividing non-log transformed values)
# i.e. mean(Group 1) - mean(Group 2)
# generally take the absolute value
# fold change = | mean(Group 1) - mean(Group 2)|

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

save.models <- function(object)
{
	fma <- gsub(" ", "", object$formula)
	num.genes <- as.characater(length(object$models))
	filename <- paste(num.genes, fma, ".RData", sep="")
	save(object$models, file)
}

# run linear model with patient clusters as only term or as co-variate 
# cluster based on interaction model CVMN filtered data
clusters <- return_cluster(resultsCVMN[[5]]$cutoffs_fdr[[1]]$'CANCER1:COPD1', resultsCVMN[[1]]$model_resids, noNA[[2]], colClus="unsupervised", type=COLUMNS)
clusters <- as.factor(clusters)
clusPhen <- noNA[[2]]
clusPhen$clusters <- clusters

covar.clusters <- list(list(list("clusters"), TRUE),
					list(list("clusters", "BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE))
					
results.clusters <- analyze_lists(covar.clusters, noNA[[1]], clusPhen)

fcs <- fold.change(noNA[[1]], clusPhen$clusters)

# 8758 genes significant for clusters at fdr < 0.05
clusters.fdr.05.ind <- results.clusters[[2]]$cutoffs_fdr[[4]]$clusters2
clusters.fdr.25.ind <- results.clusters[[2]]$cutoffs_fdr[[1]]$clusters2

clusters.fdr.05.fcs <- fcs[clusters.fdr.05.ind, 3]
clusters.fdr.25.fcs <- fcs[clusters.fdr.25.ind, 3]

clusters.05.abs.fcs <- abs(clusters.fdr.05.fcs)
clusters.25.abs.fcs <- abs(clusters.fdr.25.fcs)

clusters.fdr.05.fc.1.ind <- clusters.fdr.05.ind[which(clusters.05.abs.fcs >= .7)]
clusters.fdr.25.fc.1.ind <- clusters.fdr.25.ind[which(clusters.25.abs.fcs >= .7)]

# permute clusters to see if effect is real or random
# see if any clinical variables associate with clusters

clusters.perm2 <- sample(clusters, 657, replace=FALSE)
clusters.perm3 <- sample(clusters, 657, replace=FALSE)
clusters.perm4 <- sample(clusters, 657, replace=FALSE)

clusPhen2 <- noNA[[2]]
clusPhen2$clusters <- clusters.perm2

clusPhen3 <- noNA[[2]]
clusPhen3$clusters <- clusters.perm3

clusPhen4 <- noNA[[2]]
clusPhen4$clusters <- clusters.perm4

results.clusters.perm2 <- analyze_lists(covar.clusters, noNA[[1]], clusPhen2)
results.clusters.perm3 <- analyze_lists(covar.clusters, noNA[[1]], clusPhen3)
results.clusters.perm4 <- analyze_lists(covar.clusters, noNA[[1]], clusPhen4)

covar.clusters.rep <- list(list(list("clusters"), TRUE))
num.sig.genes.clusters <- numeric(100)
for(i in 1:100)
{
	cat("Iteration", i, "\n")
	clusters.perm.i <- sample(clusters, 657, replace=FALSE)
	clusPhen.i <- noNA[[2]]
	clusPhen.i$clusters <- clusters.perm.i
	results.clusters.perm.i <- analyze_lists(covar.clusters.rep, noNA[[1]], clusPhen.i)
	num.sig.genes.clusters[i] <- length(results.clusters.perm.i[[1]]$cutoffs_p[[1]]$clusters2)

}

save(num.sig.genes.clusters, file="numSigGeneClusters")
# num.sig.genes.clusters range from 62-3337, much lower than the 9752 (FDR < .25) found in the model:
#	gene: clusters + batch + rin + smk + py + gender + age + cancer + copd + copd*cancer

# check to see if clusters differ on RIN, BATCH
t.test(clusPhen$RIN[clusPhen$clusters==1], clusPhen$RIN[clusPhen$clusters==2])
	# no difference on RIN
	
t.test(clusPhen$BATCH[clusPhen$clusters==1], clusPhen$BATCH[clusPhen$clusters==2])

# chi square tests?
clusPhen$BATCH[clusPhen$clusters==1]
clusPhen$BATCH[clusPhen$clusters==2]


# look at gene set from Joe
# load file, variable rmeans
load("../../GeneSets/2014_03_02_sig_in_both_selected_genes_rmeans.RData")
joe.genes <- names(rmeans[rmeans > 0])
inds <- match(joe.genes, rownames(noNA[[1]])[results.clusters[[2]]$cutoffs_fdr[[1]]$clusters2])
inds[!is.na(inds)]


#clusters.genes <- rownames(noNA[[1]][results.clusters[[2]]$cutoffs_fdr[[10]]$clusters2,])

#overlap.genes <- clusters.genes[match(joe.genes, clusters.genes)]
#overlap.genes <- overlap.genes[!is.na(overlap.genes)]


results.joe <- analyze_lists(covar.clusters, joe.data, 
# look at gene set(s) from Grant
