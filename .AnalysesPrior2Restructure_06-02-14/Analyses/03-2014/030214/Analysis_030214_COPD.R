# Allegro Analysis - finding subgroups within the data and covariates to explain them
# Find groups based on the splits in the interaction genes
# Next, find genes significant for those groups and see what the genes are
# Jake Kantrowitz
# 02-28-14

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


# calculate fold change
# already working with log data
# thus, fold change calculated by subtraction (dividing non-log transformed values)
# i.e. mean(Group 1) - mean(Group 2)
# generally take the absolute value
# fold change = | mean(Group 1) - mean(Group 2)|

fold.change <- function(data, groups)
{
	g1.ind <- which(groups==1)
	g2.ind <- which(groups==2)
	num.genes <- nrow(data)
	fcs <- matrix(nrow=num.genes, ncol=3)
	
	for(gene in 1:num.genes)
	{
		fcs[gene, 1] <- mean(data[gene, g1.ind])
		fcs[gene, 2] <- mean(data[gene, g2.ind])
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
