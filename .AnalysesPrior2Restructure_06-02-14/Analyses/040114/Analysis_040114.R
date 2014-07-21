# COPD*CANCER interaction analysis in ~700 patients from Allegro Cohort (IDE & CLIA)
# Jake Kantrowitz
# kantro@bu.edu, jacob.kantrowitz@gmail.com
# 03-31-14

# Goal(s) for the day:
# Run SVA to account for COPD+CANCER (2 separate phenotypes)
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

## (1) remove patients without cancer designation (i.e. cancer==na)
noNA.ind <- union(which(is.na(data_phenotype$CANCER)), which(is.na(data_phenotype$PY)))
noNA.UNF <- remove_patients(noNA.ind, data_expression, data_phenotype)
noNA.CV <- remove_patients(noNA.ind, data_filtered$cv, data_phenotype)
noNA.CVMN <- remove_patients(noNA.ind, data_filtered$cv_mn, data_phenotype)


# Run SVA
# modified return_sva in Allegro_Functions:
#	mod0 = BATCH + RIN + SMK + GENDER + AGE
#	mod1 = BATCH + RIN + SMK + GENDER + AGE + COPD + CANCER

#sva.copd.cancer <- return_sva(noNA.UNF[[1]], noNA.UNF[[2]])
#save(sva.copd.cancer, file="sva.copd.cancer.RData")
load("sva.copd.cancer.RData") # generated 03-31-14

#sva.Phen <- cbind(noNA.UNF[[2]], sva.copd.cancer$sv)
#names(sva.Phen)[13:17] <- c("sv1", "sv2", "sv3", "sv4", "sv5")
#save(sva.Phen, file="sva.Phen.Rdata")
load("sva.Phen.Rdata") # generated 03-31-14


covars.sv <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3", "COPD", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3", "COPD", "CANCER", "COPD*CANCER"), TRUE))

results.sv.CVMN <- analyze_lists(covars.sv, noNA.CVMN[[1]], sva.Phen)

# save genes and heatmap of raw data FDR < 0.25
save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[1]]$'COPD1:CANCER1',noNA.CVMN[[1]], sva.Phen, "sv123.cancer.copd.fdr.25") 

# run model for residuals to use in heatmaps
covars.sv.res <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3"),TRUE))
results.sv.CVMN.res <- analyze_lists(covars.sv.res, noNA.CVMN[[1]], sva.Phen)
res <- results.sv.CVMN.res[[1]]$model_resids

# save genes and heatmap of raw data
save.genes(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1',noNA.CVMN[[1]], sva.Phen, "sv123.cancer.copd.p.01.raw") 

# save genes and heatmap of residuals
save.genes(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1',res, sva.Phen, "sv123.cancer.copd.p.01.residual") 

# define clusters for residuals heatmap (there appear to be about 2 clusters, compared with >=5 clusters in raw data)
clusters.CC.intx <- return_cluster(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1',res, sva.Phen, n.clusters=2, type=ROWS)
# save the gene lists separately for each of the two clusters
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx==1], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.residual.cluster1.txt")
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx==2], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.residual.cluster2.txt")

# define clusters for raw heatmap (there appear to be >=4 clusters)
# check clustering parameters
clusters.CC.intx.raw <- return_cluster(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1',noNA.CVMN[[1]], sva.Phen, n.clusters=4, type=ROWS)
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx.raw==1], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.raw.cluster1.txt")
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx.raw==2], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.raw.cluster2.txt")
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx.raw==3], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.raw.cluster3.txt")
save_entrez(results.sv.CVMN[[2]]$cutoffs_p[[2]]$'COPD1:CANCER1'[clusters.CC.intx.raw==4], rownames(noNA.CVMN[[1]]), "sv123.cancer.copd.p.01.raw.cluster4.txt")


# Interactive genes are difficult to find
# Let's try looking for 'additive' genes that are significant for COPD and significant for CANCER

# find COPD genes
copd.ind <- sva.Phen$COPD==1 & sva.Phen$CANCER==0 
healthy.ind <- sva.Phen$COPD==0 & sva.Phen$CANCER==0
copd.healthy.ind <- healthy.ind | copd.ind 

save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[5]]$COPD1, noNA.CVMN[[1]][,copd.healthy.ind], sva.Phen[copd.healthy.ind,], "sva123.COPD.fdr.01.raw")
save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[5]]$COPD1, res[,copd.healthy.ind], sva.Phen[copd.healthy.ind,], "sva123.COPD.fdr.01.resid")

save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[4]]$COPD1, noNA.CVMN[[1]][,copd.healthy.ind], sva.Phen[copd.healthy.ind,], "sva123.COPD.fdr.05.raw")
save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[4]]$COPD1, res[,copd.healthy.ind], sva.Phen[copd.healthy.ind,], "sva123.COPD.fdr.05.resid")

clusters.copd.resid <- return_cluster(results.sv.CVMN[[2]]$cutoffs_fdr[[4]]$COPD1, res[,copd.healthy.ind], sva.Phen[copd.healthy.ind,])
save_entrez(results.sv.CVMN[[2]]$cutoffs_fdr[[4]]$COPD1[clusters.copd.resid==1], rownames(res), "sv123.copd.fdr.05.res.cluster1.txt")
save_entrez(results.sv.CVMN[[2]]$cutoffs_fdr[[4]]$COPD1[clusters.copd.resid==2], rownames(res), "sv123.copd.fdr.05.res.cluster2.txt")

# Just for the hell of it using Jake's IDE chart review to change the COPD status of several patients
#"1-16-0594" "1-15-0516" "1-11-0529" "1-11-0518" "1-7-0533"  "1-7-0502" "1-5-0576" 
# 1-15-0516 has already been removed in the 657 patients in sva.Phen because there was no CANCER status
emphysema <- c("1-16-0594", "1-15-0516", "1-11-0529", "1-11-0518", "1-7-0533", "1-7-0502", "1-5-0576", "1-8-0501",
				"1-21-0528", "1-21-0520", "1-21-0515", "1-21-0510", "1-20-0513", "1-11-0515")
pdx.dk <- c("1-25-0514", "1-22-0508", "1-16-0513")
				 
temp.Phen <- sva.Phen
temp.Phen$COPD[match(emphysema, temp.Phen$BARCODE)] <- 1
temp.Phen$COPD[match(pdx.dk, temp.Phen$BARCODE)] <- NA

# need to remove those patients with NA for COPD
noNA.ind2 <- which(is.na(temp.Phen$COPD))
noNA.UNF2 <- remove_patients(noNA.ind2, noNA.UNF[[1]], temp.Phen)
noNA.CV2 <- remove_patients(noNA.ind2, noNA.CV[[1]], temp.Phen)
noNA.CVMN2 <- remove_patients(noNA.ind2, noNA.CVMN[[1]], temp.Phen)

# Before re-running any analyses just re-print the heatmap with the changed COPD status
copd.ind2 <- temp.Phen$COPD==1 & temp.Phen$CANCER==0 
healthy.ind2 <- temp.Phen$COPD==0 & temp.Phen$CANCER==0
copd.healthy.ind2 <- healthy.ind2 | copd.ind2 
save.genes(results.sv.CVMN[[2]]$cutoffs_fdr[[5]]$COPD1, res[,copd.healthy.ind2], temp.Phen[copd.healthy.ind2,], "EDITsva123.COPD.fdr.01.resid")
# heatmaps turn out similarly (almost exact) because most of the altered COPD status patients also have cancer (only 1 COPD only, no cancer pt.)



# rerun analysis
covars.sv.edit <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3"), TRUE),
					list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3", "COPD", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "sv1", "sv2", "sv3", "COPD", "CANCER", "COPD*CANCER"), TRUE))

# remember that the Surrogate Variables included in these models were calculated based on the original non-edited COPD values
results.sv.CVMN.edit <- analyze_lists(covars.sv.edit, noNA.CVMN2[[1]], noNA.CVMN2[[2]])

res2 <- results.sv.CVMN.edit[[1]]$model_resids
copd.ind3 <- noNA.CVMN2[[2]]$COPD==1 & noNA.CVMN2[[2]]$CANCER==0
healthy.ind3 <- noNA.CVMN2[[2]]$COPD==0 & noNA.CVMN2[[2]]$CANCER==0
copd.healthy.ind3 <- copd.ind3 | healthy.ind3
# save residual data for FDR < 0.01 for COPD in model 3 including copd*cancer term
save.genes(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[5]]$COPD1, res2[,copd.healthy.ind3], noNA.CVMN2[[2]][copd.healthy.ind3,], "EDIT.reanalyze.sva123.COPD.fdr.01.resid")
# save raw data for FDR < 0.01 for COPD in model 3 including copd*cancer term
save.genes(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[5]]$COPD1, noNA.CVMN2[[1]][,copd.healthy.ind3], noNA.CVMN2[[2]][copd.healthy.ind3,], "EDIT.reanalyze.sva123.COPD.fdr.01.raw")


cancer.ind3 <- noNA.CVMN2[[2]]$CANCER==1 & noNA.CVMN2[[2]]$COPD==0
cancer.healthy.ind3 <- cancer.ind3 | healthy.ind3
# save raw data for FDR < 0.1 for Cancer in model 3, which includes an interaction term
save.genes(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[3]]$CANCER1, noNA.CVMN2[[1]][,cancer.healthy.ind3], noNA.CVMN2[[2]][cancer.healthy.ind3,], "EDIT.reanalyze.sva123.CANCER.fdr.1.raw")
save.genes(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[3]]$CANCER1, res2[,cancer.healthy.ind3], noNA.CVMN2[[2]][cancer.healthy.ind3,], "EDIT.reanalyze.sva123.CANCER.fdr.1.resid")

clusters.cancer.resid.fdr.1 <- return_cluster(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[3]]$CANCER1, res2[,cancer.healthy.ind3], noNA.CVMN2[[2]][cancer.healthy.ind3,])
save_entrez(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[3]]$CANCER1[clusters.cancer.resid.fdr.1==1], rownames(res2), "EDIT.reanalyze.sv123.cancer.fdr.1.res.cluster1.txt")
save_entrez(results.sv.CVMN.edit[[3]]$cutoffs_fdr[[3]]$CANCER1[clusters.cancer.resid.fdr.1==2], rownames(res2), "EDIT.reanalyze.sv123.cancer.fdr.1.res.cluster2.txt")


# some helpful examples of running covariate list models, saving genes (txt and heatmap) and returning clusters
#results.UNF <- analyze_lists(covars, data.noNA[[1]], data.noNA[[2]])
#save.genes(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], "cvmnFILT.noNA.smk+copd.p.01")
#cs <- return_cluster(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], colClus="unsupervised", type=ROWS)