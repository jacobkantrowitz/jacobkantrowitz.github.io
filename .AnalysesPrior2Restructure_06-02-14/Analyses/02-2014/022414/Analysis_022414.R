# Allegro Analysis - finding subgroups within the data and covariates to explain them
# Jake Kantrowitz
# 02-20-14

init_dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
setwd(init_dir)
load("data_phenotype_updated.RData")


#cat(class(data_phenotype$CANCER),"\n")
# Questions:
# 1. should I use filtered data?
#	I've shown that using cv or cv_mn filtered data is effective
#	However, only need to find subgroups, don't actually care about genes


# set covariates to use for heatmap residuals
covars_res <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE"), FALSE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE)
				)

pres <- !(is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY))
cat(pres,"\n")
# run the model for residuals for heatmap viewing
#results_res <- analyze_lists(covars_res, data_filtered$cv_mn[,pres], data_phenotype[pres,])

# store the residuals for easy access in writing heatmaps
#resids <- results_res[[1]]$model_resids


# set covariates for basic model to look for subgroups
# likely will need to look at different methods of viewing the data to find the best way
# to make the groups
covars1 <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), FALSE),
				list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), TRUE)
				)

#results_mod <- analyze_lists(covars1, data_filtered$cv_mn[,pres], data_phenotype[pres,])
results_mod3 <- analyze_lists(covars1, data_filtered$cv_mn, data_phenotype)

covars2 <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), TRUE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER"), TRUE)
				)
#results_mod2 <- analyze_lists(covars2, data_filtered$cv_mn[,pres], data_phenotype[pres,])

# find the COPD*CANCER p or fdr cutoff giving a few hundred genes
# store the indices
#inds <- results_mod[[1]]$cutoffs_p[[2]]$'COPD1:CANCER1'

# define healthy patients and create heatmap of interaction genes
#hh <- data_phenotype$HEALTHY==TRUE
#pdf("intrxGenesPY_healthy.pdf")
#generate_heatmap(inds,resids[,hh], data_phenotype[hh,], colClus="unsupervised")
#dev.off()

# define COPD only patients and create heatmap of interaction genes
#copd <- data_phenotype$CANCER==0 & data_phenotype$COPD==1
#pdf("intrxGenesPY_copd.pdf")
#generate_heatmap(inds,resids[,copd], data_phenotype[copd,], colClus="unsupervised")
#dev.off()

# define Cancer only patients and create heatmap of interaction genes
#cancer <- data_phenotype$CANCER==1 & data_phenotype$COPD==0
#pdf("intrxGenesPY_cancer.pdf")
#generate_heatmap(inds,resids[,cancer], data_phenotype[cancer,], colClus="unsupervised")
#dev.off()

# define patients with both and create heatmap of interaction genes
#cc <- data_phenotype$CANCER==1 & data_phenotype$COPD==1
#pdf("intrxGenesPY_copdca.pdf")
#generate_heatmap(inds,resids[,cc], data_phenotype[cc,], colClus="unsupervised")
#sdev.off()

# generate a heatmap
#pdf("Exploratory_copd+cancer_SV_p01_312genes_BRSGASSS_CVMN.pdf")
#generate_heatmap(inds_cvmn_sv_p01, res_sv, data_phenotype_COPDCA_sva, colClus=data_phenotype_COPDCA_sva$indicator)
#dev.off()



