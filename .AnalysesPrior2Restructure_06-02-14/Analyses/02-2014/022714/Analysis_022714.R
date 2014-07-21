# Allegro Analysis - using CLIA/IDE Study Section Status as a covariate
# Jake Kantrowitz
# 02-27-14

init.dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("data_phenotype_updated.RData")
load("barcode.section.RData")
cat("Loaded OK\n")
setwd(init.dir)

data_phenotype$StudySection <- barcode.section$'Study Section'[match(data_phenotype$BARCODE, barcode.section$Barcode)]

# set covariates to use for analysis
# analyze only cancer patients - only non-COPD patients
covars_ss <- list(list(list("StudySection"), FALSE),
				list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE"), TRUE)
				)

covars_ccc <- list(list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE),
					list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD"), TRUE),
					list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "CANCER", "COPD", "COPD*CANCER"), TRUE)
				)
				
covars.copd <- list(list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE)
				)

covars.cancer <- list(list(list("BATCH", "RIN", "StudySection", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE)
				)




ind.NA <- which((is.na(data_phenotype$CANCER) | is.na(data_phenotype$PY)))
noNA <- remove_patients(ind.NA, data_expression, data_phenotype)
noNA.cv <- remove_patients(ind.NA, data_filtered$cv, data_phenotype)
noNA.cvmn <- remove_patients(ind.NA, data_filtered$cv_mn, data_phenotype)

# run some models
resultsUNF_ss <- analyze_lists(covars, noNA[[1]], noNA[[2]])

resultsUNF_ccc <- analyze_lists(covars_ccc, noNA[[1]], noNA[[2]])

resultsUNF_res <- analyze_lists(covars_ss, noNA[[1]], noNA[[2]])


resultsCV_ss <- analyze_lists(covars_ss, noNA.cv[[1]], noNA.cv[[2]])
resultsCV_ccc <- analyze_lists(covars_ccc, noNA.cv[[1]], noNA.cv[[2]])

#######################################################
########### run only COPD vs healthy patients
ind.CANCER.cv <- which(noNA.cv[[2]]$CANCER==1)
noCANCER.cv <- remove_patients(ind.CANCER.cv, noNA.cv[[1]], noNA.cv[[2]])
resultsCV.copd <- analyze_lists(covars.copd, noCANCER.cv[[1]], noCANCER.cv[[2]])

resultsCV.copd.resid <- analyze_lists(covars_ss, noCANCER.cv[[1]], noCANCER.cv[[2]])

# run COPD genes for Cancer
# find genes loosely significant for COPD @ p < 0.05
genes.copd.p.05.ind <- resultsCV.copd[[1]]$cutoffs_p[[1]]$COPD1
# filter genes not significant for COPD @ p < 0.05
genes.cv.copd <- keep_genes(genes.copd.p.05.ind, data_filtered$cv)


# run cancer model on genes significant for COPD @ p < 0.05
resultsCV.copd.cancer <- analyze_lists(covars.cancer

########### run only Cancer vs healthy patients
ind.COPD.cv <- which(noNA.cv[[2]]$COPD==1)
noCOPD.cv <- remove_patients(ind.COPD.cv, noNA.cv[[1]], noNA.cv[[2]])
resultsCV.cancer <- analyze_lists(covars.cancer, noCOPD.cv[[1]], noCOPD.cv[[2]])

# run Cancer genes for COPD
# find genes loosely significant for cancer @ p < 0.05
genes.cancer.p.05.ind <- resultsCV.cancer[[1]]$cutoffs_p[[1]]$CANCER1
# filter genes not significant for cancer @ p < 0.05
genes.cv.cancer <- keep_genes(genes.cancer.p.05.ind, data_filtered$cv)
# run copd model on genes significant for cancer @ p < 0.05

#######################################################
# generate some heatmaps

generate_heatmap(resultsUNF_ccc[[4]]$cutoffs_p[[3]]$'CANCER1:COPD1', resultsUNF_res[[1]]$model_resids, noNA[[2]])

generate_heatmap(resultsUNF_ccc[[4]]$cutoffs_p[[2]]$'CANCER1:COPD1', resultsUNF_res[[1]]$model_resids, noNA[[2]])



#compute fold change for COPD genes to then map them
library(gtools)

