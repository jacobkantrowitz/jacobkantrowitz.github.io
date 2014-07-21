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
						
#results.cancer <- analyze_lists(covars.cancer, cancer.CVMN[[1]], cancer.CVMN[[2]])
#results.can.gender0 <- analyze_lists(covars.can.gender, can.gender0.CVMN[[1]], can.gender0.CVMN[[2]])
#results.can.gender1 <- analyze_lists(covars.can.gender, can.gender1.CVMN[[1]], can.gender1.CVMN[[2]])

# try removing the covariates/confounds
covars.cancer1 <- list(list(list("CANCER", "GENDER", "CANCER*GENDER"),TRUE))
#results.cancer1 <- analyze_lists(coavrs.cancer1, cancer.CVMN[[1]], cancer.CVMN[[2]])


# some helpful examples of running covariate list models, saving genes (txt and heatmap) and returning clusters
#results.UNF <- analyze_lists(covars, data.noNA[[1]], data.noNA[[2]])
#save.genes(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], "cvmnFILT.noNA.smk+copd.p.01")
#cs <- return_cluster(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], colClus="unsupervised", type=ROWS)