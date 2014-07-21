# COPD*CANCER interaction analysis in ~700 patients from Allegro Cohort (IDE & CLIA)
# Jake Kantrowitz
# kantro@bu.edu, jacob.kantrowitz@gmail.com
# 03-24-14

######################################################################

init.dir <- getwd()
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Functions_v01.r")
source("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/data_phenotype_updated.RData")
cat("Loaded OK\n")

cat("\nLoading cluster data for 657 patients; removed 43 w/o cancer or PY values\n")
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Analyses/030214/clusters.RData")
load("/net/ibx109/ibx109fs1/protected/projects/pulmarray/Allegro/COPD_Cancer/Analyses/031114/clusPhen.RData")
cat("Loaded OK\n")

save.packages()

## (1) remove patients without cancer designation (i.e. cancer==na)
noNA.ind <- union(which(is.na(data_phenotype$CANCER)), which(is.na(data_phenotype$PY)))
noNA.UNF <- remove_patients(noNA.ind, data_expression, data_phenotype)
noNA.CV <- remove_patients(noNA.ind, data_filtered$cv, data_phenotype)
noNA.CVMN <- remove_patients(noNA.ind, data_filtered$cv_mn, data_phenotype)

# working with WGCNA to develop modules within the Allegro micro arrays
# load WGCNA package
# will likely NEED TO CHANGE VERSION OF R TO GET WGCNA TO LOAD
shared.library.path <- file.path("/unprotected/projects/cbmhive/R_packages", getRversion());
.libPaths(shared.library.path, .libPaths());
library(WCGNA)
# the following setting is important, do not omit (for WGCNA work)
options(stringsAsFactors = FALSE)

# We transpose the data so that the rows correspond to samples and the columns correspond to genes
tran.noNA.UNF <- t(noNA.UNF[[1]])

collectGarbage()

# Following code from "FemaleLiver-02-networkConstr-blockwise.pdf" tutorial from WGCNA website
# 2.c automatic block-wise network construction and module detection
# 2.c.1 choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from=12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(tran.noNA.UNF, powerVector=powers, verbose=5)








covars <- list(list(list("RIN", "SMK","GENDER", "AGE"), FALSE),
				list(list("RIN", "SMK", "GENDER", "AGE", "CANCER", "COPD"), FALSE),
				list(list("RIN", "SMK", "GENDER", "AGE", "CANCER", "COPD", "CANCER*COPD"), FALSE),
				list(list("RIN", "SMK", "GENDER", "AGE", "CANCER", "COPD", "SMK*COPD"), TRUE)
				)

# some helpful examples of running covariate list models, saving genes (txt and heatmap) and returning clusters
#results.UNF <- analyze_lists(covars, data.noNA[[1]], data.noNA[[2]])
#save.genes(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], "cvmnFILT.noNA.smk+copd.p.01")
#cs <- return_cluster(results.UNF[[1]]$cutoffs_p[[2]]$'SMK2:COPD1', data.noNA[[1]], data.noNA[[2]], colClus="unsupervised", type=ROWS)