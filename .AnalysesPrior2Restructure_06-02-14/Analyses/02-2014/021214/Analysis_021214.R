# Allegro Analysis
# Jake Kantrowitz
# 02-12-14

init_dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
setwd(init_dir)

# Find 




covars_copd <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE"), TRUE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER"), TRUE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), TRUE),
				list(list("BATCH", "RIN", "GENDER", "AGE"), FALSE),
				list(list("COPD", "CANCER", "COPD*CANCER"), FALSE),
				list(list("SMK", "COPD", "CANCER", "COPD*CANCER"), FALSE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "indicator"), FALSE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD"), FALSE)
				)




results_cov_copd <- analyze_lists(covars_copd, data_expression, data_phenotype)

results_cov_filtCV_copd <- analyze_lists(covars_copd, data_filtered$cv, data_phenotype)

#results_cov_filtCVMN_copd <- analyze_lists(covars_copd, data_filtered$cv_mn, data_phenotype)




