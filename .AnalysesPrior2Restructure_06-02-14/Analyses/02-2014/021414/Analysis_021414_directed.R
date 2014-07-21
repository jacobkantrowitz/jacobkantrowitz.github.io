# Allegro Analysis - preparing for the SILCC call
# Jake Kantrowitz
# 02-12-14

init_dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
setwd(init_dir)


# set covariates to use for heatmap residuals
covars_res <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE"), TRUE)
				)

# looking at 91/98 genes from Katie's signature
# modeled below covariate list
# 20-30 genes significant for COPD*CANCER
# heatmap looks the same with split groups within the diseased patients etc.
# really need to look at what these groups are
covars_katie <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), TRUE)
				)
				
				
results_res <- analyze_lists(covars_res, data_expression, data_phenotype)

resids <- results_res[[1]]$model_resids

# set covariates to use in initial COPD-specific model
covars_copd <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD"), TRUE)
				)
				
# set the data to include only healthy participants and participants with COPD only (no cancer)
pts1 <- data_phenotype$CANCER==0
expr_copdHealthy <- data_expression[, pts1]
phen_copdHealthy <- data_phenotype[pts1,]

# analyze the data with the above covariates
#results_copd <- analyze_lists(covars_copd, expr_copdHealthy, phen_copdHealthy)

# find the genes 'significant' for COPD at p < 0.05 (likely includes false positives)
# this is meant as a loose filter for COPD-specific genes
inds_copd <- results_copd[[1]]$cutoffs_p[[1]]$COPD1	

# set covariates to use in second step analysis: COPD-specific genes involved in lung carcinogenesis
covars_copdca <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "CANCER"), TRUE)
				)
				
# set the data to include only patients with COPD only and with cancer + COPD
# set the data to include only those genes filtered from above
pts2 <- (data_phenotype$CANCER==1 & data_phenotype$COPD==1) | (data_phenotype$CANCER==0 & data_phenotype$COPD==1)
expr_copdCancer <- data_expression[inds_copd,pts2]
phen_copdCancer <- data_phenotype[pts2,]

#results_copdca <- analyze_lists(covars_copdca, expr_copdCancer, phen_copdCancer)
			
				
#inds_cvmn_sv_p01 <- results_cov_sv_filtCVMN[[2]]$cutoffs_p[[2]]$'CANCER1:COPD1'
#res_sv <- results_cov_sv_filtCVMN2[[1]]$model_resids

#pdf("Exploratory_copd+cancer_SV_p01_312genes_BRSGASSS_CVMN.pdf")
#generate_heatmap(inds_cvmn_sv_p01, res_sv, data_phenotype_COPDCA_sva, colClus=data_phenotype_COPDCA_sva$indicator)
#dev.off()



