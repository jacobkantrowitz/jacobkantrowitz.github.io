# Allegro Analysis - preparing for the SILCC call
# Jake Kantrowitz
# 02-12-14

init_dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
setwd(init_dir)



covars_copd <- list(list(list("BATCH", "RIN", "SMK", "GENDER", "AGE"), TRUE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER"), TRUE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"), TRUE),
				list(list("BATCH", "RIN", "GENDER", "AGE"), FALSE),
				list(list("COPD", "CANCER", "COPD*CANCER"), FALSE),
				list(list("SMK", "COPD", "CANCER", "COPD*CANCER"), FALSE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "indicator"), FALSE),
				list(list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD"), FALSE)
				)
				
covars_sv <- list(list(list("BATCH" ,"RIN", "SMK", "GENDER", "AGE", "V12", "V13", "V14"), FALSE),
				list(list("BATCH" ,"RIN", "SMK", "GENDER", "AGE", "CANCER", "COPD", "V12", "V13", "V14"), FALSE),
				list(list("BATCH" ,"RIN", "SMK", "GENDER", "AGE", "CANCER", "COPD", "CANCER*COPD", "V12", "V13", "V14"),FALSE),
				list(list("BATCH" ,"RIN", "SMK", "GENDER", "AGE"),TRUE),
				list(list("V12", "V13", "V14", "V15"), FALSE)
				)



#results_cov_copd <- analyze_lists(covars_copd, data_expression, data_phenotype)

#results_cov_filtCV_copd <- analyze_lists(covars_copd, data_filtered$cv, data_phenotype)

#results_cov_filtCVMN_copd <- analyze_lists(covars_copd, data_filtered$cv_mn, data_phenotype)

#results_cov_sv_filtCVMN <- analyze_lists(covars_sv, data_filtered$cv_mn, data_phenotype_COPDCA_sva)
results_cov_sv_filtCVMN2 <- analyze_lists(covars_sv, data_filtered$cv_mn, data_phenotype_COPDCA_sva)


### for results_cov_sv_filtCVMN
# SVA built with COPD + CANCER + with and without COPD*CANCER

# 312 genes p < 0.01
# 11074 genes run; 111 by chance
inds_cvmn_sv_p01 <- results_cov_sv_filtCVMN[[2]]$cutoffs_p[[2]]$'CANCER1:COPD1'

# 177 genes p < 0.005
# 11074 genes run; 56 by chance
inds_cvmn_sv_p005 <- results_cov_sv_filtCVMN[[2]]$cutoffs_p[[3]]$'CANCER1:COPD1'

# 28 genes fdr < 0.2
# 11074 genes run
inds_cvmn_sv_f2 <- results_cov_sv_filtCVMN[[2]]$cutoffs_fdr[[1]]$'CANCER1:COPD1'

res_sv <- results_cov_sv_filtCVMN2[[1]]$model_resids

pdf("Exploratory_copd+cancer_SV_p01_312genes_BRSGASSS_CVMN.pdf")
generate_heatmap(inds_cvmn_sv_p01, res_sv, data_phenotype_COPDCA_sva, colClus=data_phenotype_COPDCA_sva$indicator)
dev.off()

pdf("Exploratory_copd+cancer_SV_p005_177genes_BRSGASSS_CVMN.pdf")
generate_heatmap(inds_cvmn_sv_p005, res_sv, data_phenotype_COPDCA_sva, colClus=data_phenotype_COPDCA_sva$indicator)
dev.off()

pdf("Exploratory_copd+cancer_SV_f2_28genes_BRSGASSS_CVMN.pdf")
generate_heatmap(inds_cvmn_sv_f2, res_sv, data_phenotype_COPDCA_sva, colClus=data_phenotype_COPDCA_sva$indicator)
dev.off()


### for results_cov_copd and results_cov_filtCVMN_copd
#generate heatmap of ~300 genes with residuals from model 1 or model 2 using copd*cancer genes from model 3

# generate heatmap with unfiltered models

#summary(results_cov_copd[[3]]$cutoffs_p[[1]])
# 299 genes with p < 0.01 where 197 expected by chance
#inds_unf_p01 <- results_cov_copd[[3]]$cutoffs_p[[2]]$'COPD1:CANCER1'
# 171 genes with p < 0.005 where 99 expected by chance
#inds_unf_p005 <- results_cov_copd[[3]]$cutoffs_p[[3]]$'COPD1:CANCER1'
#res1_unf <- results_cov_copd[[1]]$model_resids

# generate heatmap with residuals from BATCH + RIN + SMK + GENDER + AGE
#pdf("Exploratory_copd+cancer_p01_291genes_BRSGA_unfilt.pdf")
#generate_heatmap(inds_unf_p01, res1_unf, data_phenotype, colClus=data_phenotype$indicator)
#dev.off()

# generate heatmap with residuals from BATCH + RIN + SMK + GENDER + AGE
#pdf("Exploratory_copd+cancer_p005_171genes_BRSGA_unfilt.pdf")
#generate_heatmap(inds_unf_p005, res1_unf, data_phenotype, colClus=data_phenotype$indicator)
#dev.off()


# generate heatmap with CV filtered models


# DONE
# generate heatmap with CV_MEAN filtered models
#check for the numbers of significant genes at varying p/fdr values
#summary(results_cov_copd[[3]]$cutoffs_p[[1]])

#summary(results_cov_filtCV_copd[[3]]$cutoffs_p[[1]])

# summary(results_cov_filtCVMN_copd[[3]]$cutoffs_p[[1]])
# above summary shows that 219 genes significant at p <0.01 for COPD*CANCER
# only expect 11074 * 0.01 = 111 genes by chance

#inds_cvmn <- results_cov_filtCVMN_copd[[3]]$cutoffs_p[[2]]$'COPD1:CANCER1'
#res1_cvmnv <- results_cov_filtCVMN_copd[[1]]$model_resids
#res2_cvmnv <- results_cov_filtCVMN_copd[[2]]$model_resids

# generate heatmap with residuals from BATCH + RIN + SMK + GENDER + AGE
#pdf("Exploratory_copd*cancer_p01_219genes_BRSGA_CVMN.pdf")
#generate_heatmap(inds_cvmn, res1_cvmnv, data_phenotype, colClus=data_phenotype$indicator)
#dev.off()

# generate heatmap with residuals from BATCH + RIN + SMK + GENDER + AGE + CANCER + COPD
#pdf("Exploratory_copd*cancer_p01_219genes_BRSGACC_CVMNx.pdf")
#generate_heatmap(inds_cvmn, res2_cvmnv, data_phenotype, colClus=data_phenotype$indicator)
#dev.off()