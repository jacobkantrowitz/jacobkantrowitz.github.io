# Allegro Analysis - finding subgroups within the data and covariates to explain them
# Jake Kantrowitz
# 02-20-14

init.dir <- getwd()
setwd("../../")
source("Allegro_Analysis_v02.r")
cat("\nLoading updated phenotype information\n")
load("data_phenotype_updated.RData")
cat("Loaded OK\n")
setwd(init.dir)


# set covariates to use for analysis
# analyze only cancer patients - only non-COPD patients
covars <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "COPD"), TRUE),
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

# want comparisons for 1:3, 2:3, 3:4
cancer.unfiltered <- list(); copd.unfiltered <- list(); cancopd.unfiltered <- list()
cancer.cv <- list(); copd.cv <- list(); cancopd.cv <- list()
cancer.cvmn <- list(); copd.cvmn <- list(); cancopd.cvmn <- list()

numCV <- 14764
numCVMN <- 11074
for(mod in 1:length(resultsUNF[[4]]$models))
{
	cancer.unfiltered[[mod]] <- anova(resultsUNF[[1]]$models[[mod]], resultsUNF[[3]]$models[[mod]])
	copd.unfiltered[[mod]] <- anova(resultsUNF[[2]]$models[[mod]], resultsUNF[[3]]$models[[mod]])
	cancopd.unfiltered[[mod]] <- anova(resultsUNF[[3]]$models[[mod]], resultsUNF[[4]]$models[[mod]])

	if(mod <=numCV)
	{
		cancer.cv[[mod]] <- anova(resultsCV[[1]]$models[[mod]], resultsCV[[3]]$models[[mod]])
		copd.cv[[mod]] <- anova(resultsCV[[2]]$models[[mod]], resultsCV[[3]]$models[[mod]])
		cancopd.cv[[mod]] <- anova(resultsCV[[3]]$models[[mod]], resultsCV[[4]]$models[[mod]])
	
		if(mod <= numCVMN)
		{
			cancer.cvmn[[mod]] <- anova(resultsCVMN[[1]]$models[[mod]], resultsCVMN[[3]]$models[[mod]])
			copd.cvmn[[mod]] <- anova(resultsCVMN[[2]]$models[[mod]], resultsCVMN[[3]]$models[[mod]])
			cancopd.cvmn[[mod]] <- anova(resultsCVMN[[3]]$models[[mod]], resultsCVMN[[4]]$models[[mod]])
		}
	}
}


# pull out F-test p values
cancer.unfiltered.p <- numeric(); copd.unfiltered.p <- numeric(); cancopd.unfiltered.p <- numeric()
cancer.cv.p <- numeric(); copd.cv.p <- numeric(); cancopd.cv.p <- numeric()
cancer.cvmn.p <- numeric(); copd.cvmn.p <- numeric(); cancopd.cvmn.p <- numeric()


for(i in 1:length(resultsUNF[[4]]$models))
{
	cancer.unfiltered.p[i] <- cancer.unfiltered[[i]]$'Pr(>F)'[[2]]
	copd.unfiltered.p[i] <- copd.unfiltered[[i]]$'Pr(>F)'[[2]]
	cancopd.unfiltered.p[i] <- cancopd.unfiltered[[i]]$'Pr(>F)'[[2]]
	
	if(i <=numCV)
	{
		cancer.cv.p[i] <- cancer.cv[[i]]$'Pr(>F)'[[2]]
		copd.cv.p[i] <- copd.cv[[i]]$'Pr(>F)'[[2]]
		cancopd.cv.p[i] <- cancopd.cv[[i]]$'Pr(>F)'[[2]]

	}
			if(i <= numCVMN)
			{
				cancer.cvmn.p[i] <- cancer.cvmn[[i]]$'Pr(>F)'[[2]]
				copd.cvmn.p[i] <- copd.cvmn[[i]]$'Pr(>F)'[[2]]
				cancopd.cvmn.p[i] <- cancopd.cvmn[[i]]$'Pr(>F)'[[2]]

			}
}

covars_smk <- list(list(list("SMK"), TRUE))
results_smk <- analyze_lists(covars_smk, noNA[[1]], noNA[[2]])