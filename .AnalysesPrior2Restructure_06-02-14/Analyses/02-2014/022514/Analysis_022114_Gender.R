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
covars <- list(list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER"), TRUE),
					list(list("BATCH", "RIN", "SMK", "PY", "GENDER", "AGE", "CANCER", "GENDER*CANCER"), TRUE)
				)

############## Run with unfiltered data

ind.COPD <- which(data_phenotype$COPD==1)
noCOPD <- remove_patients(ind.COPD, data_expression, data_phenotype)

ind.NA <- which((is.na(noCOPD[[2]]$CANCER) | is.na(noCOPD[[2]]$PY)))
noNA <- remove_patients(ind.NA, noCOPD[[1]], noCOPD[[2]])

resultsUNF <- analyze_lists(covars, noNA[[1]], noNA[[2]])


############## Run with CV filtered data
ind.COPD <- which(data_phenotype$COPD==1)
noCOPD <- remove_patients(ind.COPD, data_filtered$cv, data_phenotype)

ind.NA <- which((is.na(noCOPD[[2]]$CANCER) | is.na(noCOPD[[2]]$PY)))
noNA <- remove_patients(ind.NA, noCOPD[[1]], noCOPD[[2]])

resultsCV <- analyze_lists(covars, noNA[[1]], noNA[[2]])

############## Run with CVMN filtered data
ind.COPD <- which(data_phenotype$COPD==1)
noCOPD <- remove_patients(ind.COPD, data_filtered$cv_mn, data_phenotype)

ind.NA <- which((is.na(noCOPD[[2]]$CANCER) | is.na(noCOPD[[2]]$PY)))
noNA <- remove_patients(ind.NA, noCOPD[[1]], noCOPD[[2]])

resultsCVMN <- analyze_lists(covars, noNA[[1]], noNA[[2]])

# want comparisons for 1:3, 2:3, 3:4
cancer.unfiltered <- list(); gender.unfiltered <- list(); cangen.unfiltered <- list()
cancer.cv <- list(); gender.cv <- list(); cangen.cv <- list()
cancer.cvmn <- list(); gender.cvmn <- list(); cangen.cvmn <- list()

numCV <- 14764
numCVMN <- 11074
for(mod in 1:length(resultsUNF[[4]]$models))
{
	cancer.unfiltered[[mod]] <- anova(resultsUNF[[1]]$models[[mod]], resultsUNF[[3]]$models[[mod]])
	gender.unfiltered[[mod]] <- anova(resultsUNF[[2]]$models[[mod]], resultsUNF[[3]]$models[[mod]])
	cangen.unfiltered[[mod]] <- anova(resultsUNF[[3]]$models[[mod]], resultsUNF[[4]]$models[[mod]])

	if(mod <=numCV)
	{
		cancer.cv[[mod]] <- anova(resultsCV[[1]]$models[[mod]], resultsCV[[3]]$models[[mod]])
		gender.cv[[mod]] <- anova(resultsCV[[2]]$models[[mod]], resultsCV[[3]]$models[[mod]])
		cangen.cv[[mod]] <- anova(resultsCV[[3]]$models[[mod]], resultsCV[[4]]$models[[mod]])
	
		if(mod <= numCVMN)
		{
			cancer.cvmn[[mod]] <- anova(resultsCVMN[[1]]$models[[mod]], resultsCVMN[[3]]$models[[mod]])
			gender.cvmn[[mod]] <- anova(resultsCVMN[[2]]$models[[mod]], resultsCVMN[[3]]$models[[mod]])
			cangen.cvmn[[mod]] <- anova(resultsCVMN[[3]]$models[[mod]], resultsCVMN[[4]]$models[[mod]])
		}
	}
}


# pull out F-test p values
cancer.unfiltered.p <- numeric(); gender.unfiltered.p <- numeric(); cangen.unfiltered.p <- numeric()
cancer.cv.p <- numeric(); gender.cv.p <- numeric(); cangen.cv.p <- numeric()
cancer.cvmn.p <- numeric(); gender.cvmn.p <- numeric(); cangen.cvmn.p <- numeric()


for(i in 1:length(resultsUNF[[4]]$models))
{
	cancer.unfiltered.p[i] <- cancer.unfiltered[[i]]$'Pr(>F)'[[2]]
	gender.unfiltered.p[i] <- gender.unfiltered[[i]]$'Pr(>F)'[[2]]
	cangen.unfiltered.p[i] <- cangen.unfiltered[[i]]$'Pr(>F)'[[2]]
	
	if(i <=numCV)
	{
		cancer.cv.p[i] <- cancer.cv[[i]]$'Pr(>F)'[[2]]
		gender.cv.p[i] <- gender.cv[[i]]$'Pr(>F)'[[2]]
		cangen.cv.p[i] <- cangen.cv[[i]]$'Pr(>F)'[[2]]

	}
			if(i <= numCVMN)
			{
				cancer.cvmn.p[i] <- cancer.cvmn[[i]]$'Pr(>F)'[[2]]
				gender.cvmn.p[i] <- gender.cvmn[[i]]$'Pr(>F)'[[2]]
				cangen.cvmn.p[i] <- cangen.cvmn[[i]]$'Pr(>F)'[[2]]

			}
}

covars_resid <- list(list(list("BATCH", "RIN", "SMK", "AGE", "PY"), TRUE))

results_resid <- analyze_lists(covars_resid, noNA[[1]], noNA[[2]])





