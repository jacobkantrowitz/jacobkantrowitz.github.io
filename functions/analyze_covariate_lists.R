# analyze_covariate_lists accepts lists of covariates that can be analyzed using 'return_lm'
# the covariates are checked for presence/absence with 'check_covariates'
# the covariates are also checked for NAs that would prevent the models from running properly
# analyze_covariate_lists takes an expressionSet object and a list of lists of covariates
# that are represented in the phenotype data of the expressionSet object
# analyze custom lists of covariates
analyze_covariate_lists <- function(exprData, covariates)
{
	results <- 0
	
	test <- check_covariates(exprData, covariates)
	if(test==0)
	{
		presentNAs <- printNACovariates(exprData, covariates)
	}
	else
	{
		cat("\nSome of the covariates provided were not found in the phenotype data.\n",
			"Please edit the indicated covariates and resubmit.\n\n")

		return(1)
	}
	
	# if all covariates check out OK --> proceed with analysis
	if(presentNAs==0)
	{
		# if 'covariates' is a single list of covariates, run lm on data_to_analyze
		# using the covariates from the list
		if(is.character(covariates[[1]]))
		{
			model <- paste(covariates, collapse=" + ")
			#cat("Running one model:", model, "\n")
			results <- return_lm(exprData, covariates)
			return(results)
		}
		
		# else the first item of covariates is another list i.e. covariates is a list of lists
		else if(is.list(covariates[[1]]))
		{
			results <- list()
			cat("Multiple models to run:", length(covariates), "models\nRunning Models\n\n")
			for(model in 1:length(covariates))
			{
				cat("Running model #", model, "\n", sep="")
				temp <- return_lm(exprData, covariates[[model]])
				results[[paste(covariates[[model]], collapse="+")]] <- temp
				cat("Done running model #", model, "\n\n", sep="")
			}
		}
	}
	
	# else some covariates are not available --> cannot create all models
	else
	{
		cat("\nSome of the covariates provided were found to have NAs in the phenotype data.\n",
		"Please remove the covariates or their NAs and resubmit.\n\n")
	
	}
	
	return(results)
}