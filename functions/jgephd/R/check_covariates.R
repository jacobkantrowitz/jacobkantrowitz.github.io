# check_covariates determines whether all of the given covariates exist within the 
# phenotype data (pData) of the given expressionSet object
# check_covariates takes an expressionSet object and a list of covariates
# check_covariates returns a 1/0 to indicate true/false (why 1/0 and not TRUE/FALSE?)

# editing this function to also perform the NA checking for covariates
# this function will now (as of 2014-06-04) return FAIL if either of these two conditions are met
#   1. one of the phenotypes listed in covariates is not present
#   2. one of the phenotypes listed has NAs present in the given exprData object
check_covariates <- function(exprData, covariates)
{
	FAIL = 1
	PASS = 0
	failed <- numeric()
	
	# if the supplied 'covariates' list represents one model
	if(is.character(covariates[[1]]))
	{
		isolate_covariates <- strsplit(as.character(covariates), "*", fixed=TRUE)
		
		all_covariates <- character()
		for(covs in isolate_covariates)
		{
			for(cov in covs)
				all_covariates <- append(all_covariates, cov)
		}
		
		all_covariates <- unique(all_covariates)
		cat("Unique covariates:", all_covariates, "\n")
		
		# each char item in 'covariates' must be from the names of 'phenotype'
		covs_available <- number_of(is.na(match(all_covariates, varLabels(exprData))))
		# if at least one supplied covariate is not in the phenotype data
		if(covs_available)
		{
			cat("Some covariates are not traits or are spelled incorrectly:\n",
					as.character(all_covariates[is.na(match(all_covariates, varLabels(exprData)))]),
					"\n")
			model <- paste(covariates, collapse=" + ")
			cat("Model to edit:", model, "\n\n")
			return(FAIL)
				
		}
		# else all the supplied covariates are in the phenotype data
		else
		{
			model <- paste(covariates, collapse=" + ")
			cat("Model OK:", model, "\n\n")

      # now also check that each of the covariates (that we've shown are present)
      # do not have NAs in them which would cause problems when running the linear model
			presentNAs <- printNACovariates(exprData, all_covariates)
			if(presentNAs==0)
			{
		  	return(PASS)
			}
      else # there are NAs in the data that will cause problems
      {
        return(FAIL)
      }
		}
	}
	
	# else if the supplied 'covariates' list represents multiple models
	else if(is.list(covariates[[1]]))
	{
		cat("\nMultiple models to check:", length(covariates), "models\nChecking Models\n\n")
		for(covs in 1:length(covariates))
		{
			cat("Checking model #", covs, ":", paste(covariates[[covs]], collapse=" + "), "\n", sep="")
			test <- check_covariates(exprData, covariates[[covs]])
			
			if(test)
			{
				failed <- append(failed, covs)
			}
			#cat("Done check model #", covs, "\n\n")
		}
		
		if(length(failed)==0)
		{
			cat("All models passed: OK\n")
			return(PASS)
		}
		else
		{
			cat("\nModels failed: ", length(failed), "/", length(covariates), "\n\n", sep="")
			cat("Edit models: ", paste(failed, collapse=", "), "\n\n", sep="")
			return(FAIL)
		}
	}
	
	# else the first item of covariates is neither a char nor a list i.e. wrong data type
	else
	{
		cat("\nERROR: The covariates data supplied is not of the right class.\n")
		cat("Covariates of class", class(covariates), "\n\n")
		return(FAIL)
	}

}