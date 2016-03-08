# analyze_lists calls analyze_covariate_lists with a list of lists of covariates
# the difference between this function and analyze_covariate_lists is that here
# each list has a boolean (TRUE/FALSE) and the covariate list is only modeled if
# the boolean is TRUE
# analyze_lists takes a list (see below for the form of that list) and an expressionSet object
# analyze_lists returns the result of running analayze_covariate_lists on the select list
# of covariate lists accepted here with the TRUE

# NOTE: this may be out-of-date based on what I expect my workflow to look like going forward
analyze_lists <- function(listBool, exprData)
{
	# listBool is a list with n items
	# for each of n items there are two items 
	# item 1: list of covariates to be sent to analyze_covariate_lists
	# item 2: boolean indicating whether or not to send item 1
	
	# create a list of covariate lists to analyze
	analyze <- list()

	cat("Received ",  length(listBool), " models\nChecking models\n", sep="")
	# for each item i, if it's boolean is TRUE, append it to a list of covariate lists to analyze
	for(analysis in listBool)
	{
		complete <- analysis[[2]]
		covars <- analysis[[1]]

		if(complete==TRUE)
		{
			analyze <- append(analyze, list(covars))
		}
	}
	
	cat("Running ", length(analyze), " of ", length(listBool), " models\n\n", sep="")
	if(length(analyze) > 0)
	{
		return(analyze_covariate_lists(exprData, analyze))
	}
}