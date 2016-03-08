# calculate_total_levels returns the number of different covariate levels in a given set
# of covariates. Every level of each factor is counted; interaction terms are counted; 
# individual continuous covariates are counted 
calculate_total_levels <- function(exprData, covariates)
{
	total <- 0
	for(cov in covariates)
	{
		if(grepl("*", cov, fixed=TRUE))
		{
			total <- total + 1
		}
		else
		{
			col <- match(cov, varLabels(exprData))
			if(is.factor(pData(exprData)[,col]))
			{
				temp <- length(levels(pData(exprData)[,col])) - 1
				#print(temp)
				total <- total + temp
			}
			else
			{
				total <- total + 1
			}
		}
	}
	return(total)
}