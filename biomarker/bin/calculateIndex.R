# Author: Joshua Campbell
# Date: Unknown
# Purpose: Script to calculate the indexes for the biomarker pipeline (the order in which models are run)

## Calculates all possible set of index permutations given the parameter lists
calculateIndex <- function(params) {

	## Calculate the number of parameters for each variable
	param.max = as.numeric(lapply(params, length))
	
	## Total iterations
	total = prod(param.max)
	num.params = length(params)
	
	## Set up a index matrix
	index <- matrix(0, nrow=total, ncol=num.params)
	
	for (i in 1:num.params) {
		if(i == 1) {
			e = 1	
		} 
		else {
			e = prod(param.max[1:(i-1)])
		}
	
		num = total/prod(param.max[1:i])
		
		index[,i] = rep(rep(1:param.max[i], each=e), num)
	}
	
	return(index)
}	
