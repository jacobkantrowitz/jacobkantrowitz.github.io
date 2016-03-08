# return_cutoff_lists generates lists of indices for which the p-value or fdr-value
# given is less than the prescribed cutoff (e.g. fdr < 0.25, fdr < 0.2, fdr < 0.1, etc.)
# return_cutoff_lists takes a matrix of values (p or fdr) and returns lists of indices
# into that matrix (row only, i.e. indicating the gene) that fit the criteria of having a
# p or fdr value below the given cutoff
return_cutoff_lists <- function(values, corrected)
{
	temp <- list()
	#cutoffs <- c(0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
	
	if(corrected)
	{
		cutoffs <- f_cutoffs
	}
	else
	{
		cutoffs <- p_cutoffs
	}
	
	for(value in 1:length(cutoffs))
	{
		temp[[value]] <- list()
		names(temp)[[value]] <- toString(cutoffs[[value]])
		for(column in 1:ncol(values))
		{
			temp[[value]][[column]] <- which(values[,column] < cutoffs[value])
		}
		names(temp[[value]]) <- colnames(values)
	}
	return(temp)
}