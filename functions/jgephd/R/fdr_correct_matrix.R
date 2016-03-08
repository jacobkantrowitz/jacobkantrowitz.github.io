# fdr_correct_matrix takes a matrix of p-values and returns a matrix of fdr-corrected values
fdr_correct_matrix <- function(ps)
{
	temp <- ps
	for(column in 1:dim(ps)[2])
	{
		temp[,column] <- p.adjust(temp[, column], method="fdr")
	}
	return(temp)
}
