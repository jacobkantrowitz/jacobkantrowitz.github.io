# removePercentES takes an expressionSet, a percentage (0.XX) and a function for filtering
# the gene data from the expressionSet. The genes are sorted in ascending order based
# on the method provided and then the bottom XX% are removed. The expressionSet object is
# returned with the genes removed from the expression matrix.
removePercent <- function(exprData, percent, method)
{
	# calculate the function (e.g. mean, coefficient of variance) for each gene
	methodApplyData <- apply(exprs(exprData),1,method)
	
	# sort and index the resulting function (e.g. mean) values
	sortMethodData <- sort(methodApplyData, index.return = TRUE)
	
	# find the index in the sorted indices for the cutoff
	# i.e. if removing 25% of 10,000 genes, then the cutoff would be 2500
	minInd <- percent*length(sortMethodData$ix)
	
	# find the index of the gene that will be last to be included (i.e. last in the sort)
	# i.e. if have 10,000 genes then this will be 10,000
	maxInd <- length(sortMethodData$ix)
	
	# pull the subset of gene indices that will be kept
	# i.e. if keeping 2500/10000 genes these will be sorted.index[2500:10000]
	keepInd <- sortMethodData$ix[minInd:maxInd]
	
	# subset the original expressionSet and include only those genes above the cutoff,
	# removing the desired percentage of the genes based on the sort method provided
	exprs(exprData) <- exprs(exprData)[keepInd,]
	
	#keep_data$ix <- sort_method_data$ix[min_ind:max_ind]
	#keep_data$x <- sort_method_data$x[min_ind:max_ind]
	
	return(exprData)

}