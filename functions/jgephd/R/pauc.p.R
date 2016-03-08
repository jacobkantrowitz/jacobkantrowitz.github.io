# pauc.p is a function written by Marc Lenburg and downloaded from the pulmonomics wiki
# this function has not been updated to work with expressionSet data
pauc.p <- function(d, class.vector, p = 0.1, n = 10)
{
	# d = data matrix with samples in columns
	# class.vector = class annotation for columns of d (see help for rowpAUCs function in genefilter for more info)
	# p = 1 - specificity value at which to determine pAUC (e.g. p = 0.1 is integral of sensitivity over the range of 90 - 100% specificity)
	# n = times to permute class labels and calculate pAUCs for each gene in dataset (used to estimate null distribution). So if the dataset has 20,000 genes, n = 10 will give you 200,000 pAUCs based on 10 different permutations of the true class label.
	require(genefilter)
	actual <- area(rowpAUCs(d, class.vector, p)) # this is how to calculate pAUC. This is the only useful thing in this function
	permuted <- NULL
	for (i in 1:n)
	{
		# this is a crude attempt to estimate the null distribution of pAUCs (globally rather than per gene)
		permuted <- c(permuted, area(rowpAUCs(d, sample(class.vector, length(class.vector)), p)))
		cat("*")
	}
	cat("\n")

	empiric.p <- function(val, random.values)
	{
		return(sum(random.values > val))
	}
	
	# the next step is not optimized for speed, but determines the empiric p-value for each observed pAUC based on the estimated null distribution
	p.val <- apply(cbind(actual), 1, function(x) empiric.p(x, permuted))
	p.val <- p.val / length(permuted)
	return(p.val)

}