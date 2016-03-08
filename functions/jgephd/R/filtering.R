# Filtering uses removePercent to remove a percentage of a data set's features
# based on a using a mean and/or coefficient of variance method for sorting
# Filtering takes an expressionSet object and optionally the percents to remove
# from the bottom of the data set (e.g. the features with the lowest coefficient of variance)
filtering <- function(dataToFilter, percent1=0.25, percent2=0.25)
{
	cat("Filtering gene expression data...............\n")
	
	dataFiltered <- dataToFilter
	# MULTIPLE FILTERING METHODS

	# 1. BY COEFFICIENT OF VARIANCE
	cv <- removePercent(dataToFilter, percent1, coefficient.of.variance)
	cv.mn <- removePercent(cv, percent2, mean)

	# 2. BY MEAN
  	mn <- removePercent(dataToFilter, percent1, mean)
  	mn.cv <- removePercent(mn, percent2, coefficient.of.variance)
  	
	# 3. BY SUBJECT EXPRESSION LEVELS
  	# STILL TO IMPLEMENT
	
	cat("Completed filtering\n")
	
	# Return unfiltered data along with differently filtered data sets
	return(list("unfiltered"=dataToFilter, "cv"=cv, "cv.mn"=cv.mn, "mn"=mn, "mn.cv"=mn.cv))
}