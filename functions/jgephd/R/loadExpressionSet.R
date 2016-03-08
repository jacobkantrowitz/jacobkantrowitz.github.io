# loadExpressionSet load the expressionSet found in the global or given RDS file
# loadExpressionSet takes pathName and dataFileName but uses the global variables if none given
# loadExpressionSet checks that the loaded object is an expressionSet before returning it
loadExpressionSet <- function(pathName=dataDir, fileName=dataFileName)
{
	cat("Loading data............... ", pathName, fileName, "\n", sep="")
	tmpName <- paste(pathName, fileName, sep="")
	tmpObj <- readRDS(tmpName)
	
	cat("Data loaded successfully\nChecking object type\n")
	if(class(tmpObj)=="ExpressionSet")
	{
		cat("Object is of class ExpressionSet\n")
		return(tmpObj)
	}
	else
	{
		cat("ERROR: Object is NOT of class ExpressionSet\n\nSTOP\n")
	}
}
