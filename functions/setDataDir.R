# setDataDir sets the directory where the expressionSet RDS exists
# setDataDir takes the directory where the data of interest resides
# setDataDir sets the global variable 'dataDir'
setDataDir <- function(pathName)
{
	dataDir <<- pathName
}