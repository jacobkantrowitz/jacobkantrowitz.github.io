# setDataFilename sets the filename of the expressionSet RDS file
# this file will be sought in the directory specified by the global variable dataDir
# setDataFilename sets the global variable dataFileName
setDataFileName <- function(fileName)
{
	dataFileName <<- fileName
}