# Author: Joe Perez-Rogers
# Date: 2014-05-07
# Script to reload all packages and functions from the current project bin
# Usage: reloadFunctions()

reloadFunctions <- function(){
	source("../../bin/sourceDirectory.R")
	sourceDir("../../bin/")
}
