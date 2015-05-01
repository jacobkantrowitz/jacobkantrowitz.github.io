# Author: Joe Perez-Rogers (adapted from ?source examples in R documentation
# Date: 2014-05-07
# Function to source all R scripts in a directory
# Usage: sourceDir("/path/to/my/scripts.R")

sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(trace) cat(nm,":")           
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}

