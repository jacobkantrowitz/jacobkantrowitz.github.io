# save.packages saves the installed R packages at the time of running
# this can be used to verify versions of packages used in a certain analysis
# the package information is saved into a file with the date of run
#' Saves the list of currently installed packages to to the current directory to facilitate reproducibility
#' 
#' @examples
#' save.packages()
save.packages <- function()
{
	cat("\nSaving packages\n")
	int <- installed.packages()
	dt <- as.character(Sys.Date())
	write.table(int, file=paste(dt,"_packages.csv", sep=""), sep=",")
	cat("Saved OK\n")
}