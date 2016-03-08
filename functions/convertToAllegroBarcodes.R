# Author: Joe Perez-Rogers
# Date: 2014-05-07
# Script to convert my barcodes to Allegro's barcodes (e.g. convert 1-1-0001 to 10001)
# Usage: convertToAllegroBarcodes(my.barcodes)
# Inputs: vector of barcodes (factors or characters)
# Outputs: vector of updated barcodes (characters)

convertToAllegroBarcodes <- function(x){
	as.character(as.numeric(unlist(lapply(x,
		function(i){
		paste(unlist(strsplit(as.character(i),"-"))[2],unlist(strsplit(as.character(i),"-"))[3],sep="")
		}
	))))
}


