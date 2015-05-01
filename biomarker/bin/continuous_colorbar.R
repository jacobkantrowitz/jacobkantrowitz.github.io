continuous.colorbar <- function(var, low.color="blue", mid.color="none", high.color="red", type=c("IQR", "uniform"), iqr.scalar=1.5) {
	
	## This function will generate a continuous color bar based on a continuous variable
	## Written by Josh Campbell, 8-23-09
	## Requires the "gplots" package

	## Parameters #############################################################################
	
	# var:		The continuous variable
	# low.color:	Color for the lowest value.  Default is blue.
	# mid.color:	Not required.  Color for the mid range values.  Defaults to no 
	#			middle color.		
	# high.color:	Color for the highest value.  Default is red.
	# type:		This determines how the color bar is made.  If type is "uniform",
	#			the color bar will be generated over the entire range of the
	#			continous variable (from the minimum value to the maximum).  If 
	#			type is "IQR", the color bar will be generated from the 1st 
	#			Quartile - (IQR * iqr.scalar) to the 3rd Quartile + (IQR * iqr.scalar)
	#			where IQR is the inner-quartile range.  Any value outside of this
	#			range will be given the highest (or lowest) color possible.  
	#			"IQR" is the default because it is more robust against outliers.
	# iqr.scalar:	Only used if "IQR" is chosen for "type".  

	##########################################################################################

		
	## Load the gplots package
	library(gplots, verbose=F)
	
	## Convert variable to numeric datatype
	var <- as.numeric(var)

	## Generates colorbar using the gplots function "colorpanel"
	if(mid.color != "none") {
		col.bar <- colorpanel(length(var), low=low.color, mid=mid.color, high=high.color)
	}
	else {
		col.bar <- colorpanel(length(var), low=low.color, high=high.color)
	}
	
	## Sets the parameter type if it was not explicitly set
	if(length(type > 1)) {
		type <- type[1]
	}
	
	## Obtain the summary of the variable which has max, min, etc.
	s <- summary(var)

	## Set up the color index.  If "uniform" is selected, the color index
	## will span the entire range of the data (from max to min).  If "IQR"
	## is chosen, the color index will span from 
	## 1st Quartile - (IQR * iqr.scalar) to 3rd Quartile + (IQR * iqr.scalar)
	col.index <- c()
	if (type=="uniform") {
		col.index <- seq(from=s[1], to=s[6], length=length(var))
	}
	if (type=="IQR") {
		iqr <- iqr.scalar*(s[5] - s[2])
		low <- s[2] - iqr
		high <- s[5] + iqr
		col.index <- seq(from=low, to=high, length=length(var))
	}

	## This loop will find which color index each value is closest to and
	## assign  to that value the color corresponding to that index
	var.color <- rep(mid.color, length(var))
	for (i in 1:length(var)) {
		var.sub <- abs(var[i] - col.index)
		var.sub.min <- which(var.sub == min(var.sub))
		var.color[i] <- col.bar[var.sub.min]
	}

	return(var.color)
}
