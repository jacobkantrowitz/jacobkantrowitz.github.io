# print_status is a visualization tool for representing a progress bar while running large
# analyses; a * is printed for each 2% of an analysis completed. 
# print_status takes two inputs, 'current', and 'total', which represent the current
# level of completion compared to the total level
# (e.g. currently running gene 1098 of 19684 total genes)
print_status <- function(current, total)
{
	pct <- round(current/total * 100)
	prior_pct <- round((current-1)/total*100)
	if((pct %% 2)==0 & (prior_pct %% 2)==1)
	{
		cat("*", sep="")
	}
	
}