# number_of is simply a shortened call for 'length(which(x))', a frequently used
# set of commands to determine how many items meet some criteria
# number_of returns the number of items in the given object that are TRUE
number_of <- function(bools)
{

	return(length(which(bools)))
}