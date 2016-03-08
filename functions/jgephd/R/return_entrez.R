# return_entrez takes a list of Affy rownames (i.e. those ending with "_at" and 
# removes the "_at", returning a list of entrez gene IDs
return_entrez <- function(rownms)
{
	temp <- sub("_at", "", rownms)
	return(temp)
}