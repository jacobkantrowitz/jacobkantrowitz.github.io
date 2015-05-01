# save_entrez saves a file with the entrez IDs of the genes of interest
# save_entrez takes a vector of indices, a vector of rownames, and a filename
# the indices must refer to the rownames and will be used to pull the gene entrez IDs
# from the rownames; these IDs will be saved in the file with the given filename
save_entrez <- function(indices, rownms, filename)
{
	entrez_codes <- return_entrez(rownms[indices])
	lapply(entrez_codes, write, filename, append=TRUE)
	cat("Saved ", length(indices), " gene entrez codes to\n", getwd(), "/",filename, "\n\n", sep="")
}

save_geneSymbols <- function(model, filename){
  lapply(model$geneSymbols, write, filename, append=TRUE)
  cat("Saved ", length(model$geneSymbols), " gene Symbols to\n", getwd(), "/",filename, "\n\n", sep="")
  
  
}
