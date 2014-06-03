# return_lm calculates and returns linear models based on the genes from an expressionSet
# and the covariates in the separate 'covariates' list
# return_lm takes two items, an expressionSet object and a list of covariates
return_lm <- function(exprData, covariates)
{
	number_of_genes <- featureNumber(exprData)
	number_of_covs <- length(covariates)
	
	# CREATE SPACE FOR MODELS 
	models <- vector("list", number_of_genes)

	# CREATE SPACE FOR EXTRACTED P-VALUES FROM GENE MODELS
	number_of_levels <- calculate_total_levels(exprData, covariates)
	model_ps <- matrix(nrow=number_of_genes, ncol=number_of_levels)
	
	# CREATE SPACE FOR MODEL RESIDUALS
	model_resids <- exprs(exprData)
	model_resids[,] <- 0

	# RUN MODELS
	model <- paste(covariates, collapse=" + ")
	cat("Running one model:", model, "\nRunning model for", number_of_genes, "genes\n")

	# Print status bar
	cat(rep("*", 50),"\n", sep="")
	
	for(gene in 1:number_of_genes)
	{	
		gene_name <- featureNames(exprData)[gene]
		#cat("Running models for gene #", gene, " name:", gene_name, "\n")
		print_status(gene, number_of_genes)
		
		models[[gene]] <- lm(as.formula(paste("exprs(exprData)[", gene, ",] ~", model, sep="")), data=pData(exprData))
		models[[gene]]$gene_name <- gene_name
				
		number_of_coeffs <- nrow(summary(models[[gene]])$coefficients)
		
		#cat("number of coeffs", number_of_coeffs, "\n", "dim model_ps", dim(model_ps), "\n")
		#print(summary(models[[gene]])$coefficients[2:number_of_coeffs,4])
		#print(summary(models[[gene]])$coefficients)
		
		#often problematic line due to NAs present in one of the pData fields
		#sometimes also problem due to lack of variability of one phenotype type with the others
		model_ps[gene,] <- summary(models[[gene]])$coefficients[2:number_of_coeffs,4]

		model_resids[gene,] <- summary(models[[gene]])$residuals
		
	}
	# Print a new line after finishing the progress bar
	cat("\n")
	
	covariate_names <- rownames(summary(models[[gene]])$coefficients)[2:number_of_coeffs]
	model_fdr <- fdr_correct_matrix(model_ps)
	colnames(model_fdr) <- covariate_names
	colnames(model_ps) <- covariate_names
	
	model_formula <- paste("gene ~ ", paste(covariates, collapse=" + "))
	
	cutoffs_p <- return_cutoff_lists(model_ps, UNCORRECTED)
	cutoffs_fdr <- return_cutoff_lists(model_fdr, CORRECTED)
	
	results <- list("models"=models,
					"model_ps"=model_ps,
					"cutoffs_p"=cutoffs_p,
					"model_fdr"=model_fdr,
					"cutoffs_fdr"=cutoffs_fdr,
					"model_resids"=model_resids,
					#"data"=data_temp,
					#"phenotype"=phenotype,
					"covariates"=covariates,
					"formula"=model_formula)
					
	return(results)

}