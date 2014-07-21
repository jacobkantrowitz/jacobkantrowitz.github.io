# Script for analysis of Allegro data (708); COPD vs CANCER
# Jake Kantrowitz, BUSM MD/PhD Candidate 2019, Molecular and Translational Med.
# jacob.kantrowitz@gmail.com, kantro@bu.edu
# (617)-584-2393
# 9/30/13
# We are looking for an interaction between cancer and COPD in the Allegro dataset
# This is a continuation of the work I did over the summer of 2013 as a rotation student

# THINGS TO DO 
# Perform analysis of copd/cancer vs 3 other groups
# Perform analysis of copd/cancer vs all else
# Perform analysis of copd/cancer vs cancer
# Filter with different methods (mean, coef_var, none, etc.)

# Use SVA to find structure in data, correlate with possible covariates and then
# 	include those covariates that have high correlation with a surrogate variable
# 	in subsequent analyses.

# Alternative method for including covariates: just use AIC plot to find model with 
# 	best fit... seems more honest and sound. This isn't really possible though because
#	the data is a matrix and not just one variable. There isn't a single 'response'
#	variable we're trying to explain. This is why SVA is actually good.


# Define COPD in new ways and run analysis with each definition/term
#	1. Allegro CRF definition (what I used initially)
#	2. FEV1/FVC ratio as continuous metric instead of binary operator (40% have PFTs)
#	3. Binary COPD by PFTs (Lower Limit of Normal or GOLD)
#	4. FEV1 a continuous ratio
#	5. (FEV1/FVC) / LLN --> offers some sort of normalization for better comparison between patients
# What else?

# Want genes

# FOR EACH ANALYSIS, SAVE IN A SEPARATE DIRECTORY:
#	- separate script - is this necessary? maybe better to have modifiable script with good annotations/instructions
#	- list of genes of interest with p-values, fdr-values, higher/lower indicator
#	- heatmap of genes of interest
#	- readme file with filtering method, number genes analyzed, what else?



#########################################################
############## LOAD NECESSARY LIBRARIES #################
#########################################################

# load bioinformatics packages
library(affy)

# set package library to load heatmap3
.libPaths(c(sprintf("/unprotected/projects/cbmhive/R_packages/R-%s", getRversion()), .libPaths()))
library(heatmap3)
library(sva)

#########################################################
############ DEFINE GLOBAL VARIABLES  ###################
#########################################################

# SAVE INITIAL DIRECTORY
initial_directory <- getwd()

# DEFINE DATA DIRECTORY AND FILENAME
data_directory <- "/protected/projects/pulmarray/Biollegro/RDS"
data_filename <- "Bronc_708_First_CEL.rds"

# DEFINE THE NUMBER OF GENES INCLUDED IN CURRENT ANALYSIS
number_of_genes <- 0

# DEFINE CUTOFFS TO USE FOR FDR AND P-VALUE SIGNIFICANCE
f_cutoffs <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
p_cutoffs <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)

CORRECTED <- 1
UNCORRECTED <- 0


########## Color Schemes ########## 

bluered <- colorRampPalette(c("blue","white","red"))(256)

# Define color schemes, based on factors
cancer_colors <- c("0"="black", "1"="red");
copd_colors <- c("0"="bisque", "1"="blue");
gender_colors <- c("0"="blue", "1"="pink");
smoking_colors <- c("1"="grey", "2"="white");
batch_colors <- c("12"="red", "13"="blue", "14"="orange", "15"="green", "16"="black");
copdca_colors <- c("0"="green", "1"="yellow", "2"="orange", "3"="red")


#########################################################
############ DEFINE HELPER FUNCTIONS  ###################
#########################################################

# remove_percent takes a set of data, a percentage (0.XX) and a function for filtering
# the dataset; the data set is sorted in ascending order based on the provides filtering
# function and then the bottom XX% are removed; the expression dataset is returned 
# with the XX% removed
remove_percent <- function(expr, percent, method)
{

	method_data <- apply(expr,1,method)
	sort_method_data <- sort(method_data, index.return = TRUE)
	min_ind <- percent*length(sort_method_data$ix)
	max_ind <- length(sort_method_data$ix)
	
	keep_ind <- sort_method_data$ix[min_ind:max_ind]
	keep_data <- expr[keep_ind,]
	
	#keep_data$ix <- sort_method_data$ix[min_ind:max_ind]
	#keep_data$x <- sort_method_data$x[min_ind:max_ind]
	
	return(keep_data)

}


coef_of_var <- function(x)
{
	# coef_of_var is a simple function that returns the coefficient of variance
	# for a variable X (e.g. a data frame)

	sd(x) / abs(mean(x))
}


save_analysis <- function(models)
{
	# FOR EACH ANALYSIS, SAVE IN A SEPARATE DIRECTORY:
	#	- separate script
	#	- list of genes of interest with p-values, fdr-values, higher/lower indicator
	#	- heatmap of genes of interest
	#	- readme file with filtering method, number genes analyzed, what else?


}


init <- function()
{
	# SET THE WORKING DIRECTORY TO THE DATA DIRECTORY
	setwd(data_directory)
	
	# PRINT CURRENT OPERATIONS
	cat("Loading data...", data_directory, data_filename, "\n")
	
	# LOAD DATA
	data_loaded <- readRDS(data_filename)
	
	# EXTRACT THE PHENOTYPES
	data_phenotype <- pData(data_loaded)
	
	# EXTRACT THE EXPRESSION LEVELS
	data_expression <- exprs(data_loaded)
	
	# DEFINE THE NUMBER OF GENES IN THE EXPRESSION DATA
	assign("number_of_genes", dim(data_expression)[1], envir=.GlobalEnv)

	# SET THE OBJECTS TO RETURN AS ATTRIBUTES OF THE data_loaded OBJECT
	attributes(data_loaded) <- list("expression"=data_expression, "phenotype"=data_phenotype)

	return(data_loaded)
}


prepare_data <- function(data_to_analyze)
{
	# COERCE DATA TO APPROPRIATE TYPES
	
	# Column 1 - BATCH
	data_to_analyze$BATCH <- factor(data_to_analyze$BATCH)
	
	# Column 3 - CANCER
	data_to_analyze$CANCER <- factor(data_to_analyze$CANCER)
	
	# Column 4 - SMK
	data_to_analyze$SMK <- factor(data_to_analyze$SMK)
	
	# Column 5 - COPD
	data_to_analyze$COPD <- factor(data_to_analyze$COPD)
	
	# Column 6 - GENDER
	data_to_analyze$GENDER <- factor(data_to_analyze$GENDER)

	# Create a new indicator function such that I = 1 when COPD=1 and CA=1, 0 otherwise
	data_to_analyze$COPDCA <- factor(data_to_analyze$COPD == 1 & data_to_analyze$CANCER == 1)
	
	# Create a new indicator function that I = 1 when COPD=0 and CA=0, 0 otherwise
	data_to_analyze$HEALTHY <- factor(data_to_analyze$COPD == 0 & data_to_analyze$CANCER == 0)
	
	data_to_analyze$indicator <- 0
	data_to_analyze$indicator[which(data_to_analyze$CANCER==1)] <- 1
	data_to_analyze$indicator[which(data_to_analyze$COPD==1)] <- 2
	data_to_analyze$indicator[which(data_to_analyze$COPDCA==TRUE)] <- 3
	data_to_analyze$indicator <- factor(data_to_analyze$indicator)
	

	return(data_to_analyze)
}


filtering <- function(data_to_filter)
{
	cat("Filtering data...\n")
	
	data_filtered <- data_to_filter
	# MULTIPLE FILTERING METHODS

	# 1. BY COEFFICIENT OF VARIANCE
	cv <- remove_percent(data_to_filter, .25, coef_of_var)
	cv_mn <- remove_percent(cv, .25, mean)

	# 2. BY MEAN
  	mn <- remove_percent(data_to_filter, .25, mean)
  	mn_cv <- remove_percent(mn, .25, coef_of_var)
  	
	# 3. BY SUBJECT EXPRESSION LEVELS
  	# STILL TO IMPLEMENT

	# SET THE FILTERED SETS AS ATTRIBUTES OF UNFILTERED SET
	attributes(data_filtered) <- list("cv"=cv, "cv_mn"=cv_mn, "mn"=mn, "mn_cv"=mn_cv)
	
	# RETURN UNFILTERED DATA WITH DIFFERENTLY FILTERED SETS AS ATTRIBUTES
	return(data_filtered)
}


return_lm <- function(data_to_analyze, phenotype, covariates)
{
	number_of_genes <- dim(data_to_analyze)[1]
	number_of_covs <- length(covariates)
	
	# CREATE SPACE FOR MODELS 
	models <- vector("list", number_of_genes)

	# CREATE SPACE FOR EXTRACTED P-VALUES FROM GENE MODELS
	number_of_levels <- calculate_total_levels(phenotype, covariates)
	model_ps <- matrix(nrow=number_of_genes, ncol=number_of_levels)
	
	# CREATE SPACE FOR MODEL RESIDUALS
	model_resids <- data_to_analyze
	model_resids[,] <- 0

	# RUN MODELS
	model <- paste(covariates, collapse=" + ")
	cat("Running one model:", model, "\nRunning model for", number_of_genes, "genes\n")

	# Print status bar
	cat(rep("*", 50),"\n", sep="")
	
	for(gene in 1:number_of_genes)
	{	
		gene_name <- rownames(data_to_analyze)[gene]
		#cat("Running models for gene #", gene, " name:", gene_name, "\n")
		print_status(gene, number_of_genes)
		
		models[[gene]] <- lm(paste("data_to_analyze[", gene, ",] ~", model, sep=""), data=phenotype)
		models[[gene]]$gene_name <- gene_name
				
		number_of_coeffs <- nrow(summary(models[[gene]])$coefficients)
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
	
	results <- list("models"=models,
					"model_ps"=model_ps,
					"model_fdr"=model_fdr,
					"model_resids"=model_resids,
					"data"=data_to_analyze,
					"phenotype"=phenotype,
					"covariates"=covariates,
					"formula"=model_formula)
					
	return(results)

}

#results <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs)


# Check that covariates exist within the given phenotype data
check_covariates <- function(phenotype, covariates)
{
	FAIL = 1
	PASS = 0
	failed <- numeric()
	
	# if the supplied 'covariates' list represents one model
	if(is.character(covariates[[1]]))
	{
		isolate_covariates <- strsplit(as.character(covariates), "*", fixed=TRUE)
		
		all_covariates <- character()
		for(covs in isolate_covariates)
		{
			for(cov in covs)
				all_covariates <- append(all_covariates, cov)
		}
		
		all_covariates <- unique(all_covariates)
		cat("Unique covariates:", all_covariates, "\n")
		
		# each char item in 'covariates' must be from the names of 'phenotype'
		covs_available <- number_of(is.na(match(all_covariates, names(phenotype))))
		# if at least one supplied covariate is not in the phenotype data
		if(covs_available)
		{
			cat("Some covariates are not traits or are spelled incorrectly:\n",
					as.character(covariates[is.na(match(all_covariates, names(phenotype)))]),
					"\n")
			model <- paste(covariates, collapse=" + ")
			cat("Model to edit:", model, "\n\n")
			return(FAIL)
				
		}
		# else all the supplied covariates are in the phenotype data
		else
		{
			model <- paste(covariates, collapse=" + ")
			cat("Model OK:", model, "\n\n")
			return(PASS)
		}
	}
	
	# else if the supplied 'covariates' list represents multiple models
	else if(is.list(covariates[[1]]))
	{
		cat("\nMultiple models to check:", length(covariates), "models\nChecking Models\n\n")
		for(covs in 1:length(covariates))
		{
			cat("Checking model #", covs, ":", paste(covariates[[covs]], collapse=" + "), "\n", sep="")
			test <- check_covariates(phenotype, covariates[[covs]])
			
			if(test)
			{
				failed <- append(failed, covs)
			}
			#cat("Done check model #", covs, "\n\n")
		}
		
		if(length(failed)==0)
		{
			cat("All models passed: OK\n")
			return(PASS)
		}
		else
		{
			cat("\nModels failed: ", length(failed), "/", length(covariates), "\n\n", sep="")
			cat("Edit models: ", paste(failed, collapse=", "), "\n\n", sep="")
			return(FAIL)
		}
	}
	
	# else the first item of covariates is neither a char nor a list i.e. wrong data type
	else
	{
		cat("\nERROR: The covariates data supplied is not of the right class.\n")
		cat("Covariates of class", class(covariates), "\n\n")
		return(FAIL)
	}

}


# analyze custom lists of covariates
analyze_covariate_lists <- function(data_to_analyze, phenotype, covariates)
{
	results <- 0
	
	test <- check_covariates(phenotype, covariates)
	
	# if all covariates check out OK --> proceed with analysis
	if(test==0)
	{
		# if 'covariates' is a single list of covariates, run lm on data_to_analyze
		# using the covariates from the list
		if(is.character(covariates[[1]]))
		{
			model <- paste(covariates, collapse=" + ")
			#cat("Running one model:", model, "\n")
			results <- return_lm(data_to_analyze, phenotype, covariates)
			return(results)
		}
		
		# else the first item of covariates is another list i.e. covariates is a list of lists
		else if(is.list(covariates[[1]]))
		{
			results <- list()
			cat("Multiple models to run:", length(covariates), "models\nRunning Models\n\n")
			for(model in 1:length(covariates))
			{
				cat("Running model #", model, "\n", sep="")
				temp <- return_lm(data_to_analyze, phenotype, covariates[[model]])
				results[[paste(covariates[[model]], collapse="+")]] <- temp
				cat("Done running model #", model, "\n\n", sep="")
			}
		}
	}
	
	# else some covariates are not available --> cannot create all models
	else
	{
		cat("\nSome of the covariates provided were not found in the phenotype data.\n",
		"Please edit the indicated covariates and resubmit.\n\n")
	
	}
	
	return(results)
}

#analyze_covariate_lists(data_to_analyze, data_phenotype, covs2)


fdr_correct_matrix <- function(ps)
{
	temp <- ps
	for(column in 1:dim(ps)[2])
	{
		temp[,column] <- p.adjust(temp[, column], method="fdr")
	}
	return(temp)
}


#print_status <- function(current, total)
#{
#	pct <- round(current/total * 100)
#	cat(pct,"% of models run\n", sep="")
#}

print_status <- function(current, total)
{
	pct <- round(current/total * 100)
	prior_pct <- round((current-1)/total*100)
	if((pct %% 2)==0 & (prior_pct %% 2)==1)
	{
		cat("*", sep="")
	}
	
}


return_sva <- function(data_expression, data_phenotype)
{

	# There can be no NA values in the phenotype data used to create the models for SVA

	# Build models for inclusion in SVA
	mod1 <- model.matrix(~COPD + CANCER + COPD*CANCER, data=data_phenotype)
	mod0 <- model.matrix(~1, data=data_phenotype)
	
	# Construct the surrogate variables and return them
	
	sva_obj <- sva(data_expression, mod1, mod0)
	
	return(sva_obj)
}

number_of <- function(bools)
{

	return(length(which(bools)))
}


return_entrez <- function(rownms)
{
	temp <- sub("_at", "", rownms)
	return(temp)
}

save_entrez <- function(indices, rownms, filename)
{
	entrez_codes <- return_entrez(rownms[indices])
	lapply(entrez_codes, write, filename, append=TRUE)
	cat("Saved ", length(indices), " gene entrez codes to\n", getwd(), "/",filename, "\n\n", sep="")
}

return_cutoff_lists <- function(values, corrected)
{
	temp <- list()
	#cutoffs <- c(0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
	
	if(corrected)
	{
		cutoffs <- f_cutoffs
	}
	else
	{
		cutoffs <- p_cutoffs
	}
	
	for(value in 1:length(cutoffs))
	{
		temp[[value]] <- list()
		names(temp)[[value]] <- toString(cutoffs[[value]])
		for(column in 1:ncol(values))
		{
			temp[[value]][[column]] <- which(values[,column] < cutoffs[value])
		}
		names(temp[[value]]) <- colnames(values)
	}
	return(temp)
}

fp_cutoffs <- function(model_object)
{
	f_cutoffs <- return_cutoff_lists(model_object$model_fdr, CORRECTED)
	p_cutoffs <- return_cutoff_lists(model_object$model_ps, UNCORRECTED)
	
	return(list("fdr_cutoffs"=f_cutoffs, "p_cutoffs"=p_cutoffs))	
}


calculate_total_levels <- function(phenotype, covariates)
{
	total <- 0
	for(cov in covariates)
	{
		if(grepl("*", cov, fixed=TRUE))
		{
			total <- total + 1
		}
		else
		{
			col <- match(cov, names(phenotype))
			if(is.factor(phenotype[,col]))
			{
				temp <- length(levels(phenotype[,col])) - 1
				#print(temp)
				total <- total + temp
			}
			else
			{
				total <- total + 1
			}
		}
	}
	return(total)
}


remove_genes <- function(inds, data)
{
	cat("Removing ", length(inds), " genes from the data.\n",
		(nrow(data)-length(inds)), "genes remain.\n\n")
		
	remain_ind <- which(is.na(match(1:nrow(data), inds)))
	return(data[remain_ind,])
}

keep_genes <- function(inds, data)
{
	cat("Keeping ", length(inds), " genes from the data.\n",
		(nrow(data)-length(inds)), "genes removed.\n\n")
		
	remain_ind <- which(!is.na(match(1:nrow(data), inds)))
	return(data[remain_ind,])
}

remove_patients <- function(inds, data, phenotype)
{
	cat("Removing ", length(inds), " patients from the data.\n",
		(ncol(data)-length(inds)), "patients remain.\n\n")
		
	remain_ind <- which(is.na(match(1:ncol(data), inds)))
	temp_data <- data[,remain_ind]
	temp_phenotype <- phenotype[remain_ind,]
	return(list(temp_data,temp_phenotype))
}

max_genes_heatmap <- 200
min_genes_heatmap <- 6


## HEATMAP FUNCTIONS
generate_heatmaps <- function(model_object)
{
	clabels <- cbind(
		"Indicator" = copdca_colors[model_object$phenotype$indicator],
		"COPD Status" = copd_colors[model_object$phenotype$COPD],
		"Smoking Status" = smoking_colors[model_object$phenotype$SMK],
		"Cancer Status" = cancer_colors[model_object$phenotype$CANCER],
		"Gender" = gender_colors[model_object$phenotype$GENDER],
		"Batch" = batch_colors[model_object$phenotype$BATCH]
	)

	cutoffs <- fp_cutoffs(model_object)
	
	type_names <- strsplit(names(cutoffs), "_")

	#cat("type names:", paste(type_names, collapse=" "), "\n")
	
	for(type in 1:length(cutoffs)) # go through fdr and p
	{
		value_names <- names(cutoffs[[type]]) # names of the cutoff values
		#cat("value names:", paste(value_names, collapse=" "), "\n")
		
		for(level in 1:length(cutoffs[[type]])) # go through all cutoff levels
		{
			covar_names <- names(cutoffs[[type]][[level]])
			#cat("covar names:", paste(covar_names, collapse=" "), "\n")
			
			for(covar in 1:length(cutoffs[[type]][[level]])) # go through all covariates
			{
				gene_index <- cutoffs[[type]][[level]][[covar]]
				number_sig_genes <- length(gene_index)
				
				if((number_sig_genes <= max_genes_heatmap) & (number_sig_genes >= min_genes_heatmap))
				{
					map_title <- paste(model_object$formula, "\n", number_sig_genes, "genes significant for ", covar_names[[covar]], "\n",
						type_names[[type]][[1]], "<", value_names[[level]], "\n\n", collapse="")
					filename <- paste(initial_directory, "/", number_sig_genes, "Genes", type_names[[type]][[1]], value_names[[level]], covar_names[[covar]], ".pdf", sep="")
					
					cat(map_title)
					
					pdf(filename)
					# generate a heatmap
					#heatmap3(model_object$data[gene_index,], col = bluered, hclustfun=function(d) hclust(d, method="average"),row.dendrogram = FALSE, col.clustering = "semisupervised", ColSideColors = clabels, main = map_title)
					heatmap3(model_object$data[gene_index,], col = bluered, col.clustering = "semisupervised", ColSideColors = clabels, main = map_title)
					
					
					dev.off()
				}
				else
				{
					# do nothing; too many genes
				}
			}
		}
	}
}


#generate_heatmaps(results_resid)

#"Heatmap for X genes FDR 0.05 Model: gene ~ X + Y + Z"




cat("Start analysis\n\n")

#########################################################
################# START ANALYSIS ########################
#########################################################

# Define lists of covariates to use for different model analyses
covars <- list(list("BATCH", "RIN", "SMK", "GENDER", "AGE"),
				list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER"),
				list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER"))

covs_resid2 <- list("BATCH", "RIN", "GENDER", "AGE")
covs_copdca_intrx <- list("COPD", "CANCER", "COPD*CANCER")
covs_pulm <- list("SMK", "COPD", "CANCER", "COPD*CANCER")
covs_interesting <- list("SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER")

# Load and prep the data for analysis
data_new <- init()

# Extract the expression and phenotype data from the loaded data
data_phenotype <- attributes(data_new)$phenotype
data_expression <- attributes(data_new)$expression

# Remove patients with phenotype data missing (i.e. NAs)
# THIS SHOULD BE MOVED INTO THE init() FUNCTION
remove_pts <- numeric()
for(trait in 1:length(data_phenotype))
{
	remove_pts <- append(remove_pts, which(is.na(data_phenotype[,trait])))

}
# also remove patients who are never smokers i.e. SMK==3
remove_pts <- append(remove_pts, which(data_phenotype$SMK==3))

# remove the duplicate indices
remove_pts <- unique(remove_pts)

# create a true/false vector with false indicating a patient index to remove
remove_pts <- is.na(match(1:dim(data_phenotype)[1], remove_pts))
# remove the appropriate patient data from the phenotype data
data_phenotype <- data_phenotype[remove_pts,]
# remove the appropriate patient data from the expression data
data_expression <- data_expression[, remove_pts]


# Prepare the phenotype data for use in models i.e. make columns into factors
data_phenotype <- prepare_data(data_phenotype)

# Filter the data and save different filtered sets as attributes
# Filtered sets:
#	a. filtered_cv
#	b. filtered_cv_mn
#	c. filtered_mn
#	d. filtered_mn_cv
data_filtered <- attributes(filtering(data_expression)) # can adjust this function so it just returns a list...

# define number of genes remaining in once or twice filtered data sets
number_of_genes_once <- dim(data_filtered$cv)[1]
number_of_genes_twice <- dim(data_filtered$cv_mn)[1]


# RUN LINEAR MODELS

# set the data (filtered or unfiltered from above) to use in the models
#	1. data_filtered$cv		filtered by coefficient of variance
#	2. data_filtered$cv_mn	filtered by cv and then by mean
#	3. data_filtered$mn		filtered by mean
#	4. data_filtered$mn_cv	filtered by mn and then by cv
#	5. data_expression		unfiltered data (all genes included)
#	6. other filtering methods
data_to_analyze <- data_expression

# Analysis
# Remove effects due to BATCH, RIN, SMK, GENDER, AGE
# Model: gene ~ BATCH + RIN + SMK + GENDER + AGE
covs_resid <- list("BATCH", "RIN", "SMK", "GENDER", "AGE")  
results_resid <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs_resid)

# p-value and FDR cutoffs
cutoffs_resid_p <- return_cutoff_lists(results_resid$model_ps, UNCORRECTED)
cutoffs_resid_fdr <- return_cutoff_lists(results_resid$model_fdr, CORRECTED)


# Analysis
# Model: BATCH + RIN + SMK + GENDER + AGE + SMK + COPD + CANCER
covs_total <- list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER")
results_total <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs_total)

# p-value and FDR cutoffs
cutoffs_total_p <- return_cutoff_lists(results_total$model_ps, UNCORRECTED)
cutoffs_total_fdr <- return_cutoff_lists(results_total$model_fdr, CORRECTED)


# Analysis
# Check for an interaction term in the total model with all available phenotypes
# Model: BATCH + RIN + SMK + GENDER + AGE + SMK + COPD + CANCER + CANCER*COPD
covs_intrx <- list("BATCH", "RIN", "SMK", "GENDER", "AGE", "COPD", "CANCER", "COPD*CANCER")
results_intrx <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs_intrx)

# p-value and FDR cutoffs
cutoffs_intrx_p <- return_cutoff_lists(results_intrx$model_ps, UNCORRECTED)
cutoffs_intrx_fdr <- return_cutoff_lists(results_intrx$model_fdr, CORRECTED)


# Analysis - invalid: can't remove effects of model terms without also affecting dfs and ps
covs_remove <- list("BATCH", "RIN")
results_remove <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs_remove)

results_remove_fdr_cutoffs <- return_cutoff_lists(results_intrx$model_fdr, CORRECTED)



remove_inds <- union(results_total_p_cutoffs[[1]]$BATCH13,
						union(results_total_p_cutoffs[[1]]$BATCH14,
						union(results_total_p_cutoffs[[1]]$BATCH15,
						union(results_total_p_cutoffs[[1]]$BATCH16,
						results_total_p_cutoffs[[1]]$RIN))))
remove_inds <- unique(remove_inds)
						
keep_inds <- union(results_total_p_cutoffs[[1]]$AGE,
						union(results_total_p_cutoffs[[1]]$GENDER1,
						union(results_total_p_cutoffs[[1]]$SMK2,
						union(results_total_p_cutoffs[[1]]$AGE,
						union(results_total_p_cutoffs[[1]]$COPD1,
						results_total_p_cutoffs[[1]]$CANCER1)))))
keep_inds <- unique(keep_inds)

# in remove_inds only keep the indices that are not also in keep_inds
remove_inds <- remove_inds[is.na(match(remove_inds, keep_inds))]

cleaned_data <- remove_genes(remove_inds, data_to_analyze)
results_interesting <- analyze_covariate_lists(cleaned_data, data_phenotype, covs_interesting)
results_interesting_fdr_cutoffs <- return_cutoff_lists(results_interesting$model_fdr, CORRECTED)
results_interesting_p_cutoffs <- return_cutoff_lists(results_interesting$model_ps, UNCORRECTED)

results_key <- analyze_covariate_lists(results_resid$model_resids, results_resid$phenotype, covs_copdca_intrx)


results_resid2 <- analyze_covariate_lists(data_to_analyze, data_phenotype, covs_resid2)
results_pulm <- analyze_covariate_lists(results_resid2$model_resids, results_resid2$phenotype, covs_pulm)






#results_covars <- analyze_covariate_lists(data_to_analyze, data_phenotype, covars)

#results_total_p_cutoffs <- return_cutoff_lists(results_covars[[2]]$model_ps, UNCORRECTED)

# run model: CANCER + COPD + CANCER*COPD on residuals
# residuals from model: BATCH + RIN + SMK + GENDER + AGE
# Should remove those genes that were significant for the terms in the original model
#resids <- results_resid$model_resids
#RIN_resids <- remove_genes(results_resid_fdr_cutoffs[[11]]$RIN, resids)
#SMK_resids <- keep_genes(results_resid_fdr_cutoffs[[4]]$SMK, resids)

#results_smk <- analyze_covariate_lists(SMK_resids, data_phenotype, covs_copdca_intrx)
#results_rin <- analyze_covariate_lists(RIN_resids, data_phenotype, covs_copdca_intrx)

#results_smk_fdr_cutoffs <- return_cutoff_lists(results_smk$model_fdr)
#results_rin_fdr_cutoffs <- return_cutoff_lists(results_rin$model_fdr)

#results_rin_ps_cutoffs <- return_cutoff_lists(results_rin$model_ps, UNCORRECTED)


# TRY REMOVING SUBJECTS WITH RIN < 3 or 4
# THERE ARE 13 PATIENTS OF 700 WITH RIN <= 3
#RIN_cutoff <- 3
#inds <- which(data_phenotype$RIN <= RIN_cutoff)

#temp <- remove_patients(inds, data_to_analyze,data_phenotype)
#data_to_analyze <- temp[[1]]
#data_phenotype <- temp[[2]]


#match(rownames(SMK_resids[results_smk_fdr_cutoffs[[1]]$'COPD1:CANCER1',]),rownames(RIN_resids[results_rin_fdr_cutoffs[[1]]$'COPD1:CANCER1',]))
#rownames(SMK_resids[results_smk_fdr_cutoffs[[1]]$'COPD1:CANCER1',])



#resids <- run_for_residuals(data_to_analyze, data_phenotype)
#results <- analyze2(resids, data_phenotype)


#copdca_intrx_model <- analyze_covariate_lists(copdca_resid_2_analyze, data_phenotype, covs_copdca_intrx)





# Repeat analysis with healthy patients removed
#healthy_ind <- which(data_phenotype$HEALTHY==TRUE)
#diseased_ind <- union(which(data_phenotype$COPD==1), which(data_phenotype$CANCER==1))
#data_diseased <- remove_patients(healthy_ind, data_expression, data_phenotype)
#data_diseased_expr  <- data_diseased[[1]]
#data_diseased_phen <- data_diseased[[2]]
#results_covars_diseased <- analyze_covariate_lists(data_diseased_expr, data_diseased_phen, covs_resid)


#data_diseased_expr2 <- data_expression[,is.na(match(1:700,healthy_ind))]
#data_diseased_phen2 <- data_phenotype[is.na(match(1:700, healthy_ind)),]
#results_covars_diseased <- analyze_covariate_lists(data_diseased_expr2, data_diseased_phen2, covs_resid)


#blech <- analyze_covariate_lists(data_expression[,diseased_ind], data_phenotype[diseased_ind,], covs_resid)


# find genes significant for only RIN or BATCH and not GENDER, AGE, SMK, COPD, or CANCER



# Analysis
# Total model 
# Find genes fdr corrected < 0.05 significant for COPD or CANCER in Total model
# Re-analyze those genes for only pulm covariates (smk, cancer, copd, intrx)
copd_genes <- which(results_total$model_fdr[,9] < 0.05)
copd_data <- keep_genes(copd_genes, results_total$data)
results_copd <- analyze_covariate_lists(copd_data, results_total$phenotype, covs_pulm)




# Analysis
# Filtering with CV
# set data_to_analyze; set data_phenotype
data_to_analyze <- data_filtered$cv; data_phenotype <- 























####################################################
############# OLD CODE ############################
####################################################

####################################
####################################
# Run the linear models for data_to_analyze
#1results <- analyze(data_to_analyze, data_phenotype)

# Extract the models and p-values from the results
#2models <- results[seq(1,10,2)]
#3ps <- results[seq(2,10,2)]


# fdr correct the p-values from the linear models
#4fdrs <- fdr_correct(ps)
####################################
####################################



# have run the below analysis but have not re-coded it for this function-based version
# it was run as part of the v00 script
# STOPPED RECODING HERE FOR NOW: 2:57, monday november 25th
# 
# find the means for each of the groups (1/0) for the genes of interest
#model_COPDCA_fdr <- p.adjust(model_COPDCA_ps, method="fdr")
#fdr05_ind <- which(model_COPDCA_fdr < 0.05)
#num_genes_interest <- length(fdr05_ind)

#names_genes_interest <- row.names(data_to_analyze)[fdr05_ind]
#COPDCA_interesting <- data.frame(row.names=names_genes_interest)
#COPDCA_interesting$ind <- fdr05_ind
#COPDCA_interesting$p_vals <- model_COPDCA_ps[fdr05_ind]
#COPDCA_interesting$fdr_vals <- model_COPDCA_fdr[fdr05_ind]





#intx_ps <- models$COPD_CA_interaction_ps[,3]
#interesting <- which(intx_ps < 0.0005)
#interesting_loose <- which(intx_ps < 0.005) # gets 265 genes for the interaction term
#
#aics <- matrix(ncol=3,nrow=length(interesting))
#
#for(gene in 1:length(interesting))
#{
#	print(interesting[gene])
#	aics[gene,1] = AIC(model_COPD_CA[[interesting[gene]]])
#	aics[gene,2] = AIC(model_COPD_CA_interaction[[interesting[gene]]])
#	
#	aics[gene,3] = aics[gene,1] > aics[gene,2]
#
#
#}
#
#
#aics2 <- matrix(ncol=3, nrow=length(intx_ps))
#for(gene in 1:length(intx_ps))
#{
#	#print(gene)
#	#print(intx_ps[gene])
#	aics2[gene,1] = AIC(model_COPD_CA[[gene]])
#	aics2[gene,2] = AIC(model_COPD_CA_interaction[[gene]])
#	
#	aics2[gene,3] = aics2[gene,1] > aics2[gene,2]
#}
#
#aics_loose <- matrix(ncol=3, nrow=length(interesting_loose))
#for(gene in 1:length(interesting_loose))
#{
#	aics_loose[gene,1] = AIC(model_COPD_CA[[interesting_loose[gene]]])
#	aics_loose[gene,2] = AIC(model_COPD_CA_interaction[[interesting_loose[gene]]])
#	aics_loose[gene,3] = aics_loose[gene,1] > aics_loose[gene,2]
#
#}
#
## There are 3928 genes 'better explained' with interaction model, based on AIC values of the
## two different models. Seems to be mostly explained by the p-values (or at least related)
## i.e. the p-values for the interaction models, of the models better explained by the 
## interaction model are the smallest p-values from the set (the p value ranges are non-overlapping)
#intrx_better_ind <- which(aics2[,3] > 0)
#intrx_better_bol <- aics2[,3] > 0
#
#
## need to remove the NA cancer and COPD patients, of which there are 3
## need to figure out how to work with the exprSet data type
#nas_removed_ind <- !(is.na(phenotype.data$CANCER) | is.na(phenotype.data$COPD) | is.na(phenotype.data$RIN))
#expression.data.na <- expression.data[,nas_removed_ind]
#phenotype.data.na <- phenotype.data[nas_removed_ind,]
#
## number of patients with copd and cancer
#length(which(phenotype.data.na$CANCER==1 & phenotype.data.na$COPD==1))
## number of patients with cancer but not copd
#length(which(phenotype.data.na$CANCER==1 & phenotype.data.na$COPD==0))
## number of patients with copd but not cancer
#length(which(phenotype.data.na$CANCER==0 & phenotype.data.na$COPD==1))
## number of patients with neither copd nor cancer
#length(which(phenotype.data.na$CANCER==0 & phenotype.data.na$COPD==0))
#
## define a looser p value with which to pull 'interesting' genes from the 
#interesting_loose <- which(intx_ps < 0.005)
##
#
#
#colnames(copd_results) <- c("estimate_inter", "estimate_copd", "estimate_cancer", "estimate_sect", "std_error_inter", "std_error_copd", "std_error_cancer", "std_error_sect", "t_inter", "t_copd", "t_cancer", "t_sect", "p_inter", "p_copd", "p_cancer", "p_sect")
#copd_results_data <- data.frame(copd_results)
#
## CORRECT FOR MULTIPLE COMPARISONS USING FDR
## HOW TO CORRECT FOR MULTIPLE MODELS???
#print("Correcting for multiple comparisons...")
#copd_results_fdr_corrected <- p.adjust(copd_results_data$p_copd, method="fdr")
#
#
#copd_smk = table(phenotype.data$SMK, phenotype.data$COPD)
#chi_copd_smk = chisq.test(copd_smk)
#
#copd_can = table(phenotype.data$COPD, phenotype.data$CANCER)
#chi_copd_can = chisq.test(copd_cancer)
#
#can_smk = table(phenotype.data$SMK, phenotype.data$SMK)
#chi_can_smk = chisq.test(can_smk)
## Smoking and COPD are correlated - no way!
#
#
## FIND GENES OF INTEREST
## DIFFERENTIATE BY GENES UP AND GENES DOWN REGULATED
#pos_sect_ind = which(copd_results_data$t_sect > 0)
#neg_sect_ind = which(copd_results_data$t_sect < 0)
#
#t_sect_fdr_pos = copd_results_fdr_corrected[pos_sect_ind]
#t_sect_fdr_neg = copd_results_fdr_corrected[neg_sect_ind]
#
#
## PRODUCE HEAT MAPS
#
#
### CREATE .RNK FILE FOR GSEA UP AND DOWN REGULATION
### ENTREE for _at, will need to manually download
#

# Run SVA to find covariates of interest by correlating SVA variables to possible clinical variables
#cat("Running SVA to find possible covariates of interest: model COPD + CA + COPD*CA\n")
#sva_obj <- return_sva(data_expression, data_phenotype)

# Now find out if the variables from SVA are correlated with any of the phenotype traits
#sva_lm_ps <- matrix(nrow=sva_obj$n.sv, ncol=dim(data_phenotype)[2])
#for(sva_var in 1:sva_obj$n.sv)
#{
#	for(trait in 1:length(data_phenotype))
#	{
#		mod1 <- lm(sva_obj$sv[,sva_var] ~ data_phenotype[,trait])
#		sva_lm_ps[sva_var,trait] <- summary(mod1)$coefficients[2,4]
		#sva_lm[sva_var,trait] <- cor(sva_obj$sv[,sva_var], as.numeric(data_phenotype[,trait]))
#	}
#}

# Want a vector of length(number of traits in phenotype) with a 1 if correlated with a sv, 0 otherwise
#sva_lm_ps_fdr <- p.adjust(sva_lm_ps, method="fdr")
#sva_lm_ps_fdr <- matrix(sva_lm_ps_fdr, nrow=38, ncol=10, byrow=FALSE)
#cor_trait <- numeric(10)
#for(trait in 1:length(data_phenotype))
#{
#	num <-  length(which(sva_lm_ps_fdr[,trait] < 0.05))
#	if(num > 0)
#	{
#		cor_trait[trait] <- 1
#	}
#	else
#	{
#		cor_trait[trait] <- 0
#	}
#}
# all the traits correlate with at least one of the surrogate variables - oi!
setwd(initial_directory)


