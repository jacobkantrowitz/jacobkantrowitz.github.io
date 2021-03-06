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


run_for_residuals <- function(data_to_analyze, phenotype)
{
	# run a model with all phenotype terms (SMK, AGE, GENDER, RIN, BATCH)
	# return the residuals from the models
	number_of_genes <- dim(data_to_analyze)[1]
	number_of_patients <- dim(data_to_analyze)[2]
	models <- vector("list", number_of_genes)
	resids <- data_to_analyze
	resids[,] <- 0
	
	# Run models
	cat("Running linear models...\n")
	
	for(gene in 1:number_of_genes)
	{
		gene_name <- rownames(data_to_analyze)[gene]
		#cat("Running models for gene #", gene, " name:", gene_name, "\n")
		print_status(gene, number_of_genes)

		models[[gene]] <- lm(data_to_analyze[gene,] ~ SMK + RIN + BATCH + GENDER + AGE, na.action="na.exclude", data=phenotype)
		models[[gene]]$gene_name <- gene_name
		resids[gene,] <- summary(models[[gene]])$residuals
	}
	
	return(resids)
}



analyze <- function(data_to_analyze, phenotype)
{

	number_of_genes <- dim(data_to_analyze)[1]
	
	# CREATE SPACE FOR GENE MODELS 
	COPD <- vector("list", number_of_genes)
	CA <- vector("list", number_of_genes)
	COPD_CA <- vector("list", number_of_genes)
	COPD_CA_interaction <- vector("list", number_of_genes)
	COPDCA <- vector("list", number_of_genes)

	# CREATE SPACE FOR EXTRACTED P-VALUES FROM GENE MODELS
	COPD_ps <- matrix(nrow=number_of_genes, ncol=1)
	CA_ps <- matrix(nrow=number_of_genes, ncol=1)
	COPD_CA_ps <- matrix(nrow=number_of_genes, ncol=2)
	COPD_CA_interaction_ps <- matrix(nrow=number_of_genes, ncol=3)
	COPDCA_ps <- matrix(nrow=number_of_genes, ncol=1)
	
	# RUN MODELS
	cat("Running linear models...\n")

	for(gene in 1:number_of_genes)
	{	
		gene_name <- rownames(data_to_analyze)[gene]
		#cat("Running models for gene #", gene, " name:", gene_name, "\n")
		print_status(gene, number_of_genes)
	
		# COPD MODEL
		COPD[[gene]] <- lm(data_to_analyze[gene,] ~ COPD, na.action="na.exclude", data=phenotype)
		COPD[[gene]]$gene_name <- gene_name
		COPD_ps[gene] <- summary(COPD[[gene]])$coefficients[2,4]

		# CANCER MODEL
		CA[[gene]] <- lm(data_to_analyze[gene,] ~ CANCER, na.action="na.exclude", data=phenotype)
		CA[[gene]]$gene_name <- gene_name
		CA_ps[gene] <- summary(CA[[gene]])$coefficients[2,4]
	
		# COPD + CANCER MODEL
		COPD_CA[[gene]] <- lm(data_to_analyze[gene,] ~ COPD + CANCER, na.action="na.exclude", data=phenotype)
		COPD_CA[[gene]]$gene_name <- gene_name
		COPD_CA_ps[gene,] <- summary(COPD_CA[[gene]])$coefficients[2:3,4]
	
		# COPD + CANCER + COPD*CANCER MODEL
		#COPD_CA_interaction[[gene]] <- lm(data_to_analyze[gene,] ~ COPD + CANCER + SMK + RIN + BATCH + GENDER + AGE + COPD*CANCER, na.action="na.exclude", data=phenotype)
		COPD_CA_interaction[[gene]] <- lm(data_to_analyze[gene,] ~ CANCER + COPD + SMK + RIN + BATCH + GENDER + AGE + COPD*CANCER, na.action="na.exclude", data=phenotype)
		COPD_CA_interaction[[gene]]$gene_name <- gene_name
		COPD_CA_interaction_ps[gene,] <- summary(COPD_CA_interaction[[gene]])$coefficients[2:4,4]

		# COPDCA MODEL - uses COPDCA as the only predictor; 1 is COPD&CA, 0 is everything else
		# if the goal is to look at the interaction between COPD and cancer, i.e. find out
		# what is present/changed only in this context, then possibly can collapse the 
		# groups to be COPD&CA vs everything else (i.e. CA, COPD, or nothing)
		COPDCA[[gene]] <- lm(data_to_analyze[gene,] ~ COPDCA, na.action="na.exclude", data=phenotype)
		COPDCA[[gene]]$gene_name <- gene_name
		COPDCA_ps[gene,] <- summary(CA[[gene]])$coefficients[2,4]		
		
		
		# models to consider
		# model: gene_expression ~ COPD + RIN(maybe) 
		# model: gene_expression ~ RIN (to consider the genes most affected by RIN and remove them)
			# can at least correct for RIN
			
		# model: COPDCA vs COPD
	}
	
	results <- list("COPD"=COPD, "COPD_ps"=COPD_ps,
								"CA"=CA, "CA_ps"=CA_ps,
								"COPD_CA"=COPD_CA, "COPD_CA_ps"=COPD_CA_ps,
								"COPD_CA_interaction"=COPD_CA_interaction, "COPD_CA_interaction_ps"=COPD_CA_interaction_ps,
								"COPDCA"=COPDCA, "COPDCA_ps"=COPDCA_ps)

	return(results)
}


fdr_correct <- function(ps_list)
{
	temp_copy <- ps_list
	number_of_ps <- length(ps_list)
	
	for(ps in 1:number_of_ps)
	{
		temp <- ps_list[[ps]]
		
		for(column in 1:dim(temp)[2])
		{
			temp[,column] <- p.adjust(temp[, column], method="fdr")
		}
		ps_list[[ps]] <- temp
	}

	return(ps_list)
}


print_status <- function(current, total)
{
	pct <- round(current/total * 100)
	cat(pct,"% of models run\n")
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

generate_cutoff_lists <- function(values)
{
	temp <- list()
	cutoffs <- c(0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
	for(i in cutoffs)
	{
		# NEED TO FINISH
	}

}





#########################################################
################# START ANALYSIS ########################
#########################################################


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

resids <- run_for_residuals(data_to_analyze, data_phenotype)

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





# pull out the p-values from all the different models
# put pvalues for each model into model-specific matrix
# added these lines to the above for-loop
#for(gene in 1:num_genes)
#{
#	model_COPD_ps[gene] <- summary(model_COPD[[gene]])$coefficients[2,4]
#	model_COPD_CA_ps[gene,] <- summary(model_COPD_CA[[gene]])$coefficients[2:3,4]
#	model_COPD_CA_interaction_ps[gene,] <- summary(model_COPD_CA_interaction[[gene]])$coefficients[2:4,4]
#	model_COPD_ps[gene,] <- summary(model_COPD[[gene]])$coefficients[2,4]
#}


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
setwd(initial_directory)


