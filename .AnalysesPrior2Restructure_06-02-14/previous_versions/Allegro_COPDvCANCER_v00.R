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
# What else?

# FOR EACH ANALYSIS, SAVE IN A SEPARATE DIRECTORY:
#	- separate script
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

#########################################################
############ DEFINE GLOBAL VARIABLES  ###################
#########################################################

# SAVE INITIAL DIRECTORY
init_directory <- getwd()

# DEFINE DATA DIRECTORY
data_directory <- "/protected/projects/pulmarray/Biollegro/RDS"

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
	# LOAD DATA

}

prep_data(data_to_analyze)
{
	# COERCE DATA TO APPROPRIATE TYPES

}



analyze <- function(data_to_analyze)
{	
	# 1. 

	# 2. 
	
	# 3. 

}


#########################################################
################# START ANALYSIS ########################
#########################################################

# store initial working directory; return to initial working directory at script end
initial <- getwd()
setwd("/protected/projects/pulmarray/Biollegro/RDS")
dataloc <- getwd()

print("Reading in data...")
data <- readRDS("Bronc_708_First_CEL.rds")
phenotype.data <- pData(data)
expression.data <- exprs(data)
num_genes <- dim(expression.data)[1]

# COERCE categorical phenotypes into factors
# Column 1 - BATCH
phenotype.data$BATCH <- factor(phenotype.data$BATCH)
# Column 3 - CANCER
phenotype.data$CANCER <- factor(phenotype.data$CANCER)
# Column 4 - SMK
phenotype.data$SMK <- factor(phenotype.data$SMK)
# Column 5 - COPD
phenotype.data$COPD <- factor(phenotype.data$COPD)
# Column 6 - GENDER
phenotype.data$GENDER <- factor(phenotype.data$GENDER)

# Create a new indicator function such that I = 1 when COPD=1 and CA=1, 0 otherwise
phenotype.data$COPDCA = factor(phenotype.data$COPD == 1 & phenotype.data$CANCER == 1)


# FILTER DATA
print("Filtering data...")
  # DIFFERENT FILTERING METHODS

  # 1) BY COEFFICIENT OF VARIANCE

	cv_filt_data <- remove_percent(expression.data, .25, coef_of_var)
	cv_mn_filt_data <- remove_percent(cv_filt_data, .25, mean)

  # 2) BY MEAN
  
  	mn_filt_data <- remove_percent(expression.data, .25, mean)
  	mn_cv_filt_data <- remove_percent(mn_filt_data, .25, coef_of_var)
  	
  # 3) BY SUBJECT EXPRESSION LEVELS
  
  	# STILL TO IMPLEMENT...

num_genes_filt_mn <- dim(mn_filt_data)[1]
num_genes_filt_cv <- dim(mn_cv_filt_data)[1]  	

# RUN LINEAR MODELS

# set the data (filtered or unfiltered from above) to use in the models
#	1. cv_filt_data		filtered by coefficient of variance
#	2. cv_mn_filt_data	filtered by cv and then by mean
#	3. mn_filt_data		filtered by mean
#	4. mv_cv_filt_data	filtered by mn and then by cv
#	5. expression.data	unfiltered data (all genes included)

data_to_analyze = cv_filt_data
num_genes = dim(data_to_analyze)[1]

# create a list to store the results from the linear models
#copd_results <- c()

# create space to hold the multiple models for each gene
model_COPD <- vector("list", num_genes)
model_CA <- vector("list", num_genes)
model_COPD_CA <- vector("list", num_genes)
model_COPD_CA_interaction <- vector("list", num_genes)
model_COPDCA <- vector("list", num_genes)

# create space to hold the extracted p-values from the gene models
model_COPD_ps <- matrix(nrow=num_genes, ncol=1)
model_CA_ps <- matrix(nrow=num_genes, ncol=1)
model_COPD_CA_ps <- matrix(nrow=num_genes, ncol=2)
model_COPD_CA_interaction_ps <- matrix(nrow=num_genes, ncol=3)
model_COPDCA_ps <- matrix(nrow=num_genes, ncol=1)

print("Modeling data")
for(gene in 1:num_genes)
{	
	gene_name <- rownames(data_to_analyze)[gene]
	cat("Running models for gene #", gene, " name:", gene_name, "\n")
	
	# COPD MODEL
	model_COPD[[gene]] <- lm(data_to_analyze[gene,] ~ COPD, na.action="na.exclude", data=phenotype.data)
	model_COPD[[gene]]$gene_name <- gene_name
	
	# CANCER MODEL
	model_CA[[gene]] <- lm(data_to_analyze[gene,] ~ CANCER, na.action="na.exclude", data=phenotype.data)
	model_CA[[gene]]$gene_name <- gene_name
	
	# COPD + CANCER MODEL
	model_COPD_CA[[gene]] <- lm(data_to_analyze[gene,] ~ COPD + CANCER, na.action="na.exclude", data=phenotype.data)
	model_COPD_CA[[gene]]$gene_name <- gene_name
	
	# COPD + CANCER + COPD*CANCER MODEL
	model_COPD_CA_interaction[[gene]] <- lm(data_to_analyze[gene,] ~ COPD + CANCER + COPD*CANCER, na.action="na.exclude", data=phenotype.data)
	model_COPD_CA_interaction[[gene]]$gene_name <- gene_name



	# SOME OLD CODE - undecided as to whether to delete or not
	#models[[gene]] <- lm(expression.data[gene,] ~ phenotype.data$COPD*phenotype.data$CANCER, na.action="na.exclude")
	#models[[gene]]$gene_name <- rownames(expression.data)[gene]
	#gene_temp_info <- c(summary(models[[gene]])$coefficients[,1], summary(models[[gene]])$coefficients[,2],summary(models[[gene]])$coefficients[,3], summary(models[[gene]])$coefficients[,4])
	#copd_results <- rbind(copd_results, gene_temp_info)
	#rownames(copd_results)[gene] <- rownames(expression.data)[gene]


	# Extract the p-values from each of the different models
	model_COPD_ps[gene] <- summary(model_COPD[[gene]])$coefficients[2,4]
	model_CA_ps[gene,] <- summary(model_CA[[gene]])$coefficients[2,4]
	model_COPD_CA_ps[gene,] <- summary(model_COPD_CA[[gene]])$coefficients[2:3,4]
	model_COPD_CA_interaction_ps[gene,] <- summary(model_COPD_CA_interaction[[gene]])$coefficients[2:4,4]

	# models to consider
	# model: gene_expression ~ COPD + RIN(maybe) 
	# model: gene_expression ~ RIN (to consider the genes most affected by RIN and remove them)
		# can at least correct for RIN
		
}

for (gene in 1:num_genes)
{
	gene_name <- rownames(data_to_analyze)[gene]
	cat("Running models for gene #", gene, " name:", gene_name, "\n")

	# COPDCA MODEL - uses COPDCA as the only predictor; 1 is COPD&CA, 0 is everything else
	model_COPDCA[[gene]] <- lm(data_to_analyze[gene,] ~ COPDCA, na.action="na.exclude", data=phenotype.data)
	model_COPDCA[[gene]]$gene_name <- gene_name
	model_COPDCA_ps[gene,] <- summary(model_CA[[gene]])$coefficients[2,4]
}


# find the means for each of the groups (1/0) for the genes of interest
model_COPDCA_fdr <- p.adjust(model_COPDCA_ps, method="fdr")
fdr05_ind <- which(model_COPDCA_fdr < 0.05)
num_genes_interest <- length(fdr05_ind)

names_genes_interest <- row.names(data_to_analyze)[fdr05_ind]
COPDCA_interesting <- data.frame(row.names=names_genes_interest)
COPDCA_interesting$ind <- fdr05_ind
COPDCA_interesting$p_vals <- model_COPDCA_ps[fdr05_ind]
COPDCA_interesting$fdr_vals <- model_COPDCA_fdr[fdr05_ind]





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


intx_ps <- model_COPD_CA_interaction_ps[,3]
interesting <- which(intx_ps < 0.0005)

aics <- matrix(ncol=3,nrow=length(interesting))

for(gene in 1:length(interesting))
{
	print(interesting[gene])
	aics[gene,1] = AIC(model_COPD_CA[[interesting[gene]]])
	aics[gene,2] = AIC(model_COPD_CA_interaction[[interesting[gene]]])
	
	aics[gene,3] = aics[gene,1] > aics[gene,2]


}



colnames(copd_results) <- c("estimate_inter", "estimate_copd", "estimate_cancer", "estimate_sect", "std_error_inter", "std_error_copd", "std_error_cancer", "std_error_sect", "t_inter", "t_copd", "t_cancer", "t_sect", "p_inter", "p_copd", "p_cancer", "p_sect")
copd_results_data <- data.frame(copd_results)

# CORRECT FOR MULTIPLE COMPARISONS USING FDR
# HOW TO CORRECT FOR MULTIPLE MODELS???
print("Correcting for multiple comparisons...")
copd_results_fdr_corrected <- p.adjust(copd_results_data$p_copd, method="fdr")


copd_smk = table(phenotype.data$SMK, phenotype.data$COPD)
chi_copd_smk = chisq.test(copd_smk)

copd_can = table(phenotype.data$COPD, phenotype.data$CANCER)
chi_copd_can = chisq.test(copd_cancer)

can_smk = table(phenotype.data$SMK, phenotype.data$SMK)
chi_can_smk = chisq.test(can_smk)
# Smoking and COPD are correlated - no way!


# FIND GENES OF INTEREST
# DIFFERENTIATE BY GENES UP AND GENES DOWN REGULATED
pos_sect_ind = which(copd_results_data$t_sect > 0)
neg_sect_ind = which(copd_results_data$t_sect < 0)

t_sect_fdr_pos = copd_results_fdr_corrected[pos_sect_ind]
t_sect_fdr_neg = copd_results_fdr_corrected[neg_sect_ind]


# PRODUCE HEAT MAPS


## CREATE .RNK FILE FOR GSEA UP AND DOWN REGULATION
## ENTREE for _at, will need to manually download

setwd(initial)


