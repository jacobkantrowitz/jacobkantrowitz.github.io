
#########################################################
############## LOAD NECESSARY LIBRARIES #################
#########################################################

# load bioinformatics packages
library(affy)

# set package library to load heatmap3
.libPaths(c(sprintf("/unprotected/projects/cbmhive/R_packages/R-%s", getRversion()), .libPaths()))
library(heatmap3)
#library(sva)
#library(sm)

#########################################################
############ DEFINE GLOBAL VARIABLES  ###################
#########################################################

# DEFINE CUTOFFS TO USE FOR FDR AND P-VALUE SIGNIFICANCE
f_cutoffs <- c(0.25, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
p_cutoffs <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)

CORRECTED <- 1
UNCORRECTED <- 0

ROWS <- 1
COLUMNS <- 2


# Variables that will be needed for analysis
# directory where the expressionSet RDS exists
dataDir <- ""
# filename of the expressionSet RDS file
dataFileName <- ""

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

# New Methods to Write so my scripts work with expressionSet data (biobase)
# setDataDir - set the directory where the expressionSet RDS exists
setDataDir <- function(pathName)
{
	dataDir <<- pathName
}

# setDataFilename - set the filename of the expressionSet RDS file
setDataFileName <- function(fileName)
{
	dataFileName <<- fileName
}

# loadData - load the expressionSet found in the RDS file
# check that it is an expressionSet object before returning it
loadExpressionSet <- function(pathName=dataDir, fileName=dataFileName)
{
	cat("Loading data............... ", pathName, fileName, "\n", sep="")
	tmpName <- paste(pathName, fileName, sep="")
	tmpObj <- readRDS(tmpName)
	
	cat("Data loaded successfully\nChecking object type\n")
	if(class(tmpObj)=="ExpressionSet")
	{
		cat("Object is of class ExpressionSet\n")
		return(tmpObj)
	}
	else
	{
		cat("ERROR: Object is NOT of class ExpressionSet\n\nSTOP\n")
	}
}

# return the number of samples from an expressionSet object (biobase)
sampleNumber <- function(exprSet)
{
	return(dim(exprSet)[2])
}

traitNumber <- function(exprSet)
{
	return(length(varLabels(exprSet)))
}

featureNumber <- function(exprSet)
{
	return(dim(exprSet)[1])
}

coercePhenotypeFields <- function(exprData, type)
{
	# COERCE DATA TO APPROPRIATE TYPES
	lbls <- varLabels(exprData)
	exprData2 <- exprData
	
	# for each phenotype fields in exprData
	for(i in 1:traitNumber(exprData))
	{
		# if the label being checked is the same label in the table at this index
		if(lbls[i]==type[i,2])
		{
			# if the class of the label being checked is different than what it should be
			if(class(pData(exprData)[,i])!=type[i,3])
			{
				# change the class of the label being checked to what it should be
				
				if(type[i,3]=="factor")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to factor\n", sep="")
					pData(exprData2)[,i] <- as.factor(pData(exprData2)[,i])
				}
				else if(type[i,3]=="numeric")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to numeric\n", sep="")
					pData(exprData2)[,i] <- as.numeric(pData(exprData2)[,i])
				}
				else if(type[i,3]=="integer")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to integer\n", sep="")
					pData(exprData2)[,i] <- as.integer(pData(exprData2)[,i])
				}
				else if(type[i,3]=="character")
				{
					cat("Changing trait ", i, "-", as.character(type[i,2]), " to character\n", sep="")
					pData(exprData2)[,i] <- as.character(pData(exprData2)[,i])
				}
				else
				{
					cat("ERROR: no class type name: ", type[i,3], "\n", sep="")
				}
				
			}

		}
	
	}
	return(exprData2)
}	

# removePercentES takes an expressionSet, a percentage (0.XX) and a function for filtering
# the gene data from the expressionSet. The genes are sorted in ascending order based
# on the method provided and then the bottom XX% are removed. The expressionSet object is
# returned with the genes removed from the expression matrix.
removePercentES <- function(exprData, percent, method)
{
	# calculate the function (e.g. mean, coefficient of variance) for each gene
	methodApplyData <- apply(exprs(exprData),1,method)
	
	# sort and index the resulting function (e.g. mean) values
	sortMethodData <- sort(methodApplyData, index.return = TRUE)
	
	# find the index in the sorted indices for the cutoff
	# i.e. if removing 25% of 10,000 genes, then the cutoff would be 2500
	minInd <- percent*length(sortMethodData$ix)
	
	# find the index of the gene that will be last to be included (i.e. last in the sort)
	# i.e. if have 10,000 genes then this will be 10,000
	maxInd <- length(sortMethodData$ix)
	
	# pull the subset of gene indices that will be kept
	# i.e. if keeping 2500/10000 genes these will be sorted.index[2500:10000]
	keepInd <- sortMethodData$ix[minInd:maxInd]
	
	# subset the original expressionSet and include only those genes above the cutoff,
	# removing the desired percentage of the genes based on the sort method provided
	exprs(exprData) <- exprs(exprData)[keepInd,]
	
	#keep_data$ix <- sort_method_data$ix[min_ind:max_ind]
	#keep_data$x <- sort_method_data$x[min_ind:max_ind]
	
	return(exprData)

}

# calculates the coefficient of variance
coefficient.of.variance <- function(x)
{
	# coefficient.of.variance is a simple function that returns the coefficient of variance
	# for a variable X (e.g. a data frame)

	sd(x) / abs(mean(x))
}


filteringES <- function(dataToFilter, percent1=0.25, percent2=0.25)
{
	cat("Filtering gene expression data...............\n")
	
	dataFiltered <- dataToFilter
	# MULTIPLE FILTERING METHODS

	# 1. BY COEFFICIENT OF VARIANCE
	cv <- removePercentES(dataToFilter, percent1, coefficient.of.variance)
	cv.mn <- removePercentES(cv, percent2, mean)

	# 2. BY MEAN
  	mn <- removePercentES(dataToFilter, percent1, mean)
  	mn.cv <- removePercentES(mn, percent2, coefficient.of.variance)
  	
	# 3. BY SUBJECT EXPRESSION LEVELS
  	# STILL TO IMPLEMENT
	
	cat("Completed filtering\n")
	
	# Return unfiltered data along with differently filtered data sets
	return(list("unfiltered"=dataToFilter, "cv"=cv, "cv.mn"=cv.mn, "mn"=mn, "mn.cv"=mn.cv))
}


# remove_percent takes a set of data, a percentage (0.XX) and a function for filtering
# the dataset; the data set is sorted in ascending order based on the provides filtering
# function and then the bottom XX% are removed; the expression dataset is returned 
# with the XX% removed
#remove_percent <- function(expr, percent, method)
#{
#
#	method_data <- apply(expr,1,method)
#	sort_method_data <- sort(method_data, index.return = TRUE)
#	min_ind <- percent*length(sort_method_data$ix)
#	max_ind <- length(sort_method_data$ix)
#	
#	keep_ind <- sort_method_data$ix[min_ind:max_ind]
#	keep_data <- expr[keep_ind,]
#	
#	#keep_data$ix <- sort_method_data$ix[min_ind:max_ind]
#	#keep_data$x <- sort_method_data$x[min_ind:max_ind]
#	
#	return(keep_data)
#
#}





save_analysis <- function(models)
{
	# FOR EACH ANALYSIS, SAVE IN A SEPARATE DIRECTORY:
	#	- separate script
	#	- list of genes of interest with p-values, fdr-values, higher/lower indicator
	#	- heatmap of genes of interest
	#	- readme file with filtering method, number genes analyzed, what else?


}


#init <- function()
#{
#	# SET THE WORKING DIRECTORY TO THE DATA DIRECTORY
#	setwd(data_directory)
#	
#	# PRINT CURRENT OPERATIONS
#	cat("Loading data...", data_directory, data_filename, "\n")
#	
#	# LOAD DATA
#	data_loaded <- readRDS(data_filename)
#	
#	# EXTRACT THE PHENOTYPES
#	data_phenotype <- pData(data_loaded)
#	
#	# EXTRACT THE EXPRESSION LEVELS
#	data_expression <- exprs(data_loaded)
#	
#	# DEFINE THE NUMBER OF GENES IN THE EXPRESSION DATA
#	assign("number_of_genes", dim(data_expression)[1], envir=.GlobalEnv)
#
#	# SET THE OBJECTS TO RETURN AS ATTRIBUTES OF THE data_loaded OBJECT
#	attributes(data_loaded) <- list("expression"=data_expression, "phenotype"=data_phenotype)
#
#	return(data_loaded)
#}

# New Methods to Write so my scripts work with expressionSet data (biobase)
# setDataDir - set the directory where the expressionSet RDS exists
# setDataFilename - set the filename of the expressionSet RDS file
# loadData - load the expressionSet found in 

# Variables that will be needed for analysis
# dataDir
# dataFilename


#prepare_data <- function(data_to_analyze)
#{
#	# COERCE DATA TO APPROPRIATE TYPES
#	
#	# Column 1 - BATCH
#	data_to_analyze$BATCH <- factor(data_to_analyze$BATCH)
#	
#	# Column 3 - CANCER
#	data_to_analyze$CANCER <- factor(data_to_analyze$CANCER)
#	
#	# Column 4 - SMK
#	data_to_analyze$SMK <- factor(data_to_analyze$SMK)
#	
#	# Column 5 - COPD
#	data_to_analyze$COPD <- factor(data_to_analyze$COPD)
#	
#	# Column 6 - GENDER
#	data_to_analyze$GENDER <- factor(data_to_analyze$GENDER)
#
#	# Create a new indicator function such that I = 1 when COPD=1 and CA=1, 0 otherwise
#	data_to_analyze$COPDCA <- factor(data_to_analyze$COPD == 1 & data_to_analyze$CANCER == 1)
#	
#	# Create a new indicator function that I = 1 when COPD=0 and CA=0, 0 otherwise
#	data_to_analyze$HEALTHY <- factor(data_to_analyze$COPD == 0 & data_to_analyze$CANCER == 0)
#	
#	data_to_analyze$indicator <- 0
#	data_to_analyze$indicator[which(data_to_analyze$CANCER==1)] <- 1
#	data_to_analyze$indicator[which(data_to_analyze$COPD==1)] <- 2
#	data_to_analyze$indicator[which(data_to_analyze$COPDCA==TRUE)] <- 3
#	data_to_analyze$indicator <- factor(data_to_analyze$indicator)
#	
#
#	return(data_to_analyze)
#}


#filtering <- function(data_to_filter)
#{
#	cat("Filtering data...\n")
#	
#	data_filtered <- data_to_filter
#	# MULTIPLE FILTERING METHODS
#
#	# 1. BY COEFFICIENT OF VARIANCE
#	cv <- remove_percent(data_to_filter, .25, coef_of_var)
#	cv_mn <- remove_percent(cv, .25, mean)
#
#	# 2. BY MEAN
#  	mn <- remove_percent(data_to_filter, .25, mean)
#  	mn_cv <- remove_percent(mn, .25, coef_of_var)
#  	
#	# 3. BY SUBJECT EXPRESSION LEVELS
#  	# STILL TO IMPLEMENT
#
#	# SET THE FILTERED SETS AS ATTRIBUTES OF UNFILTERED SET
#	attributes(data_filtered) <- list("cv"=cv, "cv_mn"=cv_mn, "mn"=mn, "mn_cv"=mn_cv)
#	
#	# RETURN UNFILTERED DATA WITH DIFFERENTLY FILTERED SETS AS ATTRIBUTES
#	return(data_filtered)
#}


return_lm <- function(data_temp, phenotype, covariates)
{
	number_of_genes <- dim(data_temp)[1]
	number_of_covs <- length(covariates)
	
	# CREATE SPACE FOR MODELS 
	models <- vector("list", number_of_genes)

	# CREATE SPACE FOR EXTRACTED P-VALUES FROM GENE MODELS
	number_of_levels <- calculate_total_levels(phenotype, covariates)
	model_ps <- matrix(nrow=number_of_genes, ncol=number_of_levels)
	
	# CREATE SPACE FOR MODEL RESIDUALS
	model_resids <- data_temp
	model_resids[,] <- 0

	# RUN MODELS
	model <- paste(covariates, collapse=" + ")
	cat("Running one model:", model, "\nRunning model for", number_of_genes, "genes\n")

	# Print status bar
	cat(rep("*", 50),"\n", sep="")
	
	for(gene in 1:number_of_genes)
	{	
		gene_name <- rownames(data_temp)[gene]
		#cat("Running models for gene #", gene, " name:", gene_name, "\n")
		print_status(gene, number_of_genes)
		
		models[[gene]] <- lm(as.formula(paste("data_temp[", gene, ",] ~", model, sep="")), data=phenotype)
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
analyze_covariate_lists <- function(data_temp, phenotype, covariates)
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
			results <- return_lm(data_temp, phenotype, covariates)
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
				temp <- return_lm(data_temp, phenotype, covariates[[model]])
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

analyze_lists <- function(listBool, data, phen)
{
	# listBool is a list with n items
	# for each of n items there are two items 
	# item 1: list of covariates to be sent to analyze_covariate_lists
	# item 2: boolean indicating whether or not to send item 1
	
	# create a list of covariate lists to analyze
	analyze <- list()

	cat("Received ",  length(listBool), " models\nChecking models\n", sep="")
	# for each item i, if it's boolean is TRUE, append it to a list of covariate lists to analyze
	for(analysis in listBool)
	{
		complete <- analysis[[2]]
		covars <- analysis[[1]]

		if(complete==TRUE)
		{
			analyze <- append(analyze, list(covars))
		}
	}
	
	cat("Running ", length(analyze), " of ", length(listBool), " models\n\n", sep="")
	if(length(analyze) > 0)
	{
		return(analyze_covariate_lists(data, phen, analyze))
	}
}


fdr_correct_matrix <- function(ps)
{
	temp <- ps
	for(column in 1:dim(ps)[2])
	{
		temp[,column] <- p.adjust(temp[, column], method="fdr")
	}
	return(temp)
}


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
	#mod1 <- model.matrix(~ 1 + indicator, data=data_phenotype)
	#mod1 <- model.matrix(~ 1 + COPD, data=data_phenotype)
	#mod0 <- model.matrix(~ 1, data=data_phenotype)
	
	mod1 <- model.matrix(~ 1 + BATCH + RIN + SMK + GENDER + AGE + COPD + CANCER, data=data_phenotype)
	mod0 <- model.matrix(~ 1 + BATCH + RIN + SMK + GENDER + AGE, data=data_phenotype)

	
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


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# STOPPED HERE FOR RECODING

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

keep_patients <- function(inds, data, phenotype)
{
	temp.data <- data[,inds]
	temp.phen <- phenotype[inds,]
	return(list(temp.data, temp.phen))


}

max_genes_heatmap <- 350
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
					heatmap3(model_object$data[gene_index,], col = bluered,  hclustfun=function(d) hclust(d, method="ward"),
					col.clustering = "semisupervised", ColSideColors = clabels, main = map_title)
					
					
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


generate_heatmap  <- function(inds, data, phen, rowClusters=NULL, mthd="average", colClus=phen$indicator,mn="Figure")
{
	clabels <- cbind(
		"Indicator" = copdca_colors[phen$indicator],
		"COPD Status" = copd_colors[phen$COPD],
		"Smoking Status" = smoking_colors[phen$SMK],
		"Cancer Status" = cancer_colors[phen$CANCER],
		"Gender" = gender_colors[phen$GENDER],
		"Batch" = batch_colors[phen$BATCH]
	)
	
	clabelsA <- cbind(
		"Cancer Status" = cancer_colors[phen$CANCER],
		"COPD Status" = copd_colors[phen$COPD],
		"Smoking Status" = smoking_colors[phen$SMK],
		"Gender" = gender_colors[phen$GENDER],
		"Batch" = batch_colors[phen$BATCH]
	)
	
	clabelsNEW <- cbind(
		"Indicator" = copdca_colors[phen$indicator],
		"Cancer Status" = cancer_colors[phen$FinalCaDXc],
		"COPD Status" = copd_colors[phen$COPD2_R7],
		"Smoking Status" = smoking_colors[phen$SMKc],
		"Gender" = gender_colors[phen$GENDERc],
		"Batch" = batch_colors[phen$BATCH]
	)
	
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
	heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabelsNEW, main = mn)

	if(!is.null(rowClusters))
	{
		heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, RowSideColors = rowClusters, main = mn)
	}

	# Uses ward for clustering; semi-supervised
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="ward"), col.clustering = "semisupervised", ColSideColors = clabelsA, main = "Figure")
	
	# MAIN - Uses average for clustering; semi-supervised
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "semisupervised", ColSideColors = clabels, main = "Figure")
	
	# Uses average for clustering; unsupervised
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")
	
	# Uses centroid for clustering; unsupervised
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="centroid"), col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")
	
	# Includes row labels when RowClusters is used; semi-supervised
	#heatmap3(data[inds,], col = bluered, hclustfun=function(d) hclust(d, method="average"), col.clustering = "semisupervised", ColSideColors = clabels, RowSideColors = rowClusters, main = "Figure")
}

# code to define clusters from heatmap dendrogram
#clabels <- cbind(
#	"Indicator" = copdca_colors[data_phenotype$indicator],
#	"COPD Status" = copd_colors[data_phenotype$COPD],
#	"Smoking Status" = smoking_colors[data_phenotype$SMK],
#	"Cancer Status" = cancer_colors[data_phenotype$CANCER],
#	"Gender" = gender_colors[data_phenotype$GENDER],
#	"Batch" = batch_colors[data_phenotype$BATCH]
#	)
#

return_cluster <- function(inds, data, phenotype, n.clusters=2, type=ROWS, mthd="average", colClus=phenotype$indicator,mn="Figure")
{
	# code to define clusters from heatmap dendrogram
	clabels <- cbind(
		"Indicator" = copdca_colors[phenotype$indicator],
		"COPD Status" = copd_colors[phenotype$COPD],
		"Smoking Status" = smoking_colors[phenotype$SMK],
		"Cancer Status" = cancer_colors[phenotype$CANCER],
		"Gender" = gender_colors[phenotype$GENDER],
		"Batch" = batch_colors[phenotype$BATCH]
		)
			
	pdf(NULL)
	# 1. Main - uses ward clustering; semi-supervised
	#res <- heatmap3(data_to_analyze[inds,], col = bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method="ward"), col.clustering = "semisupervised", ColSideColors = clabels, main = "Figure")
	
	# 2. uses average clustering; semi-supervised
	res <- heatmap3(data[inds,], col = bluered, keep.dendro=TRUE, hclustfun=function(d) hclust(d, method=mthd), col.clustering = colClus, ColSideColors = clabels, main = mn)
	
	# 3. unsupervised
	#res <- heatmap3(data_to_analyze[inds,], col = bluered, keep.dendro=TRUE, col.clustering = "unsupervised", ColSideColors = clabels, main = "Figure")	
	dev.off()
	
	# For columns, i.e. by patients
	if(type==COLUMNS)
	{
		clusters <- cutree(as.hclust(res$Colv), k=n.clusters)
	}
	# For rows, i.e. by genes
	if(type==ROWS)
	{
		clusters <- cutree(as.hclust(res$Rowv), k=n.clusters)
	}

	return(clusters)
}

fold.change <- function(temp, groups, g1=1, g2=2)
{
	g1.ind <- which(groups==g1)
	g2.ind <- which(groups==g2)
	num.genes <- nrow(temp)
	fcs <- matrix(nrow=num.genes, ncol=3)
	
	for(gene in 1:num.genes)
	{
		fcs[gene, 1] <- mean(temp[gene, g1.ind])
		fcs[gene, 2] <- mean(temp[gene, g2.ind])
		fcs[gene, 3] <- fcs[gene,1] - fcs[gene,2]
	}
	
	return(fcs)
}

save.genes <- function(inds, exprData, phen, filename, rClusters=NULL, mtd="average", cClus=phen$indicator,man="Figure")
{
	save_entrez(inds, rownames(exprData), paste(filename, ".txt", sep=""))
	pdf(paste(filename,".pdf", sep=""))
	generate_heatmap(inds, exprData, phen, rowClusters=rClusters, mthd=mtd, colClus=cClus, mn=man)
	dev.off()
	cat("Saved entrez IDs and saved heatmap for ", length(inds), " genes.\n", sep="")
}

save.packages <- function()
{
	cat("\nSaving packages\n")
	int <- installed.packages()
	dt <- as.character(Sys.Date())
	write.table(int, file=paste(dt,"_packages.csv", sep=""), sep=",")
	cat("Saved OK\n")
}


# Weighted Voting Algorithm
wv.model <- function (data, classlabel, correction = TRUE)
{
        # Implementation of weighted voting code from Golub et al., Science, 1999.
        # Adam Gower, 2008

	# INPUT
        # data: matrix or data frame that contains all probesets for weighted voting by rows and training samples by columns.
        # classlabel: vector of binary outputs (0 or 1).
        # correction: flag to use Broad Institute's correction method (see below).
	#
	# OUTPUT
	# A list with two elements:
	#     weights: a vector of weights for each probeset
	#     means:   a vector of means for each probeset across all samples
	#     The elements of both vectors are named according to the probeset names, i.e., the row names of the data matrix 

        y <- lapply(0:1, function (i) {data[,classlabel==i,drop=F]});                           # y[[1]] = class 0 (generally controls); y[[2]] = class 1 (generally cases).
        mu <- sapply(y, apply, 1, mean, na.rm=T, simplify=F);                                   # mean of each probeset across all samples in each class
        sigma <- sapply(y, apply, 1, sd, na.rm=T, simplify=F);                                  # standard deviation of each probeset across all samples in each class
        if (correction) {sigma <- mapply(pmax, sigma, lapply(mu, `*`, 0.2),SIMPLIFY=F)};        # Correction from Broad Institute: sigma = max(sigma, 0.2*mu);
                                                                                                #     i.e., %CV always >= 20% (communicated to me by Jen Beane)
        a <- (mu[[2]] - mu[[1]]) / (sigma[[1]] + sigma[[2]]);					# Compute the weights (signal to noise metric, similar to t statistic)
        g <- (mu[[1]] + mu[[2]]) / 2;								# Compute the means of each probeset across all samples
        return(list(weights=a, means=g));							# Return all computations in a named list
}

predict.wv <- function (model, data)
{
        # Implementation of weighted voting code from Golub et al., Science, 1999.
        # Adam Gower, 2008

	# INPUT
        # model: a list with two elements, named weights and means, that contains the weights and means for each probeset as generated by wv.model()
        # data: matrix or data frame that contains all probesets for weighted voting by rows and test samples by columns.
	#       NOTE: the row names of data must be named in the same manner as the weights and means in the model!
	#
	# OUTPUT
	# A list with three elements:
	#     scores:      a vector of weighted voting scores for each sample; negative = class 0, positive = class 1
	#     predictions: a vector of class predictions (0 or 1) for each sample as determined from the sign of the scores vector
	#     strengths:   a vector of prediction strengths for each sample
	
        indices <- match(names(model$means), rownames(data));                                   # Get indices of variables that are in the model
        votes <- model$weights * (data[indices,,drop=F]-model$means);                           # Determine votes: weight * (value - mean); samples now in columns
	
        # Create matrix V with samples by rows and two columns: the sums of all negative and positive weighted votes, respectively
        V <- sapply(c(-1,1), function (s) {                                                     # For each sign (negative or positive),
                 apply(votes, 2, function (x) {                                                 #     for each sample (column in votes matrix),
                         sum(x[which(sign(x)==s)])                                              #         sum votes for that sample by sign
                 })
             });

        scores <- V[,1] + V[,2];                                                                # Calculate scores: sum of negative votes and positive votes
        predictions <- as.numeric(scores > 0);                                                  # Make predictions: score < 0 = class 0; score > 0 = class 1
        strengths <- abs(scores) / (V[,2] - V[,1]);                                             # Prediction strengths: |scores| / sum |votes|
        return(list(scores=scores, predictions=predictions, strengths=strengths));           # Return all computations in a named list
}


#### Definitions of functions from the Wiki from Marc Lenburg
pauc.p <- function(d, class.vector, p = 0.1, n = 10)
{
	# d = data matrix with samples in columns
	# class.vector = class annotation for columns of d (see help for rowpAUCs function in genefilter for more info)
	# p = 1 - specificity value at which to determine pAUC (e.g. p = 0.1 is integral of sensitivity over the range of 90 - 100% specificity)
	# n = times to permute class labels and calculate pAUCs for each gene in dataset (used to estimate null distribution). So if the dataset has 20,000 genes, n = 10 will give you 200,000 pAUCs based on 10 different permutations of the true class label.
	require(genefilter)
	actual <- area(rowpAUCs(d, class.vector, p)) # this is how to calculate pAUC. This is the only useful thing in this function
	permuted <- NULL
	for (i in 1:n)
	{
		# this is a crude attempt to estimate the null distribution of pAUCs (globally rather than per gene)
		permuted <- c(permuted, area(rowpAUCs(d, sample(class.vector, length(class.vector)), p)))
		cat("*")
	}
	cat("\n")

	empiric.p <- function(val, random.values)
	{
		return(sum(random.values > val))
	}
	
	# the next step is not optimized for speed, but determines the empiric p-value for each observed pAUC based on the estimated null distribution
	p.val <- apply(cbind(actual), 1, function(x) empiric.p(x, permuted))
	p.val <- p.val / length(permuted)
	return(p.val)

}
#### End of definitions from Wiki from Marc Lenburg



