# Setup Script for COPD-Cancer interaction analysis
# Script for analysis of Allegro data; COPD vs CANCER
# Jake Kantrowitz, BUSM MD/PhD Candidate 2019, Molecular and Translational Med.
# jacob.kantrowitz@gmail.com, kantro@bu.edu
# (617)-584-2393
# 2014-06-03
# We are looking for an interaction between cancer and COPD in the Allegro dataset
# This is a continuation of the work I did over the summer of 2013 as a rotation student


##########################################################################################
##########################################################################################
##########################################################################################
# SOME NOTES:
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

##########################################################################################
##########################################################################################
##########################################################################################


#########################################################
################# START ANALYSIS ########################
#########################################################

cat("Start analysis\n\n")

if(exists("LOADED")==FALSE)
{
	LOADED <- FALSE
}


# new starting statement
if(LOADED==FALSE)
{

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
	
	# SAVE INITIAL DIRECTORY 
	initial_directory <- getwd()
	
	# DEFINE CUTOFFS TO USE FOR FDR AND P-VALUE SIGNIFICANCE
	f_cutoffs <- c(0.25, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)
	p_cutoffs <- c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
	
	CORRECTED <- 1
	UNCORRECTED <- 0
	
	ROWS <- 1
	COLUMNS <- 2
	
	#########################################################
	############ DEFINE COLOR SCHEMES  ######################
	#########################################################
	bluered <- colorRampPalette(c("blue","white","red"))(256)
	
	# Define color schemes, based on factors
	cancer_colors <- c("0"="black", "1"="red");
	copd_colors <- c("0"="bisque", "1"="blue");
	gender_colors <- c("0"="blue", "1"="pink");
	smoking_colors <- c("1"="grey", "2"="white");
	batch_colors <- c("12"="red", "13"="blue", "14"="orange", "15"="green", "16"="black");
	copdca_colors <- c("0"="green", "1"="yellow", "2"="orange", "3"="red")

	#########################################################
	################# START SETUP ###########################
	#########################################################

	source("../../functions/sourceDir.R")
	source("../../functions/loadFunctions.R")
	loadFunctions()
	# set the directory where the RDS expressionSet data is located
	setDataDir("/protected/projects/pulmarray/Biollegro/RDS/NewAnnot/")
	# set the filename of the RDS file holding the expressionSet object
	setDataFileName("Allegro_Bronch_PreQC_n884_NewAnnot_2014-05-27.rds")
	
	# load the expression data set
	exprData <- loadExpressionSet()
	
	# Define key phenotype fields that patients cannot be missing (i.e. NA)
	keyTraits <- c("RIN", "BATCH", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc")
	
	# remove patients who have at least one NA in the phenotype information keyTraits
	exprData <- cleanNAForAnalysis(exprData, keyTraits)
	
	# remove patients who are never smokers
	exprData <- exprData[, exprData$SMKc!=3]
		
	# coerce exprData phenotype fields to the appropriate class (e.g. factor, numeric, etc.)
	#newClasses <- "/protected/projects/pulmarray/Allegro/COPD_Cancer/tmpNewOrganization/newAnnotationClasses.txt"
	newClasses <- "../../other/newAnnotationClasses.txt"
	newClasses <- read.csv(newClasses, sep=",", head=TRUE)
	
	exprData <- coercePhenotypeFields(exprData, newClasses)

	# need to create filtered versions of the data
	exprDataFiltered <- filtering(exprData)
	
	# define number of genes remaining in once or twice filtered data sets
	geneNumberFilt1 <-  featureNumber(exprDataFiltered$cv)
	geneNumberFilt2 <-  featureNumber(exprDataFiltered$cv.mn)


	# change LOADED to TRUE so the above code is not rerun every time unnecessarily
	LOADED <- TRUE
	
	# return to the initial_directory so any saved items are appropriately located
	setwd(initial_directory)
	save.packages()

}

#for(i in 1:59)
#{
#	cat(i, varLabels(exprData2)[i], class(pData(exprData)[,i]), class(pData(exprData2)[,i]),sep="  ")
#	cat("\n")
#}