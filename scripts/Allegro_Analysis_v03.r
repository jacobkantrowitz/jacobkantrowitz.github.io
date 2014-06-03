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
############ DEFINE GLOBAL VARIABLES  ###################
#########################################################

# SAVE INITIAL DIRECTORY 
initial_directory <- getwd()

# DEFINE DATA DIRECTORY AND FILENAME
data_directory <- "/protected/projects/pulmarray/Biollegro/RDS"
data_filename <- "Bronch_708_First_CEL.rds"

# DEFINE THE NUMBER OF GENES INCLUDED IN CURRENT ANALYSIS
number_of_genes <- 0


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
	# set the directory where the RDS expressionSet data is located
	setDataDir("/protected/projects/pulmarray/Biollegro/RDS/NewAnnot/")
	# set the filename of the RDS file holding the expressionSet object
	setDataFileName("Allegro_Bronch_PreQC_n884_NewAnnot_2014-05-27.rds")
	
	# load the expression data set
	exprData <- loadExpressionSet()
	
	# Define key phenotype fields that patients cannot be missing (i.e. NA)
	keyTraits <- c("RIN", "BATCH", "AGEcalc", "GENDERc", "SMKc", "FinalCaDXc")
	# Find the column indices of the key phenotype fields in the expressionSet
	keyTraitsInd <- match(keyTraits, varLabels(exprData))

	# find indices of patients with key missing phenotype information
	remove.pts <- numeric()	
	for(t in keyTraitsInd)
	{
		remove.pts <- append(remove.pts, which(is.na(pData(exprData)[,t])))
	}
	
	# find indices of patients who are never smokers (smoking = 3)
	remove.pts <- append(remove.pts, which(exprData$SMKc==3))
	
	# remove any non-unique indices from remove.pts and sort numerically
	remove.pts <- sort(unique(remove.pts))
	
	# remove patients from expressionSet who are flagged above
	# subsetting expressionSets: exprData['range of features', 'range of samples']
	# create T/F vector for removing flagged samples
	remove.pts <- is.na(match(1:sampleNumber(exprData), remove.pts)) 
	exprData <- exprData[,remove.pts]
	
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