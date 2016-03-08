# return_sva runs SVA based on an as-of-now (2014-06-03) predefined model that needs
# to be updated; this function ideally will accept a list of covariates as well as a single
# covariate as the covariate of interest, along with an expressionSet object

# NOTE: right now (2014-06-03) this function has no utility as the covariate names have
# all changed with the updated Allegro phenotype information (clinical annotations)
return_sva <- function(exprData)
{

	# There can be no NA values in the phenotype data used to create the models for SVA

	# Build models for inclusion in SVA
	#mod1 <- model.matrix(~ 1 + indicator, data=data_phenotype)
	#mod1 <- model.matrix(~ 1 + COPD, data=data_phenotype)
	#mod0 <- model.matrix(~ 1, data=data_phenotype)
	
	mod1 <- model.matrix(~ 1 + BATCH + RIN + SMK + GENDER + AGE + COPD + CANCER, data=pData(exprData))
	mod0 <- model.matrix(~ 1 + BATCH + RIN + SMK + GENDER + AGE, data=pData(exprData))

	
	# Construct the surrogate variables and return them
	
	sva_obj <- sva(exprs(exprData), mod1, mod0)
	
	return(sva_obj)
}