#covar_correction_v2.R
#By Becky Kusko April 2014, rkusko@bu.edu
#This program iteratively tests covariates
#Adds the smallest pval covariate to the next model
#Then tests all remaining covariates
#until p > 0.05 or there are no more covars

#setup
library(MASS)
#[1] MASS_7.3-31

#####################PRELIM STUFF##################
cancer_annotation <- read.delim(file='annotation.txt',header=T,as.is=T,check.names=F,sep='\t',colClasses=c('factor','factor','factor','factor','factor','numeric','numeric','numeric','numeric','numeric','numeric','factor','factor','numeric','numeric','numeric','numeric','numeric','numeric','factor','factor','numeric','numeric','numeric','numeric','numeric','numeric','factor','numeric','numeric','numeric','numeric','numeric','numeric','factor','factor','factor','factor','factor','factor','factor','factor','factor','factor','factor','factor','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric','factor','factor','factor','factor','factor','factor'))
####It's important that the variables are read in correctly, or else the logit will fail

# Make sure that main variable has only two levels
if(length(levels(data_cancer[['Cancer']])) > 2) {
  data_cancer[['Cancer']] <- factor(as.character(data_cancer[['Cancer']]))
}

#File telling us what to do with which variable
#vars$var is in the same order as colnames(data)
vars <- read.delim(file='which_test.txt',header=T,as.is=T)


anova.worker <- function(X,dataset) {
  tryCatch({
	#X is a vector of indices, should be in the order wanted in the model
	#dataset is the dataset to be tested on

	#Create formula for model
	curr_data <- dataset
	covariates <- paste("curr_data[,",X,"]")
	y <- paste("curr_data[[Cancer']]")  ######Here you need to write the name of your main variable of interest
	fla <- paste(y," ~",paste(covariates,collapse="+"))

	#run model, collect statistics and report back
	model <- glm(as.formula(fla),family=binomial)
  	summ <- anova(model,test='Chisq');
	testedIndex <- X[length(X)]
	testedvar <- colnames(dataset)[as.numeric(testedIndex)]
   	c(testedvar,testedIndex,summ[['Pr(>Chi)']][length(X)+1])
	#this needs to return covariate index too
  },error=function(e){c(NA,NA)})
}

pick.covars <- function(dataset,variables,outname){
#dataset is the annotation file to be tested on
#variables is a list of names of columns from the dataset to be tested

data_ind <- which(colnames(dataset) %in% variables)
extra <- c()
extra_names <- c()

#This loop will be broken when pval > 0.05
for (i in 1:length(variables)){
	#First test to see if any covars are significant
	if  (i == 1){
	out <- lapply(X=data_ind,FUN=anova.worker,dataset=dataset)
	pval_table <- do.call(rbind,out)
	pval_table_order <- pval_table[order(as.numeric(pval_table[,3])),]
	nextT <- pval_table_order[1,]
	RemainCovar <- data_ind[-which(data_ind == nextT[2])]	
	} 
	if(i > 1){
		#check to see if we can stop
		if (as.numeric(nextT[3]) > 0.05){
			write.table(pval_table,outname,sep="\t",quote=FALSE)
			return(c(extra_names,nextT[3]))
			break
		}else{
			extra <- c(extra,nextT[2])
			extra_names <- c(extra_names,nextT[1])
			pval_table <- c()
			for (j in 1:length(RemainCovar)){
				input <- c(extra,RemainCovar[j])
				out <- anova.worker(input,dataset)
				pval_table <- rbind(pval_table,out)			
			}
			nextT <- pval_table[order(as.numeric(pval_table[,3])),][1,]
			RemainCovar <- RemainCovar[-which(RemainCovar == nextT[2])]
		}

		}
	
}


}


#######################################################################
###############ACTUALLY PICK COVARS####################################
#######################################################################

all <- vars$var[which(vars$loop == "Y")]
data_covars <- pick.covars(data_cancer,all,"all_remaining_vars_Cancer.txt")

#This table will be the list of recommended covars
write.table(data_covars,"all_covars_cancer.txt",sep="\t",quote=FALSE)


