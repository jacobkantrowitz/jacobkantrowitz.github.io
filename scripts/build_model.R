## Josh Campbell
## 9-12-11
## Calculates predictive models for a set of phenotypes in a dataset using a variety of parameters

## Load in libraries and functions needed later
lib.loc <- "~/R/x86_64-unknown-linux-gnu-library/3.0/"
source("biomarker_functions.R")
library(limma)
library(knnflex)
library(class)
library(e1071)
library(affy)
# library(hugene10stv1hsentrezgcdf)
library(hugene10sthsentrezgcdf, lib.loc=lib.loc)
library(randomForest)
library(MASS)
library(glmnet, lib.loc=lib.loc)
library(corpcor, lib.loc=lib.loc)
library(sva, lib.loc=lib.loc)



## Read in command line arguments
for (e in commandArgs()) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		temp = ta[[1]][2]
		if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
			temp = as.integer(temp)
		}
		if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
			temp = as.numeric(temp)
		}
		assign(ta[[1]][1],temp)
	} else {
		assign(ta[[1]][1],TRUE)
	}
}

prefix <- "130926"

xval.iteration = as.numeric(xval.iteration)
seed <- xval.iteration + 1000
set.seed(seed)
print(paste("Seed:", seed))

## Outfile prefixes
file.out = paste(prefix, "_predictions_", xval.iteration, ".txt", sep="")
unlink(file.out)


## Read in datasets/demographics
pheno <- read.table("/protected/projects/pulmseq/Allegro/Bronch_microRNA/edrn128/combined/120127_phenotypes_combined.txt", check.names=F, sep="\t", stringsAsFactors=F,  header=T)

rownames(pheno) = pheno$CEL
pheno$Cancer <- as.factor(pheno$Cancer)
pheno$Sex <- as.factor(pheno$Sex)
pheno$Smk_Status <- as.factor(pheno$Smk_Status)
pheno$Batch <- as.factor(pheno$Batch)
pheno$random_Cancer <- as.factor(pheno$random_Cancer)
pheno$random_Sex <- as.factor(pheno$random_Sex)
pheno$random_Smk_Status <- as.factor(pheno$random_Smk_Status)
pheno$random_Control <- as.factor(pheno$random_Control)
pheno$Site <- as.factor(pheno$Site)

## Set up train/test sets
xval.percentage <- 0.8
xval.train.num <- ceiling(xval.percentage*nrow(pheno))
xval <- sample(rownames(pheno))

train.samples <- xval[1:xval.train.num]
test.samples <- xval[(xval.train.num+1):length(xval)]


## Record which samples were used for train/test sets
write.table(matrix(c(xval.iteration, "Train_samples:", train.samples, "Test_samples:", test.samples), nrow=1), paste(prefix, "_xval_index.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=T)

## CEL file location
cel.file.loc = "/unprotected/projects/pulmarray/CEL/Allegro_mRNA_bronch_batch12thru15"

## Set up parameter list
param.list = list(
		sample.filter=c("none", "RIN5", "Hx"),
		normalization=c("rma"),
		batch.correct=c("none", "Combat"),
		gene.filter=c(0.25, 0.50, 1),
                phenotype=c("Smk_Status", "Sex", "Cancer", "random_Smk_Status", "random_Sex", "random_Cancer"),
		feat.select=c("lm", "glmnet", "rf"),
                num.feat=c(10,20,30,40,50,60,70,80,90,100),
		model.predict=c("wv", "knn", "svm", "rf", "nb", "glmnet")
                )

## Set up the paramter index
index = calculate.index(param.list)
index <- index[order(index[,1], index[,2], index[,3], index[,4], index[,5], index[,6], index[,7], index[,8]),] 


## Run through all possible models
for(i in 1:nrow(index)) {
	
	param.index = index[i,]
	var.samples = param.list[[1]][param.index[1]]
	var.norm = param.list[[2]][param.index[2]]
	var.batch.correct = param.list[[3]][param.index[3]]
	var.gene.select = param.list[[4]][param.index[4]]
	var.phenotype = param.list[[5]][param.index[5]]
	var.feat.select = param.list[[6]][param.index[6]]
	var.num.feat = param.list[[7]][param.index[7]]
	var.model.predict = param.list[[8]][param.index[8]]
	
	min.change = 1
	if(i > 1) {
		param.previous = index[i-1,]
		min.change = min(which(param.index != param.previous))
	}
	
	##########################################################
	## Biomarker Building#####################################
	##########################################################
	
	
	## Choose samples ########################################
	if(min.change <= 1) {   
		if(var.samples == "RIN5") {
			xval.filter = rownames(pheno)[which(pheno$Bronch_RIN < 5)]
		}else if (var.samples == "Hx") {
			xval.filter = rownames(pheno)[which(pheno$Hx_cancer == 1)]
		}else if (var.samples == "High_Prob") {
			xval.filter = rownames(pheno)[which(pheno$High_prob == 2)]
		}else {
			xval.filter = NULL	
		}
	
		if (!is.null(xval.filter)) {
			xval.train.samples <- setdiff(train.samples, xval.filter)
			xval.test.samples <- setdiff(test.samples, xval.filter)
		} else {
			xval.train.samples <- train.samples
			xval.test.samples <- test.samples
		}
	}
	
	## Normalize
	if(min.change <= 2) {
		s <- c(xval.train.samples, xval.test.samples)
		
		if(var.norm == "rma") {	
			
			xval.data <- exprs(justRMA(filenames=s, celfile.path=cel.file.loc, cdfname="hugene10stv1hsentrezgcdf"))
#       xval.data <- read.table("/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model_allbatches/LUNGevity_log2cpm_perc15_filter.txt")
			
		}else if (var.norm == "rma_sep") {
			s.12 <- intersect(rownames(pheno)[which(pheno$Batch == 12)], s) 
			b12 <- exprs(justRMA(filenames=s.12, celfile.path=cel.file.loc, cdfname="hugene10stv1hsentrezgcdf"))
			s.13 <- intersect(rownames(pheno)[which(pheno$Batch == 13)], s) 
			b13 <- exprs(justRMA(filenames=s.13, celfile.path=cel.file.loc, cdfname="hugene10stv1hsentrezgcdf"))
			s.14 <- intersect(rownames(pheno)[which(pheno$Batch == 14)], s)
                        b14 <- exprs(justRMA(filenames=s.14, celfile.path=cel.file.loc, cdfname="hugene10stv1hsentrezgcdf"))
			#s.15 <- intersect(rownames(pheno)[which(pheno$Batch == 15)], s)
                        #b15 <- exprs(justRMA(filenames=s.15, celfile.path=cel.file.loc, cdfname="hugene10stv1hsentrezgcdf"))

			#xval.data <- cbind(b12, b13, b14, b15)
			xval.data <- cbind(b12, b13, b14)
			rm(b12)
			rm(b13)
			rm(b14)
			
		}
		
		
		gc()

		
	}

	## Batch Correction
	if(min.change <= 3) {
		if(var.batch.correct == "Combat") {
			batch = pheno[colnames(xval.data), "Batch"]
			mod<- rep(1, ncol(xval.data)) #include site into the null model
			xval.data.batch <- ComBat(dat=xval.data, batch=batch, mod=mod, numCovs=NULL, par.prior=T)  #re-assigns xval.data
		} else {
			xval.data.batch <- xval.data
		}

	}

	#Gene filter
	if(min.change <= 4) {
		gene.mad = apply(xval.data.batch, 1, mad)
		o = order(gene.mad, decreasing = T)
		num.to.select = round(var.gene.select*nrow(xval.data.batch))
		xval.gene.data <- xval.data.batch[o[1:num.to.select],]
	}


	## Feature Ranking
	if(min.change <= 6) {	
		
		## Set test/training dataset
		xval.train.data <- xval.gene.data[,xval.train.samples]
		xval.test.data <- xval.gene.data[,xval.test.samples]
		
		## Set training/test phenotypes
		xval.pheno <- pheno[colnames(xval.gene.data),]
		xval.train.pheno <- as.data.frame(xval.pheno[xval.train.samples,])
		xval.test.pheno <- as.data.frame(xval.pheno[xval.test.samples,])
			
		train.label <- as.numeric(as.factor(xval.train.pheno[,var.phenotype]))-1
		test.label <- as.numeric(as.factor(xval.test.pheno[,var.phenotype]))-1
		
		if(var.feat.select == "lm") {
			if(var.phenotype == "Sex") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Cancer + xval.train.pheno$Smk_Status + xval.train.pheno$Sex)
			}else if (var.phenotype == "Smk_Status") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Cancer + xval.train.pheno$Sex + xval.train.pheno$Smk_Status)
			}else if (var.phenotype == "Cancer") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Smk_Status + xval.train.pheno$Sex + xval.train.pheno$Cancer)
			}else if (var.phenotype == "random_Sex") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Cancer + xval.train.pheno$Smk_Status + xval.train.pheno$random_Sex)
			}else if (var.phenotype == "random_Smk_Status") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Cancer + xval.train.pheno$Sex + xval.train.pheno$random_Smk_Status)
			}else if (var.phenotype == "random_Cancer") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Smk_Status + xval.train.pheno$Sex + xval.train.pheno$random_Cancer)
			}else if (var.phenotype == "random_Control") {
				model <- model.matrix(~ xval.train.pheno$Batch + xval.train.pheno$Age + xval.train.pheno$Bronch_RIN + xval.train.pheno$Smk_Status + xval.train.pheno$Sex + xval.train.pheno$random_Control)
			}
			
			fit <- lmFit(xval.train.data, model)	
			fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
			
			o <- order(fit.t[,ncol(fit.t)])
			ranked.features <- rownames(xval.train.data)[o]	
		} else if (var.feat.select == "rf") {

                        rf = randomForest(t(xval.train.data), as.factor(train.label), ntree=1000, importance=TRUE)
                        o = order(rf$importance[,3], decreasing=T)
                        ranked.features <- rownames(xval.train.data)[o]
		}
	}

       ## Feature number
       if(min.change <= 7) {
	       if(var.num.feat > length(ranked.features)) {

		       select.features <- ranked.features	

	       } else if (var.feat.select == "lm") {
		      select.features <- unique(c(head(ranked.features, n=round(var.num.feat/2)), tail(ranked.features, n=round(var.num.feat/2))))
			
		} else if (var.feat.select == "rf") {

                        select.features <- head(ranked.features, n=var.num.feat)

		} else if (var.feat.select == "glmnet") {

			 model <- glmnet(t(xval.train.data), as.factor(train.label), family="binomial", dfmax=var.num.feat, alpha=0.5)

                        model.features <- model$beta[,ncol(model$beta)]
                        ind = which(model.features != 0)
                        select.features <- names(model.features)[ind]

		}


		xval.train.correct.select.data <- xval.train.data[select.features,]
		xval.test.correct.select.data <- xval.test.data[select.features,]
       }
	
	
       ## Model building/predictions
       if(min.change <= 8) {
		set.seed(seed)


		if(var.model.predict == "wv") {
		       predict.model <- wv.model(xval.train.correct.select.data, train.label)
	
		       predictions.train <- predict.wv(predict.model, xval.train.correct.select.data)
			predictions.train.scores <- predictions.train$scores
			predictions.train <- predictions.train$predictions
		       predictions.test <- predict.wv(predict.model, xval.test.correct.select.data)
			predictions.test.scores <- predictions.test$scores
                        predictions.test <- predictions.test$predictions
		       
	       } else if (var.model.predict == "knn") {
	
	       	       tune.obj <- tune.knn(t(xval.train.correct.select.data), as.factor(train.label), k=seq(3,15,by=2), l=0, tunecontrol = tune.control(sampling = "cross"), cross = 5)
		       	       	       
		       predictions.train <- knn(t(xval.train.correct.select.data), t(xval.train.correct.select.data), as.factor(train.label), k=tune.obj$best.model$k, l=tune.obj$best.model$l, prob=TRUE)
		       predictions.test <- knn(t(xval.train.correct.select.data), t(xval.test.correct.select.data), as.factor(train.label), k=tune.obj$best.model$k, l=tune.obj$best.model$l, prob=TRUE)
		       
			predictions.train.scores <- attr(predictions.train, "prob")
			predictions.test.scores <- attr(predictions.test, "prob")

		       predictions.train <- as.numeric(as.character(predictions.train))
	       	       predictions.test <- as.numeric(as.character(predictions.test))
			predictions.train.scores[predictions.train==0] <- 1 - predictions.train.scores[predictions.train==0]
			predictions.test.scores[predictions.test==0] <- 1 - predictions.test.scores[predictions.test==0]


		       	       
	       } else if (var.model.predict == "svm") { 
	       
			tune.obj <- tune.svm(t(xval.train.correct.select.data), as.factor(train.label), kernel="linear", cost = c(0.001,0.01,1,10), tunecontrol = tune.control(sampling = "cross"), cross = 5, type="C")

			predict.model <- svm(t(xval.train.correct.select.data), as.factor(train.label), kernel="linear", type="C", cost=tune.obj$best.model$cost)

			predictions.train <- predict(predict.model, t(xval.train.correct.select.data), decision.values=TRUE)
			predictions.test <- predict(predict.model, t(xval.test.correct.select.data), decision.values=TRUE)

			if(predict.model$labels[1] == 1) {
				predictions.train.scores <- -attr(predictions.train, "decision.values")
				predictions.test.scores <- -attr(predictions.test, "decision.values")
			} else {
				predictions.train.scores <- attr(predictions.train, "decision.values")
				predictions.test.scores <- attr(predictions.test, "decision.values")
			}

			predictions.train <- as.numeric(as.character(predictions.train))
			predictions.test <- as.numeric(as.character(predictions.test))
 
	       } else if (var.model.predict == "rf") {
	       
	       	       predict.model <- randomForest(t(xval.train.correct.select.data), as.factor(train.label), ntree=1000)
	       	       
	       	       predictions.train <- predict(predict.model, t(xval.train.correct.select.data))
	       	       predictions.test <- predict(predict.model, t(xval.test.correct.select.data))
	       	       
			predictions.train.scores <- predict(predict.model, t(xval.train.correct.select.data), type="prob")[,2]
			predictions.test.scores <- predict(predict.model, t(xval.test.correct.select.data), type="prob")[,2]

	       	       predictions.train <- as.numeric(as.character(predictions.train))
	       	       predictions.test <- as.numeric(as.character(predictions.test))
	       	       
	       } else if (var.model.predict == "nb") {
	       
	       	       predict.model <- naiveBayes(t(xval.train.correct.select.data), as.factor(train.label))
	       	       
	       	       predictions.train <- predict(predict.model, t(xval.train.correct.select.data))
	       	       predictions.test <- predict(predict.model, t(xval.test.correct.select.data))

			predictions.train.scores <- predict(predict.model, t(xval.train.correct.select.data), type="raw")[,2]
                       predictions.test.scores <- predict(predict.model, t(xval.test.correct.select.data), type="raw")[,2]

	       	       
	       	       predictions.train <- as.numeric(as.character(predictions.train))
	       	       predictions.test <- as.numeric(as.character(predictions.test))
	       	}else if (var.model.predict == "lda") {

                       predict.model <- lda(t(xval.train.correct.select.data), as.factor(train.label))

                       predictions.train <- predict(predict.model, t(xval.train.correct.select.data))
                       predictions.test <- predict(predict.model, t(xval.test.correct.select.data))

			predictions.train.scores <- predictions.train$x
                       	predictions.test.scores <- predictions.test$x


                       predictions.train <- as.numeric(as.character(predictions.train$class))
                       predictions.test <- as.numeric(as.character(predictions.test$class))
               	}else if (var.model.predict == "glmnet") {
			tune.obj <- cv.glmnet(t(xval.train.correct.select.data), as.factor(train.label), family="binomial", alpha=0.5)
			predict.model <- glmnet(t(xval.train.correct.select.data), as.factor(train.label), family="binomial", alpha=0.5)
			
			predictions.train <- as.numeric(predict(predict.model, t(xval.train.correct.select.data), s=tune.obj$lambda.min, type="class"))
			predictions.test <- as.numeric(predict(predict.model, t(xval.test.correct.select.data), s=tune.obj$lambda.min, type="class"))

			predictions.train.scores <- as.numeric(predict(predict.model, t(xval.train.correct.select.data), s=tune.obj$lambda.min, type="response"))
                        predictions.test.scores <- as.numeric(predict(predict.model, t(xval.test.correct.select.data), s=tune.obj$lambda.min, type="response"))

		}	

       }
	

	## Calculate prediction performance
	perform.train <- performance(predictions.train, train.label, scores=predictions.train.scores)
	perform.test <- performance(predictions.test, test.label, scores=predictions.test.scores)


	num1 <- sum(train.label == 0)
	num2 <- sum(train.label == 1)
	num3 <- sum(test.label == 0)
	num4 <- sum(test.label == 1)

	write.table(matrix(c(xval.iteration, param.index, perform.train, perform.test, num1, num2, num3, num4), nrow=1), file.out, sep="\t", quote=F, row.names=F, col.names=F, append=T)
	
	perform.train = NA
	perform.test = NA
	predict.train = NA
	predict.test = NA
	predictions.test = NA
	predictions.train = NA

	if(i %% 1000 == 0) {
		print(paste("Completed ", i, " (", round((i/nrow(index))*100,2), "%) iterations", sep=""))
	}
	
}


	
	
	
print(warnings())	





