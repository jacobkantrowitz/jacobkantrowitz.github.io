## Author: Joe Perez-Rogers
## Date: 2014-11-21
## Title: Biomarker Discovery Pipeline

# set paths
jjk.libs <- "/usr3/graduate/kantro/R/x86_64-unknown-linux-gnu-library/3.1/lars/libs"
cbm.libs <- "/restricted/projectnb/cbmhive/R_packages/R-3.0.0"
.libPaths(c(.libPaths(),file.path(cbm.libs),file.path(jjk.libs)));

# read in command line arguments
arguments <- commandArgs(TRUE)

# load libraries
library(affy)
library(limma)
library(mclust)
library(lars)
library(MASS)
library(ggplot2)
library(reshape2)
library(AUC)
library(densityClust)
library(gplots)
library(heatmap3)
library(class)
library(e1071)
library(randomForest)
library(glmnet)
library(sva)

# load user-defined functions
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/bin/sourceDirectory.R")
funcs <- sourceDir("/protected/projects/pulmarray/Allegro/COPD_Cancer/biomarker/bin/",trace=FALSE)

# set the seed and iteration
if(length(arguments)==0){
  iteration <- sample(1:1000,1)
  prefix <- "debug"
} else {
  iteration <- as.numeric(arguments[1])
  prefix <- arguments[2]
}
cat(paste("Iteration: ",iteration,"\n",sep=""))
new.seed <- setRandomSeed(x=iteration,print=TRUE)

# set output file prefix and name
outfile <- paste(prefix,"_iter_",iteration,"_cv.txt",sep="")

# load data
#nasal.all <- readRDS(file="/protected/projects/pulmarray/Biollegro/root/data/new_data/rma/nasal_rma_train.rds")
#bronch.all <- readRDS(file="/protected/projects/pulmarray/Biollegro/root/data/bronch_all_141121.rds")
####################################################
####  SOME CODE BY JAKE
####################################################
# load eset 
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/scripts/AllegroSetup.R")
#source("../2015-05-03_CBM_ATS_2015/plotGSEA.R")
# fix the one patient with wonky data
# eventually this should just be saved in the RDS file
holdEset <- eset
holdEset$FEV1Pc[holdEset$FEV1Pc==89.2] <- 0.892

esetClean <- holdEset
esetClean <- cleanNAForAnalysis(esetClean, "COPD2_R7")
esetClean <- cleanNAForAnalysis(esetClean, "RATIOc")
esetClean <- cleanNAForAnalysis(esetClean, "PYc")
esetClean <- cleanNAForAnalysis(esetClean, "RIN")
esetClean <- removeFactorLevel(esetClean, "COPD2_R7", "DK")
esetClean <- removeFactorLevel(esetClean, "FinalCaDXc", "DK")
esetClean <- removeFactorLevel(esetClean, "COPD2_R7", "DK")
esetClean <- removeFactorLevel(esetClean, "GENDERc", "DK")

esetClean <- calcIndicator(esetClean, "FinalCaDXc", "COPD2_R7")
esetClean$smkindic <- as.numeric(as.character(esetClean$indicator))
esetClean$smkindic[esetClean$SMKc==2] <- esetClean$smkindic[esetClean$SMKc==2] + 4
esetClean$smkindic <- as.factor(esetClean$smkindic)

copdESet <- removeFactorLevel(esetClean, "COPD2_R7", "0")
####################################################
#### END CODE BY JAKE
####################################################

# formatting data

# making sure genes are in the same order in bronch and nasal esets
#pp <- formatData(eset=nasal.all,eset2=bronch.all)
#nasal.all <- pp$eset
#bronch.all <- pp$eset2
#rm(pp)

# filtering genes based on low expression
# JK Comment - right now genderExpFilter doesn't work because
# the filter is based on selecting gene symbols and not affy ids
# this should be easy to fix
#noisy.genes.all <- genderExpressionFilter(x=nasal.all)
#nasal.all <- nasal.all[!(row.names(exprs(nasal.all))%in%noisy.genes.all),]
#bronch.all <- bronch.all[row.names(exprs(nasal.all)),]

# Subsetting the data to include only RIN>4
# nasal.all <- nasal.all[,nasal.all$RIN>4]
# bronch.all <- bronch.all[,bronch.all$RIN>4]

# split into training/testing sets
cv.index <- splitCV(copdESet,pct=0.20)
copdESet.train <- copdESet[,-cv.index]
copdESet.test <- copdESet[,cv.index]

# defining the class variable for training and test sets
#train.class <- class2Binary(copdESet.train$SMKc,case=1,control=2)
# JJK switched the code here to make the current smokers = 1 11/12/15
#train.class <- as.factor(as.numeric(copdESet.train$SMKc)*-1 + 2)
#test.class <- class2Binary(copdESet.test$SMKc,case=1,control=2)
#test.class <- as.factor(as.numeric(copdESet.test$SMKc)*-1 + 2)
train.class <- copdESet.train$FinalCaDXc
test.class <- copdESet.test$FinalCaDXc


#creating the parameter list of all models to be tested
parameter.list <- list(
  # ff=c("cor","var","gsea200","none"), #"gsea500"
  ff=c("var","none"),
  fs=c("lm-t","lm-abst","pAUC","glmnet","rf"),
  bs=c(5,25,50,100),#seq(from=5,to=100,by=5)),
  to=c("all"),#"highrin","all"),
  cl=c("wv","nb","glmnet","knn","svm")#"glmnet","wv","nb","knn","svm","rf")
)

index <- calculate.index(parameter.list)
index <- index[order(index[,1], index[,2], index[,3], index[,4]),]

# starting a new results file
need.header <- TRUE

# run through all possible models
cat("Running through all possible models\n")

for(i in 1:nrow(index)) {
  # define model parameters
  ff.m <- parameter.list$ff[index[i,1]]
  fs.m <- parameter.list$fs[index[i,2]]
  bs.m <- parameter.list$bs[index[i,3]]
  to.m <- parameter.list$to[index[i,4]]
  cl.m <- parameter.list$cl[index[i,5]]
  
  # create a list to store all important information about a given iteration
  bmod <- list()
  bmod[["index"]] <- i
  bmod[["ff"]] <- ff.m
  bmod[["fs"]] <- fs.m
  bmod[["bs"]] <- bs.m
  bmod[["to"]] <- to.m
  bmod[["cl"]] <- cl.m
  
  # find the minimum changing parameter
  min.change <- 1
  if(i > 1) {
    parameter.previous <- index[i-1,]
    parameter.index <- index[i,]
    min.change <- min(which(parameter.index != parameter.previous))
  }
  min.change
  
  print(i)
  
  ###################################################
  # Feature Filtering (Unsupervised)                #
  ###################################################
  
  cat(paste("Feature filtering...",ff.m,"\n",sep=""))
  
  if(min.change <=1 ) {
    
    # none
    if(ff.m=="none"){
      ff.eset <- copdESet.train
      ff.genes <- row.names(exprs(ff.eset))
    }
    
    # variance
    if(ff.m=="var"){
      var.p <- 0.05
      v <- apply(exprs(copdESet.train),1,var)
      med.idx <- whichMedian(v)
      p <- apply(exprs(copdESet.train),1,function(x){
        var.test(x,exprs(copdESet.train)[med.idx,],alternative="greater")$p.value
      })
      ff.eset <- copdESet.train[p<var.p,]
      ff.genes <- row.names(exprs(ff.eset))
    }
        
  }
  
  ###################################################
  # Feature Selection (Supervised)                  #
  ###################################################
  
  cat(paste("Feature selection...",fs.m,"\n",sep=""))
  
  if(min.change <=2 ) {
    
    if(fs.m=="lm-t"){
      # linear model sorting genes by t-statistic
      model <- model.matrix(~1+train.class, data=pData(ff.eset),na.action=na.exclude)
      fit <- lmFit(exprs(ff.eset),model,na.action=na.exclude)  
      fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
      fs.eset <- ff.eset[names(sort(fit.t[,2],decreasing=T)),]
      fs.genes <- row.names(exprs(fs.eset))
      
    } else if(fs.m=="lm-abst"){
      # linear model sorting genes by absolute t-statistic
      model <- model.matrix(~1+train.class, data=pData(ff.eset),na.action=na.exclude)
      fit <- lmFit(exprs(ff.eset),model,na.action=na.exclude)  
      fit.t <- fit$coef / fit$stdev.unscaled / fit$sigma
      fs.eset <- ff.eset[names(sort(abs(fit.t[,2]),decreasing=T)),]
      fs.genes <- row.names(exprs(fs.eset))
      
    } else if(fs.m=="pAUC"){
      # computing the partial AUC of genes with 90% sensitivity
      p.out <- pauc.p(exprs(ff.eset),train.class)
      fs.eset <- ff.eset[names(sort(p.out,decreasing=F)),]
      fs.genes <- row.names(exprs(fs.eset))
      
    } else if(fs.m=="rf"){
      require(randomForest)
      # building a random forest with 5000 trees and ranking genes based on their importance score
      rf.out <- randomForest(t(exprs(ff.eset)), as.factor(train.class), ntree=5000, importance=TRUE)
      fs.eset <- ff.eset[names(sort(rf.out$importance[,3], decreasing=T)),]
      fs.genes <- row.names(exprs(fs.eset))
      
    } else if(fs.m=="glmnet"){
      # computing the lasso with elastic net (glmnet) to rank features
      model <- glmnet(t(as.matrix(exprs(ff.eset))),as.factor(train.class),family="binomial",dfmax=100,alpha=0.5)
      beta <- model$beta[,ncol(model$beta)]
      fs.eset <- ff.eset[names(sort(beta,decreasing=T)),]
      fs.genes <- row.names(exprs(fs.eset))
      
    } 
  }
  
  ###################################################
  # Select Biomarker Size                           #
  ###################################################
  
  cat(paste("Biomarker size...",bs.m,"\n",sep=""))
  
  if(min.change <= 3) {
    
    if(fs.m=="lm-t"){
      #picking the most up and most down in a balanced way - for uneven biomarker size (e.g. 5)
      # results in biomarker size = bs.m - 1
      bs.genes <- c(head(row.names(exprs(fs.eset)),round(bs.m/2)),tail(row.names(exprs(fs.eset)),round(bs.m/2)))
      bs.eset <- fs.eset[bs.genes,]
    } else {
      bs.genes <- head(row.names(exprs(fs.eset)),bs.m)
      bs.eset <- fs.eset[bs.genes,]
    }
    
  }
  
  ###################################################
  # Select Training Set For Classifier              #
  ###################################################
  
  cat(paste("Training model on...",bs.m,"\n",sep=""))
  
  if(min.change <= 4) {
    if (to.m=="highrin") {
      to.eset <- bs.eset[,bs.eset$RIN>4]
      to.class <- train.class[bs.eset$RIN>4]
    } else if (to.m=="all") {
      to.eset <- bs.eset
      to.class <- train.class
    }
  }
  
  ###################################################
  # Build Classifiers (Supervised)                  #
  ###################################################
  
  cat(paste("Building classifier...",cl.m,"\n",sep=""))
  
  if(min.change <= 5) {
    
    # clean up the workspace
    perform.train = NA
    perform.test = NA
    predict.train = NA
    predict.test = NA
    predictions.test = NA
    predictions.train = NA
    
    # define the test set
    test.eset <- copdESet.test[row.names(exprs(to.eset)),]
    
    if(cl.m== "wv") {
      
      predict.model <- wv.model(exprs(to.eset),to.class)
      
      predictions.train <- predict.wv(predict.model,exprs(to.eset))
      predictions.train.scores <- predictions.train$scores
      predictions.train <- predictions.train$predictions
      
      predictions.test <- predict.wv(predict.model,exprs(test.eset))
      predictions.test.scores <- predictions.test$scores
      predictions.test <- predictions.test$predictions
      
    } else if (cl.m == "knn") {
      
      f.levels <- levels(as.factor(to.class))
      
      tune.obj <- tune.knn(t(exprs(to.eset)),to.class,k=seq(3,15,by=2), l=1, tunecontrol = tune.control(sampling = "cross"), cross = 5)
      
      predictions.train <- knn(t(exprs(to.eset)), t(exprs(to.eset)), as.factor(to.class), k=tune.obj$best.model$k, l=tune.obj$best.model$l, prob=TRUE)
      predictions.train.scores <- attr(predictions.train, "prob")
      predictions.train <- as.numeric(as.character(predictions.train))                         
      predictions.train.scores[predictions.train==f.levels[1]] <- 1 - predictions.train.scores[predictions.train==f.levels[1]]                         
      
      predictions.test <- knn(t(exprs(to.eset)), t(exprs(test.eset)), as.factor(to.class), k=tune.obj$best.model$k, l=tune.obj$best.model$l, prob=TRUE)
      predictions.test.scores <- attr(predictions.test, "prob")
      predictions.test <- as.numeric(as.character(predictions.test))
      predictions.test.scores[predictions.test==f.levels[1]] <- 1 - predictions.test.scores[predictions.test==f.levels[1]]
      
    } else if (cl.m== "svm") { 
      
      # logistics
      f.levels <- levels(as.factor(to.class))
      class.weights<-c((sum(to.class==f.levels[1])/length(to.class)),(sum(to.class==f.levels[2])/length(to.class)))
      names(class.weights) <- f.levels
      
      # building model
      tune.obj <- tune.svm(t(exprs(to.eset)), to.class, kernel="linear", cost = c(0.001,0.01,1,10,100), tunecontrol = tune.control(sampling = "cross",cross=5), cross=5,type="C",class.weights=class.weights)
      predict.model <- svm(t(exprs(to.eset)), to.class, kernel="linear", type="C", cost=tune.obj$best.model$cost,,class.weights=class.weights)
      
      # training set
      predictions.train <- predict(predict.model, t(exprs(to.eset)), decision.values=TRUE)
      if(predict.model$labels[1] == 1) {
        predictions.train.scores <- -attr(predictions.train, "decision.values")
      } else {
        predictions.train.scores <- attr(predictions.train, "decision.values")
      }
      predictions.train <- as.numeric(as.character(predictions.train))
      
      # test set
      predictions.test <- predict(predict.model, t(exprs(test.eset)), decision.values=TRUE)
      if(predict.model$labels[1] == 1) {
        predictions.test.scores <- -attr(predictions.test, "decision.values")
      } else {
        predictions.test.scores <- attr(predictions.test, "decision.values")
      }
      predictions.test <- as.numeric(as.character(predictions.test))
      
    } else if (cl.m == "rf") {
      
      predict.model <- randomForest(t(exprs(to.eset)), to.class, ntree=5000)
      
      predictions.train <- predict(predict.model, t(exprs(to.eset)))
      predictions.train.scores <- predict(predict.model, t(exprs(to.eset)), type="prob")[,2]
      predictions.train <- as.numeric(as.character(predictions.train))
      
      predictions.test <- predict(predict.model, t(exprs(test.eset)))
      predictions.test.scores <- predict(predict.model, t(exprs(test.eset)), type="prob")[,2]
      predictions.test <- as.numeric(as.character(predictions.test))
      
    } else if (cl.m == "nb") {
      
      predict.model <- naiveBayes(t(exprs(to.eset)), to.class)
      
      predictions.train <- predict(predict.model, t(exprs(to.eset)))
      predictions.train.scores <- predict(predict.model, t(exprs(to.eset)), type="raw")[,2]
      predictions.train <- as.numeric(as.character(predictions.train))
      
      predictions.test <- predict(predict.model, t(exprs(test.eset)))
      predictions.test.scores <- predict(predict.model, t(exprs(test.eset)), type="raw")[,2]
      predictions.test <- as.numeric(as.character(predictions.test))
      
    }else if (cl.m == "lda") {
      
      predict.model <- lda(t(exprs(to.eset)), to.class)
      
      predictions.train <- predict(predict.model, t(exprs(to.eset)))
      predictions.train.scores <- predictions.train$x
      predictions.train <- as.numeric(as.character(predictions.train$class))
      
      predictions.test <- predict(predict.model, t(exprs(test.eset)))
      predictions.test.scores <- predictions.test$x
      predictions.test <- as.numeric(as.character(predictions.test$class))
      
    }else if (cl.m == "glmnet") {
      
      tune.obj <- cv.glmnet(t(exprs(to.eset)), to.class, family="binomial", alpha=0.5)
      predict.model <- glmnet(t(exprs(to.eset)), to.class, family="binomial", alpha=0.5)
      
      predictions.train <- as.numeric(predict(predict.model, t(exprs(to.eset)), s=tune.obj$lambda.1se, type="class"))
      predictions.train.scores <- as.numeric(predict(predict.model, t(exprs(to.eset)), s=tune.obj$lambda.1se, type="response"))
      
      predictions.test <- as.numeric(predict(predict.model, t(exprs(test.eset)), s=tune.obj$lambda.1se, type="class"))
      predictions.test.scores <- as.numeric(predict(predict.model, t(exprs(test.eset)), s=tune.obj$lambda.1se, type="response"))      
      
    }
    
  }
  
  ###################################################
  # Summarization                                   #
  ###################################################
  
  # evaluate the performance of the model
  perform.train <- performance(predictions.train, as.factor(to.class), scores=predictions.train.scores)
  perform.test <- performance(predictions.test, as.factor(test.class), scores=predictions.test.scores)
  perform.test.hr <- performance(predictions.test[test.eset$RIN>=4], as.factor(test.class)[test.eset$RIN>=4], scores=predictions.test.scores[test.eset$RIN>=4])
  perform.test.lr <- performance(predictions.test[test.eset$RIN<4], as.factor(test.class)[test.eset$RIN<4], scores=predictions.test.scores[test.eset$RIN<4])
  
  # adding prediction information to the output lists
  perform.train[["predictions"]] <- predictions.train
  perform.train[["scores"]] <- predictions.train.scores
  perform.test[["predictions"]] <- predictions.test
  perform.test[["scores"]] <- predictions.test.scores
  perform.test.hr[["predictions"]] <- predictions.test[test.eset$RIN>=4]
  perform.test.hr[["scores"]] <- predictions.test.scores[test.eset$RIN>=4]
  perform.test.lr[["predictions"]] <- predictions.test[test.eset$RIN<4]
  perform.test.lr[["scores"]] <- predictions.test.scores[test.eset$RIN<4]
  
  # adding prediction information and performance metrics to the bmod object
  bmod[["train"]] <- perform.train
  bmod[["test"]] <- perform.test
  bmod[["test.hr"]] <- perform.test.hr
  bmod[["test.lr"]] <- perform.test.lr
  
  # outputting the AUC for debugging purposes
  cat(paste("Training Set AUC: ",perform.train$AUC,"\n",sep=""))
  cat(paste("Test Set AUC: ",perform.test$AUC,"\n",sep=""))
  cat(paste("Test Set HR AUC: ",perform.test.hr$AUC,"\n",sep=""))
  cat(paste("Test Set LR AUC: ",perform.test.lr$AUC,"\n",sep=""))
  
  # writing the results of the iteration to a file
  if(need.header){
    write.table(matrix(c("Index","FF","FS","BS","TO","CL","TrainAUC","TestAUC","TestHrAUC","TestLrAUC"),nrow=1),file=outfile,sep="\t", quote=F, row.names=F, col.names=F, append=F)
    need.header <- FALSE
  }
  
  output.matrix <- matrix(c(bmod$index,bmod$ff,bmod$fs,bmod$bs,bmod$to,bmod$cl,bmod$train$AUC,bmod$test$AUC,bmod$test.hr$AUC,bmod$test.lr$AUC),nrow=1)
  write.table(output.matrix, file=outfile, sep="\t", quote=F, row.names=F, col.names=F, append=T)
  
  # put some space between model outputs
  cat("\n\n")
  
} # end cross-validation loop

###################################################
# Output Files                                    #
############s#######################################

# writing out the genes selected in the feature selection module
## re-reading in the nasal.all expression matrix so that we can get the full list of genes
#nasal.all <- readRDS(file="/protected/projects/pulmarray/Biollegro/root/data/train_nasal_141121.rds")
#bronch.all <- readRDS(file="/protected/projects/pulmarray/Biollegro/root/data/bronch_all_141121.rds")
#pp <- formatData(eset=nasal.all,eset2=bronch.all)
#nasal.all <- pp$eset
## creating an empty matrix of 0's that we will add +1's to for each selection
fs.genes.table <- matrix(0,nrow=nrow(exprs(holdEset)))
row.names(fs.genes.table) <- row.names(exprs(holdEset))
fs.genes.table[ff.genes,] <- 1
write.table(fs.genes.table,file=paste(prefix,"_iter_",iteration,"_fs_genes.txt",sep=""))

# writing out the genes selected in the biomarker




















