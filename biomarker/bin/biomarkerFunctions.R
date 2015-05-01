## Josh Campbell
## 4-13-11
## To use the functions in this script type: source("biomarker_functions.R")

## Load required libraries
library(gdata)

## Calculates all possible set of index permutations given the parameter lists
calculate.index <- function(params) {
  
  ## Calculate the number of parameters for each variable
  param.max = as.numeric(lapply(params, length))
  
  ## Total iterations
  total = prod(param.max)
  num.params = length(params)
  
  ## Set up a index matrix
  index <- matrix(0, nrow=total, ncol=num.params)
  
  for (i in 1:num.params) {
    if(i == 1) {
      e = 1	
    } 
    else {
      e = prod(param.max[1:(i-1)])
    }
    
    num = total/prod(param.max[1:i])
    
    index[,i] = rep(rep(1:param.max[i], each=e), num)
  }
  
  return(index)
}	

performance <- function(predictions, classlabel, scores=NULL) {
  # Calculates a variety of performance metrics for binary classifiers
  # Josh Campbell, 2011
  
  # INPUT
  # predictions: vector of predicted class labels (0's and 1's)
  # classlabel: vector of true class labels (0's and 1's)
  #
  # OUTPUT
  # A vector of performance metrics:
  #	Accuracy 				(TP + TN)/(TP + TN + FP + FN)
  #	Sensitivity				TP / (TP + FN)
  #	Specificity				TN / (FP + TN)
  #	Positive Predictive Value		TP / (TP + FP)
  #	Negative Predictive Value		TN / (TN + FN)
  #	Matthew's Correlation Coefficient	((TP x TN) - (FP x FN)) / sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN))
  #	Area Under the Curve			
  #	MAQCII metric				(0.5 * AUC) + (0.25 * (MCC+1))
  
  require(AUC)
  
  TP = sum(predictions == 1 & classlabel == 1, na.rm=T)
  TN = sum(predictions == 0 & classlabel == 0, na.rm=T)
  FP = sum(predictions == 1 & classlabel == 0, na.rm=T)
  FN = sum(predictions == 0 & classlabel == 1, na.rm=T)
  
  acc = (TP + TN)/(TP + TN + FP + FN)
  sens = TP / (TP + FN)
  spec = TN / (FP + TN)
  
  ppv = NA
  if(TP + FP != 0) {
    ppv = TP / (TP + FP)
  }
  
  npv = NA
  if(TN + FN != 0) {
    npv = TN / (TN + FN)
  }
  
  mcc.denominator <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
  if(is.na(mcc.denominator)){
    mcc = NA
  } else if(mcc.denominator == 0) {
    mcc = (TP * TN) - (FP * FN)
  } else {
    mcc = ((TP * TN) - (FP * FN)) / sqrt(mcc.denominator)
  }
  
  
  aucv = NA
  
  try({
    if(is.null(scores) & length(unique(predictions)) > 1){
      
      print(length(unique(predictions)))
      na <- is.na(classlabel) | is.na(predictions) | is.infinite(classlabel) | is.infinite(predictions)
      #roc = rocdemo.sca(classlabel[!na], predictions[!na], dxrule.sca)
      rocc <- AUC::roc(scores,as.factor(classlabel))
      #auc = AUC(roc)
      aucv <- AUC::auc(rocc)
      
    } else if (!is.null(scores) & length(unique(scores)) > 1){
      
      na <- is.na(classlabel) | is.na(scores) | is.infinite(classlabel) | is.infinite(scores)
      #roc = rocdemo.sca(classlabel[!na], predictions[!na], dxrule.sca)
      rocc <- AUC::roc(as.vector(scores),as.factor(classlabel))
      #auc = AUC(roc)
      aucv <- AUC::auc(rocc)
      
    } else {
      aucv = NA
    }
  }, TRUE)
  
  
  maqc2 = NA
  if(!is.na(aucv)) {
    maqc2 = (0.5 * aucv) + (0.25 * (mcc+1))
  }
  
  return(list(ACC=acc, SENS=sens, SPEC=spec, PPV=ppv, NPV=npv, MCC=mcc, AUC=aucv, MAQC2=maqc2))
}

combine.predictions <- function(v) {
  t1 <- sum(v == 0, na.rm=TRUE)
  t2 <- sum(v == 1, na.rm=TRUE)
  
  label <- NA
  
  if(t1 > t2) {
    label = 0
  } else if (t2 > t1) {
    label = 1
  }
  
  return(label)
}

wv.model <- function (data, classlabel, correction = TRUE) {
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
  #y <- lapply(levels(as.factor(classlabel)), function (i) {data[,classlabel==i,drop=F]});                           # y[[1]] = class 0 (generally controls); y[[2]] = class 1 (generally cases).
  mu <- sapply(y, apply, 1, mean, na.rm=T, simplify=F);                                   # mean of each probeset across all samples in each class
  sigma <- sapply(y, apply, 1, sd, na.rm=T, simplify=F);                                  # standard deviation of each probeset across all samples in each class
  if (correction) {sigma <- mapply(pmax, sigma, lapply(mu, `*`, 0.2),SIMPLIFY=F)};        # Correction from Broad Institute: sigma = max(sigma, 0.2*mu);
  #     i.e., %CV always >= 20% (communicated to me by Jen Beane)
  a <- (mu[[2]] - mu[[1]]) / (sigma[[1]] + sigma[[2]]);					# Compute the weights (signal to noise metric, similar to t statistic)
  g <- (mu[[1]] + mu[[2]]) / 2;								# Compute the means of each probeset across all samples
  return(list(weights=a, means=g));							# Return all computations in a named list
}

predict.wv <- function (model, data) {
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
  
  if(is.null(dim(V))){
    scores <- V[1]+V[2]
  } else {
    scores <- V[,1] + V[,2];                                                                # Calculate scores: sum of negative votes and positive votes
  }
  predictions <- as.numeric(scores > 0);
  if(is.null(dim(V))){
    strengths <- abs(scores) / (V[2] - V[1]);                                             # Prediction strengths: |scores| / sum |votes|
  } else {
    strengths <- abs(scores) / (V[,2] - V[,1]);                                             # Prediction strengths: |scores| / sum |votes|
  }
  
  return(list(scores=scores, predictions=predictions, strengths=strengths));           # Return all computations in a named list
}

pauc.p <- function(d, class.vector, p = 0.1, n = 10) {
  # d = data matrix with samples in columns
  # class.vector = class annotation for columns of d (see help for rowpAUCs function in genefilter for more info)
  # p = 1 - specificity value at which to determine pAUC (e.g. p = 0.1 is integral of sensitivity over the range of 90 - 100% specificity)
  # n = times to permute class labels and calculate pAUCs for each gene in dataset (used to estimate null distribution). So if the dataset has 20,000 genes, n = 10 will give you 200,000 pAUCs based on 10 different permutations of the true class label.
  require(genefilter)
  actual <- area(rowpAUCs(d, class.vector, p)) # this is how to calculate pAUC. This is the only useful thing in this function
  permuted <- NULL
  for (i in 1:n) {
    # this is a crude attempt to estimate the null distribution of pAUCs (globally rather than per gene)
    permuted <- c(permuted, area(rowpAUCs(d, sample(class.vector, length(class.vector)), p)))
    cat("*")
  }
  empiric.p <- function(val, random.values) {
    return(sum(random.values > val))
  }
  # the next step is not optimized for speed, but determines the empiric p-value for each observed pAUC based on the estimated null distribution
  p.val <- apply(cbind(actual), 1, function(x) empiric.p(x, permuted))
  p.val <- p.val / length(permuted)
  return(p.val)
}