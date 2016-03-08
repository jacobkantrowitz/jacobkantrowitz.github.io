create.table1 <- function(exprSet, groups, vars.to.summarize, tests.to.use=NULL,
                          save=F, filename=NULL, var.names=NULL, group.names=NULL,
                          binary.categories=NULL, round.ints=NULL){
  # exprSet is the eset object with fData with groups and vars.to.summarize
  # groups represent the groups by which data will be summarized
  # vars.to.summarize represent the variables that will be present in the table
  #              this can take 3 forms
  #              (1) a character vector of the colnames from the exprSet pData
  #              (2) a numeric vector indicating which columns from the pData
  #              (3) a logical vector indicating which columns from the pData
  # tests.to.use default to NULL; if not NULL, is a character vector of tests
  #              to use to test differences between groups for each variable;
  #              needs to be of length=length(vars.to.summarize)
  # save FALSE by default; if TRUE will output table 1 to file listed in filename
  # filename is the name of the file to output the table1
  # var.names is a vector of variable/column names with length=length(vars.to.summarize)
  # binary.categories is a dictionary of values for the binary categories
  #             e. g. list("GENDERc"=c("0"="Male", "1"="Female")); the names of the 
  #             items in the list should correspond to the vars.to.summarize
  # round.int is a vector of named integers values corresponding to at least
  #             one of the continuous variables and will change the number
  #             of significant digits included in the table e.g. c("RIN"=2)
  
  
  # require libraries...
  # TO DO - need to utilize the binary.categories to change the category levels
  #       - should fix the output of the categorical variables so they print clearly
  #       - need to implement saving; demogs and table1?
  #       - need to implement tests.to.use and add in the p-values results column(s)
  
  if(is.logical(vars.to.summarize)){
    if(length(vars.to.summarize)==dim(pData(target))[2]){
      vars.to.summarize <- varLabels(exprSet)[vars.to.summarize]
      print(paste("Variables pulled based on logicals:", vars.to.summarize))
    }
  } else if(is.numeric(vars.to.summarize)){
    if(all(vars.to.summarize < dim(pData(target))[2])){
      vars.to.summarize <- varLabels(exprSet)[vars.to.summarize]
      print(paste("Variables pulled based on numerics:", vars.to.summarize))
    }
  }
  
  # check whether groups and vars.to.summarize are all in pData(exprSet)
  if(!(all(groups %in% varLabels(exprSet)) & all(vars.to.summarize %in% varLabels(exprSet)))){
    print("ERROR: Not all groups and/or variables to summarize are present in expression set")
    return(NULL)
  }
  
  num.p.cols <- 1
  demogs <- data.frame(matrix(nrow=dim(exprSet)[2]))
  
  for(group in groups){
    demogs[[group]] <- exprSet[[group]]
  }
  
  for(variable in vars.to.summarize){
    demogs[[variable]] <- exprSet[[variable]]
  }
  
  demogs <- demogs[, -1]
  
  # determine the number of levels within each categorical variable so 
  # that each level can get its own row
  numrows <- 1
  for(x in 1:length(vars.to.summarize)){
    if(is.numeric(demogs[[vars.to.summarize[x]]])){
      numrows <- numrows + 1
    } else if(is.factor(demogs[[vars.to.summarize[x]]])){
      numrows <- numrows + length(levels(demogs[[vars.to.summarize[x]]]))
    }
  }
  
  if(!is.null(var.names)){
    if(length(var.names)==length(vars.to.summarize)){
      var.names <- var.names
      colnames(demogs) <- c(groups[1], var.names)
    } else {
      warning("var.names is given but not the correct number")
      return(NULL)
    }
  } else {
    var.names <- vars.to.summarize
  }
  
  
  # the table 1 will have as many columns as levels in groups and as many rows
  # as there are vars.to.summarize
  #table1 <- data.frame(matrix(nrow=length(vars.to.summarize)+1))
  table1 <- data.frame(matrix(nrow=numrows))
  for(i in 1:length(levels(demogs[[groups[1]]]))){
    row.ind <- 1
    group.2.test <- levels(demogs[[groups[1]]])[i]
    table1[row.ind, i] <- sum(demogs[[groups[1]]]==group.2.test)
    if(i==1){
      rownames(table1)[row.ind] <- "n"
    }
    row.ind <- row.ind + 1
    
    for(j in 1:length(var.names)){
      
      #group.2.test <- levels(demogs[[groups[1]]])[i]
      #var.2.test <- vars.to.summarize[j]
      var.2.test <- var.names[j]
      
      # set the round.int if a continous variable
      if(!is.null(round.ints)){
        round.int <- 0
        if(var.2.test %in% names(round.ints)){
          round.int <- round.ints[match(var.2.test, names(round.ints))]
        }
      } else{
        round.int <- 0
      }
      
      if(is.factor(demogs[[var.2.test]])){
        num.lvls <- length(levels(demogs[[var.2.test]]))
        var.summary <- summary(demogs[[var.2.test]][demogs[[groups[1]]]==group.2.test])
        table1[row.ind:(row.ind+num.lvls-1), i] <- var.summary
        if(i==1){
          rownames(table1)[row.ind:(row.ind+num.lvls-1)] <- paste(var.2.test, names(var.summary), sep=": ")
        }
        row.ind <- row.ind + num.lvls
        
      } else if(is.numeric(demogs[[var.2.test]])){
        table1[row.ind, i] <- continuous.summary(demogs[[var.2.test]][demogs[[groups[1]]]==group.2.test], round.int=round.int)
        if(i==1){
          rownames(table1)[row.ind] <- var.2.test
        }
        row.ind <- row.ind + 1
      }
    }
  }
  
  if(!is.null(group.names) & length(group.names)==length(levels(demogs[[groups[1]]]))){
    colnames(table1) <- group.names
  } else {
    colnames(table1) <- paste(groups[1], levels(demogs[[groups[1]]]), sep=": ")
  }
  
  
  #  if(!is.null(var.names) & length(var.names)==length(vars.to.summarize)){
  #    rownames(table1) <- c("n", var.names)
  #    } else {
  #      rownames(table1) <- c("n", vars.to.summarize)
  #      }
  
  if(save & !is.null(filename)){
    write.table(table1, quote=F, sep=",", file=paste("table1_", filename, ".csv", sep=""))
  }
  
  table1
  
}


# 0 is male
# 1 is current smoker



continuous.summary <- function(x, round.int=0, FUN2=sd){
  if(any(is.na(x))){
    warning("There are NAs in at least one continuous variable")
  }
  tempText <- round(mean(x, na.rm=T), digits=round.int)
  tempText <- paste(tempText, " (", round(FUN2(x, na.rm=T), digits=round.int), ")", sep="")
  
  tempText
}

binary.summary <- function(x){
  #tmpText <- paste(paste(levels(x), collapse=" "), "\n", sep="")
  tempText <- summary(x)[1]
  tempText <- paste(tempText, " (", summary(x)[2], ")", sep="")
  
  tempText
}

table.summary <- function(x, round.int=0, FUN2=sd){
  # apply the appropriate summary function based on the class of x
  if(is.factor(x) & length(levels(x))==2){
    #return(binary.summary(x))
    return(summary(x))
  } else if(is.numeric(x)){
    return(continuous.summary(x, round.int=round.int, FUN2=FUN2))
  }
}