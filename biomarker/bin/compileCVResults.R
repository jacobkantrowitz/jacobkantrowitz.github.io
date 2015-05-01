# Author: Joseph Perez-Rogers
# Date: 2014-05-25
# Purpose: This is a script to compile the results from several cross validation runs

compileCVResults <- function(prefix,max.iterations,write.global=FALSE){
	system("mkdir compiled_results")
	for(i in 1:max.iterations){
    print(i)
		file.name <- paste("cv/",prefix,"_",i,".txt",sep="")
		try(infile <- read.table(file.name, head=T, sep='\t'),TRUE)
		if(i==1){
			mat <- infile
		} else {
			mat <- rbind(mat,infile)
		}
	}
	if(write.global==TRUE){
		write.table(mat,paste("cv/",prefix,"_global_file.txt",sep=""),quote=F, row.names=F, col.names=T,append=F)
	}
	tmp.tag <- sample(1:1000,1)
	write.table(matrix(colnames(mat), nrow=1), paste("compiled_results/",prefix,"_compiled_results_TMP",tmp.tag,".txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=F)
	for(k in 1:length(table(mat$Index))){
		print(k)
		dat <- mat[which(mat$Index==k),]
		dat.subset <- dat[,8:(dim(mat)[2])]
		avg <- colMeans(dat.subset, na.rm=TRUE)
		write.table(c(dat[1,1:7],avg), paste("compiled_results/",prefix,"_compiled_results_TMP",tmp.tag,".txt",sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=T)
	}
}

