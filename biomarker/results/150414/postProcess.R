# read in command line arguments
arguments <- commandArgs(TRUE)

# set fold count and file prefix
fold.count <- as.numeric(arguments[1])
prefix <- arguments[2]
#fold.count <- 2
#prefix <- "jjk_pipeline_test_posControlSmoking"
tDir <- arguments[3]

mat <- c()
for (i in 1:fold.count) {
	if(file.exists(paste(tDir, "/cv/",prefix,"_iter_",i,"_cv.txt",sep=""))){
		infile <- read.table(paste(tDir, "/cv/",prefix,"_iter_",i,"_cv.txt",sep=""),head=T,sep="\t")
		mat <- rbind(mat,infile)
	}
}

matrix.summary <- c()
for (mod in names(table(mat$Index))) {
	ss <- subset(mat,Index==mod)
	train.na <- sum(is.na(ss$TrainAUC))
	test.na <- sum(is.na(ss$TestAUC))
	ms <- c(ss[1,1:6],colMeans(na.omit(ss[,7:ncol(ss)])))
	mr <- c(as.character(ms$Index),as.character(ms$FF),as.character(ms$FS),as.character(ms$BS),as.character(ms$TO),as.character(ms$CL),as.character(ms$TrainAUC),as.character(ms$TestAUC),as.character(ms$TestHrAUC),as.character(ms$TestLrAUC),train.na,test.na)
	matrix.summary <- rbind(matrix.summary,mr)
}
colnames(matrix.summary) <- c(colnames(mat),"TrainNA","TestNA")

write.table(mat,file=paste(tDir, "/concat/",prefix,"_concat.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(matrix.summary,file=paste(tDir, "/summary/",prefix,"_summary.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)