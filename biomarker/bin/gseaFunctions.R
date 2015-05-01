# Author: Joe Perez-Rogers
# Date: 2014-09-02
# Purpose: Functions to write gene sets and ranked lists in their appropriate formats

# writes a gene set
writeGeneSet <- function(genes,name="Default Gene Set Header",outfile="default.gmx"){
  v <- as.vector(c(name,"n/a",genes))
  m <- as.matrix(v,nrow=c(length(genes)+2))
  write.table(m,file=outfile,col.names=F,row.names=F,quote=F,sep="\t")
}

# writes a ranked list
writeRankedList <- function(genes,ranks,outfile="default.rnk"){
  m <- cbind(genes,ranks)
  write.table(m,file=outfile,col.names=F,row.names=F,quote=F,sep="\t")
}


# function to run gsea from within an R script and retrieve the gene symbols from leading edge
gseaPreranked <- function(gmx.file, rnk.file, output.label, output.dir,keep.files=F){

############################################################################
#Function to run GSEA from the command line and retreive gene symbols from the leading edge
#Written by: Joe Perez-Rogers, 10/2013
#Usage : gsea.out <- run.gsea(/path/to/my_file.gmx , /path/to/my_other_file.rnk , file_out_prefix , /path/to/output_directory)
#
#Example Files:
#gmx.file <- "/protected/projects/pulmarray/Biollegro/Combined/Analysis/Biomarker/Analysis/gsea/gmx/2013_10_15_Nasal_GSEA_fold_1_GSEA.gmx"
#rnk.file <- "/protected/projects/pulmarray/Biollegro/Combined/Analysis/Biomarker/Analysis/gsea/rnk/2013_10_15_Nasal_GSEA_fold_1_GSEA.rnk"
#output.label <- "Example_output"
#output.dir <- "/protected/projects/pulmarray/Biollegro/Combined/Analysis/Biomarker/Analysis/gsea"
#
#Current Limitations
#	>Only allows one gene set per run
#	>Must specify GSEA params within the function
############################################################################

	require(gdata)
	folders <- list.files(output.dir)
	#Check if gsea folder already exists. If it does, add v_# to the end of the folder name
	if(length(folders[which(startsWith(folders,output.label))])>0){
		file.exists <- TRUE
		i=2
		while(file.exists){
			n.output.label <- sprintf("%s_v%i",output.label,i)
			if(length(folders[which(startsWith(folders,n.output.label))])>0){
				i <- i + 1
			} else {
				file.exists <- FALSE
				output.label <- n.output.label
			}
		}
	}

	print(output.label)
	#Running GSEA via the command line
	system("module load java")
	run.cmd <- sprintf("java -cp /protected/projects/pulmarray/Biollegro/GSEA/gsea2-2.0.13.jar -Xmx1000m xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme weighted -rpt_label %s -include_only_symbols true -make_sets true -plot_top_x 20  -set_max 5000 -set_min 15 -zip_report false -out %s -gui false", gmx.file, rnk.file, output.label, output.dir)
	system(run.cmd)
	#Printing command that was used to run GSEA
	print(run.cmd)

	#Output directory
	folders <- list.files(output.dir)
	#Finding the GSEA folder for this run
	results.folder <- folders[which(startsWith(folders,output.label))]
	#Print that folder
	print(results.folder)
	#Full path to results file with leading edge genes
	gmx.infile <- read.table(gmx.file,head=F,sep="\t")
	fname <- as.character(gmx.infile[1,])
	results.file <- sprintf('%s/%s/%s',output.dir,results.folder,paste(toupper(fname),".xls",sep=""))
	#substr(basename(gmx.file),1,nchar(basename(gmx.file))-4)
	#Print that file path & name
	print(results.file)
	#Read in that file
	results <- read.table(results.file,sep="\t",head=T)
	# remove the gsea directory now that we've extracted the information we need
	if(!keep.files){
		system(sprintf('rm -rf %s/%s',output.dir,results.folder))
	}
	#Select gene symbols...labeled "PROBE" here that are enriched
	genes <- as.character(results$PROBE[which(results$CORE.ENRICHMENT=="Yes")])
	#Returning gene symbols list
	return(genes)
}

