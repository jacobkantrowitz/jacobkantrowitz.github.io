formatData <- function(eset, eset2, qcstats){
  
  # Making sure the genes are in the same order in both the bronchial and nasal expression sets
  eset2 <- eset2[row.names(exprs(eset)),]
  
  # Updating the gene names from probe IDs to gene symbols
  eset <- updateGeneNames(x=eset)
  eset2 <- updateGeneNames(x=eset2)
  
  # Creating a column to hold the CEL file name
  pData(eset2) <- cbind(pData(eset2),"CEL"=row.names(pData(eset2)))
  pData(eset) <- cbind(pData(eset),"CEL"=row.names(pData(eset)))
  
  # Create a column for RLE
  qc.stats.nasal <- read.table("/protected/projects/pulmarray/Biollegro/root/results/2015-04-01/150303_nasal_postQC_497_samples_QC_summary.txt",sep='\t',head=T,row.names=1)
  pData(eset) <- cbind(pData(eset),qc.stats.nasal[eset$CEL,c(1,2)])
  qc.stats.bronch <- read.table("/protected/projects/pulmarray/Biollegro/root/doc/300613_QCed_Bronc_RMA_QC_summary.txt",sep='\t',head=T,row.names=1)
  pData(eset2) <- cbind(pData(eset2),qc.stats.bronch[eset2$CEL,c(1,2)])
  
  # Setting the row.names of the phenotypic data and colnames of expression sets to the barcodes so that I can match up the samples for the correlation analysis
  row.names(pData(eset2)) <- as.character(eset2$BARCODE)
  row.names(pData(eset)) <- as.character(eset$BARCODE)
  colnames(exprs(eset2)) <- as.character(eset2$BARCODE)
  colnames(exprs(eset)) <- as.character(eset$BARCODE)
  
  return(list("eset"=eset,"eset2"=eset2))
  
}