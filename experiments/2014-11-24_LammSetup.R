# Generate and work with Lam (n=238) data for validation of COPD+Cancer signature (from Evans Day)
# 11-24-14

# CEL file location
celFileDir <- "/unprotected/projects/pulmarray/CEL/Lam_Airway_COPD_Samples/"
celFiles <- list.files(celFileDir)


# pheonotype file location
phenFile <- "/protected/projects/pulmarray/katie/lam_269egv16_log2rma/130930_237sampinfo.txt"
phenoData <- read.table(phenFile)
# phenoData$Type is the COPD indicator field (COPD n=83, normal=154, total=237)

# Affymetrix microarray normalization, annotation, and QC pipeline (Adam Gower)
naqc <- "/unprotected/projects/cbmhive/R/CBM_array_pipeline.R"


# need the pheno data column 1 to be the list of CEL file names
cels <- strsplit(celFiles, "_")
pullS <- function(d) { return(d[length(d)])}
cels[26] <- cels[[26]][[5]]

#pullS <- function(d) {
#  if(length(d)==2){
#    return(d[2])}
#  else{
#    return(d[length(d)])
#  }
#}
cels <- simplify2array(lapply(cels, pullS))
cels <- gsub("SL", "", cels)
cels <- gsub(".CEL", "", cels)
cels <- gsub("rehyb", "", cels)
cels[166] <- "273"
cels2 <- as.integer(cels)

phenoData$file <- celFiles[match(phenoData$s...27., cels2)]
phenoData$fileFull <- simplify2array(lapply(phenoData$file, function(d) {paste(celFileDir, d, sep="")}))

phenoData2 <- phenoData[, c(dim(phenoData)[2], 2:dim(phenoData)[2]-1 )]

write.table(phenoData2, file="lam_238_phenoData.txt", quote=FALSE)
