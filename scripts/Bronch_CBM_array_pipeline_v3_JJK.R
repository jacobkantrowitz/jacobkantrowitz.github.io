# Microarray pipeline for creating an ExpressionSet object of RMA-normalized data from Affymetrix CEL files
# Also performs QC (RLE, NUSE, PCA)
# Adam Gower
# Version 2013-06-25

# This script performs the following tasks:
#   - Reads a file of CEL filenames and sample annotation
#   - RMA-normalizes the CEL files to produce log2-transformed expression levels
#   - Automatically reads the gene symbols and descriptions for that platform from an R package
#   - Combines all of this information into an ExpressionSet object
#   - Writes the ExpressionSet to a compressed RDS file
#   - Produces QC plots as specified and stores QC metrics for future reference

# To use this script:
#   - Set any variables at the topmost section as needed
#   - Create a tab-delimited or CSV text file in which the first column contains CEL file names
#     and the remaining columns contain sample-specific annotation (age, sex, smoking status, etc.)
#   - Run the script from the command line, e.g.,
#     R-2.15.1 --vanilla < my_array_pipeline.R >& my_array_pipeline_output.log
#     NOTE: this example command will save the standard output and standard error to the file 'my_array_pipeline_output.log'.

# The output of this script includes the following files:
#   - (optional) JPGs of chipwide residual plots, from the RMA probe-level model (PLM) fit during computation of RLE & NUSE
#   - (optional) A PDF with boxplots of RLE and NUSE values, with color-coding to indicate samples that are outliers
#   - A PDF with histograms of the RMA-normalized gene expression in each sample, one page per sample
#   - A PDF with a 2D Principal Component Analysis (PCA) plot, showing the percent of variance explained by PC1 and PC2
#   - An ExpressionSet object containing the RMA-normalized gene expression data, stored in the compressed binary RDS file format
#   - An RDS file containing a matrix of RLE values for each feature and sample
#   - An RDS file containing a matrix of NUSE values for each feature and sample
# Note: RDS files can be read using the readRDS() command in R versions >= 2.13.0.

########## Set file paths ##########

# Specify location of shared packages
shared.library.path <- file.path("/unprotected/projects/cbmhive/R_packages", sprintf("R-%s", getRversion()));

# Specify location of Brainarray master package folder
brainarray.path <- "/unprotected/projects/pulmarray/BrainArray";

########## Set run-specific variables ##########

# Do you want to produce JPGs of residual heatmaps?  (Useful for finding chip defects like bubbles)
draw.residual.heatmaps <- FALSE;
# Do you want to produce a PDF of boxplots of RLE and NUSE values?  (Useful for finding outlier samples)
draw.boxplots <- FALSE;

# Change the output path and prefix as needed
output.path <- "/protected/projects/pulmarray/Biollegro/RMA/Temp"# e.g., "/protected/projects/pulmarray/My_project"
output.prefix <- "300613_QCed_Bronc_RMA"# e.g., "My_project_yymmdd"

# Change Brainarray mapping and version as needed
brainarray.mapping <- "entrezg"; # could also be "ensg", "refseq", etc.
brainarray.version <- "17.0.0" # e.g., "14.0.0"

# Change the CEL file path as needed
celfile.path <- "/protected/projects/pulmarray/Allegro/Bronch_mRNA/Data/CEL_12_13_14_15_16" # e.g., "/unprotected/projects/pulmarray/My_project/CEL"

# Supply the name of an annotation file (tab-delimited or CSV) that has the CEL file names in the first column.
# This will be used to populate the phenotypic data frame in the ExpressionSet object.
sample.annotation.filename <- "/protected/projects/pulmarray/Biollegro/RMA/scripts/Bronch_Annotation.txt"# e.g., "protected/projects/pulmarray/My_project/My_project_annotation.txt"

########## Other global variables and functions ##########

# Cutoff values for declaring a sample out of bounds by quality metrics
QC.cutoffs <- list(RLE = 0.1, NUSE = 1.05);

# The saveRDS() function was not implemented prior to R-2.13.0; otherwise, this is the function definition (from R-2.15.1)
if ((R.version$major < 2) || ((R.version$major == 2) && (compareVersion(R.version$minor, "13.0") < 0))) {
	saveRDS <- function (object, file = "", ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL) {
		if (is.character(file)) {
			if (file == "") stop("'file' must be non-empty string")
			mode <- if (ascii) "w" else "wb"
			con <-
				if (identical(compress, "bzip2")) bzfile(file, mode)
				else if (identical(compress, "xz")) xzfile(file, mode)
				else if (compress) gzfile(file, mode)
				else file(file, mode)
			on.exit(close(con))
		}
		else if (inherits(file, "connection")) {
			if (!missing(compress)) warning("'compress' is ignored unless 'file' is a file name")
			con <- file
		}
		else stop("bad 'file' argument")
		invisible(.Internal(serializeToConn(object, con, ascii, version, refhook)))
	}
}

########## Check input variables for validity, read annotation table, and determine CDF name ##########

# Check 'draw.residual.heatmaps' variable
if (!(identical(draw.residual.heatmaps, TRUE) || identical(draw.residual.heatmaps, FALSE))) {
	stop(sprintf("The 'draw.residual.heatmaps' variable must be either TRUE or FALSE"));
}
# Check 'draw.boxplots' variable
if (!(identical(draw.boxplots, TRUE) || identical(draw.boxplots, FALSE))) {
	stop(sprintf("The 'draw.boxplots' variable must be either TRUE or FALSE"));
}

# Check 'output.path' variable
if ((length(output.path) != 1) || !is.character(output.path)) {
	stop(sprintf("The 'output.path' variable must be a character vector of length one"));
}
if (!file.exists(output.path)) {
	stop(sprintf("The output path '%s' does not exist", output.path));
}
if (!file.info(output.path)$isdir) {
	stop(sprintf("The output path '%s' is not a directory", output.path));
}

# Check 'output.prefix' variable
if ((length(output.prefix) != 1) || !is.character(output.prefix)) {
	stop(sprintf("The 'output.prefix' variable must be a character vector of length one"));
}

# Check 'brainarray.version' variable
if (!file.exists(file.path(brainarray.path, brainarray.version))) {
	stop(sprintf("Version '%s' of Brainarray is not available", brainarray.version))
}

# Check 'celfile.path' variable
if ((length(celfile.path) != 1) || !is.character(celfile.path)) {
	stop(sprintf("The 'celfile.path' variable must be a character vector of length one"));
}
if (!file.exists(celfile.path)) {
	stop(sprintf("CEL file path '%s' does not exist", celfile.path));
}
if (!file.info(celfile.path)$isdir) {
	stop(sprintf("CEL file path '%s' is not a directory", celfile.path));
}

# Check 'sample.annotation.filename' variable
if ((length(sample.annotation.filename) != 1) || !is.character(sample.annotation.filename)) {
	stop(sprintf("The sample annotation filename must be a character vector of length one"));
}
if (!file.exists(sample.annotation.filename)) {
	stop(sprintf("Sample annotation file '%s' does not exist", sample.annotation.filename));
}

########## Set library paths and load packages ##########

# Add Brainarray library and shared library to the front of the search path
.libPaths(c(file.path(brainarray.path, brainarray.version), shared.library.path, .libPaths()));

# Load packages required for normalization, QC, and reading CEL files
library(affy);
library(affyio);
library(affyPLM);

# Open the annotation file
extension <- tolower(sub("\\.(.+)$", "\\1", sample.annotation.filename));
if (extension == "csv") {
	# If the file ends with .CSV or .csv, assume that a CSV file is supplied
	sample.annotation <- read.csv(sample.annotation.filename, row.names=1, check.names=FALSE, stringsAsFactors=FALSE);
} else {
	# Otherwise, assume that a tab-delimited text file is supplied
	sample.annotation <- read.delim(sample.annotation.filename, row.names=1, check.names=FALSE, stringsAsFactors=FALSE);
}

# Identify the array platform(s) of the CEL filenames provided and terminate if more than one platform is represented
cel.filenames <- rownames(sample.annotation);
if (all(file.exists(file.path(celfile.path, cel.filenames)))) {
	cel.platform <- unique(unlist(sapply(file.path(celfile.path, cel.filenames), read.celfile.header, info="full")["cdfName", ]));
	if (length(cel.platform) > 1) {
		stop("All CEL files must belong to the same array platform");
	} else {
		cat(sprintf("Automatically identified CEL file platform as: '%s'\n", cel.platform));
	}
	# Identify the CDF using grep of BrainArray package names.  Note the following:
	# - ".." is used in place of the two-character species code
	# - CDF names for Gene ST arrays may omit the "v1" for some Brainarray versions, so this is made optional in the regex
	cdf.regex <- sprintf("^%s..%scdf$", tolower(gsub("[-_\\.]", "", cel.platform)), brainarray.mapping);
	if (grepl("gene1[01]st", cdf.regex)) cdf.regex <- sub("(gene1[01]st)v1", "\\1(v1)?", cdf.regex);
	project.cdf <- grep(cdf.regex, list.files(file.path(brainarray.path, brainarray.version)), value=TRUE);
	if (length(project.cdf) == 0) {
		stop(sprintf(
			"The character string '%s' does not correspond to a known mapping for this platform and Brainarray version",
			brainarray.mapping
		));
	}
} else {
	stop(sprintf(
		"The following CEL files do not exist within '%s': '%s'",
		celfile.path, paste(cel.filenames[!file.exists(cel.filenames)], collapse="','"
	)));
}

print("passed")
print(stophere)

# Load the package corresponding to the CDF
library(project.cdf, character.only=TRUE);
# Load the corresponding feature-specific annotation database package corresponding to the CDF
project.db <- sub("cdf", ".db", project.cdf);
library(project.db, character.only=TRUE);

# Show information about packages that have been loaded, etc.
sessionInfo();

########## Create and save an ExpressionSet object ##########

# Create an AffyBatch using the CDF corresponding to the specified BrainArray mapping and version
cat(sprintf("Creating an AffyBatch object from %d CEL files.\n", length(cel.filenames)));
abatch <- ReadAffy(filenames=cel.filenames, celfile.path=celfile.path, cdfname=project.cdf);
# RMA-normalize the data to generate an ExpressionSet
cat("Performing RMA normalization.\n");
dataset <- rma(abatch);

# Create a data frame to hold the phenotypic (sample-specific) annotation
pData(dataset) <- sample.annotation;
# Extract the name of the CEL file from the CEL filenames supplied
sampleNames(dataset) <- basename(sampleNames(dataset));

# Create an empty data frame in the featureData slot of the ExpressionSet object to hold feature-specific annotation
cat("Retrieving feature-specific annotation.\n");
fData(dataset) <- data.frame(row.names = rownames(exprs(dataset)));
# Add feature-specific annotation for each database specified
databases <- c("Symbol"="SYMBOL", "Description"="GENENAME");
for (i in 1:length(databases)) {
	database.object <- eval(as.name(sub(".db", databases[i], project.db)));
	database.keys <- keys(database.object);
	fData(dataset)[[names(databases)[i]]] <- NA;
	fData(dataset)[database.keys, names(databases)[i]] <- unlist(mget(database.keys, database.object));
}

# Write the ExpressionSet to an .rds file
output.filename <- file.path(output.path, sprintf("%s_ExpressionSet.rds", output.prefix));
cat(sprintf("Writing ExpressionSet to '%s'.\n", output.filename));
saveRDS(dataset, file=output.filename);

########## Quality control ##########

# Fit a PLMset object to the AffyBatch
cat("Fitting a probe-level model to the AffyBatch.\n");
pset <- fitPLM(abatch);

# If requested, create images of chipwide residuals and write to files
if (draw.residual.heatmaps) {
	for (i in 1:nrow(pData(dataset))) {
		output.filename <- file.path(output.path, sprintf("%s_residuals.jpg", sub("^(.+)\\.(CEL|cel)(\\.gz)?$", "\\1", sampleNames(dataset)[i])));
		cat(sprintf("Drawing chipwide residual plot for '%s' to '%s'.\n", sampleNames(dataset)[i], output.filename));
		jpeg(filename=output.filename);
		image(pset, which=i, type="resids", col=pseudoPalette(low="blue", mid="white", high="red"), add.legend=TRUE);
		dev.off();
	}
}

cat("Computing RLE and NUSE metrics.\n");
QC <- list();
QC$RLE <- RLE(pset, type="values");
QC$NUSE <- NUSE(pset, type="values");

QC.medians <- list();
for (metric in names(QC)) {
	QC.medians[[metric]] <- apply(QC[[metric]], 2, median);
}

# Save matrices of QC metrics as RDS files
for (metric in names(QC)) {
	output.filename <- file.path(output.path, sprintf("%s_%s.rds", output.prefix, metric));
	saveRDS(QC[[metric]], output.filename);
}

# If requested, create a PDF of boxplots of the RLE and NUSE values
if (draw.boxplots) {
	output.filename <- file.path(output.path, sprintf("%s_RLE_NUSE.pdf", output.prefix));
	cat(sprintf("Drawing RLE and NUSE boxplots to '%s'.\n", output.filename));
	pdf(output.filename, width=11, height=8.5);
	boxplot(
		QC$RLE,
		main="RLE (Relative Log Expression)\nShould be centered on 0 (blue line)\nRed sample = out of bounds (dashed red line)",
		names=sampleNames(dataset), las=2,
		border=c("black","red")[(QC.medians$RLE > QC.cutoffs$RLE)+1]	
	);
	lines(x=c(0,ncol(QC$RLE)+1), y=rep(0,2), col="blue", lty=2);
	lines(x=c(0,ncol(QC$RLE)+1), y=rep(QC.cutoffs$RLE,2), col="red", lty=2);
	boxplot(
		QC$NUSE,
		main="NUSE (Normalized Unscaled Standard Error)\nShould be centered on 1 (blue line)\nRed sample = out of bounds (dashed red line)",
		names=sampleNames(dataset), las=2,
		border=c("black","red")[(QC.medians$NUSE > QC.cutoffs$NUSE)+1]	
	);
	lines(x=c(0,ncol(QC$NUSE)+1), y=rep(1,2), col="blue", lty=2);
	lines(x=c(0,ncol(QC$NUSE)+1), y=rep(QC.cutoffs$NUSE,2), col="red", lty=2);
	dev.off();
}

# Create a PDF of histograms of the expression values
output.filename <- file.path(output.path, sprintf("%s_RMA_histograms.pdf", output.prefix));
cat(sprintf("Drawing histograms to '%s'.\n", output.filename));
pdf(output.filename, width=11, height=8.5);
for (i in 1:nrow(pData(dataset))) {
	hist(exprs(dataset)[, sampleNames(dataset)[i]], breaks=100, main=sampleNames(dataset)[i], xlab="log2 (expression)");
}
dev.off();

# Perform PCA
if (ncol(dataset) > 2) {
	output.filename <- file.path(output.path, sprintf("%s_PCA.pdf", output.prefix));
	cat(sprintf("Drawing Principal Component Analysis (PCA) plot to '%s'.\n", output.filename));
	dataset.pca <- prcomp(na.omit(t(scale(t(exprs(dataset))))), center=FALSE, scale=FALSE);
	percent.variance <- summary(dataset.pca)$importance["Proportion of Variance", ] * 100;
	# Create a PDF of a plot of PC2 vs PC1
	pdf("PCA_Output.pdf");
	# Note: 10% extra space is alloted in order to give more room for the sample names
	plot(
		dataset.pca$rotation[,c("PC1","PC2")], pch=21, cex=2, col="black", bg="gray",
		xlim=1.1 * range(dataset.pca$rotation[,"PC1"]), ylim=1.1 * range(dataset.pca$rotation[,"PC2"]),
		xlab=sprintf("PC1 (%2.0f%%)", percent.variance["PC1"]), ylab=sprintf("PC2 (%2.0f%%)", percent.variance["PC2"]),
		main="Principal Component Analysis (PCA)\nacross all genes in all samples"
	);
	# Write the sample names over the points
	dev.off();
}

# Check for control (sex-specific and GSTT1/GSTM1) signal and draw a heatmap
control.symbols <- c("XIST","CYorf15A","DDX3Y","KDM5D","RPS4Y1","USP9Y","UTY","GSTT1","GSTM1");
control.genes <- featureNames(dataset)[match(control.symbols, fData(dataset)$Symbol)];
names(control.genes) <- control.symbols;
control.genes <- na.omit(control.genes);
output.filename <- file.path(output.path, sprintf("%s_control_gene_heatmap.pdf", output.prefix));
cat(sprintf("Drawing heatmap of control genes to '%s'.\n", output.filename));
pdf(output.filename, width=11, height=8.5);
# Do not scale the heatmap; at least some genes should have low or high expression to set colors properly
heatmap(
	exprs(dataset)[control.genes, ], col=colorRampPalette(c("blue","white","red"))(256), revC=TRUE, margins=c(20,5),
	labRow=names(control.genes),
	Rowv=NA, Colv=NA, scale="none", main="Heatmap of absolute (unscaled) expression of control genes"
);
dev.off();

# Create a data frame summarizing QC metrics and control gene expression and write to a file
output.filename <- file.path(output.path, sprintf("%s_QC_summary.txt", output.prefix));
cat(sprintf("Writing table of QC metrics and control gene expression to '%s'.\n", output.filename));
# Add median RLE and NUSE
QC.summary <- as.data.frame(QC.medians);
colnames(QC.summary) <- paste(colnames(QC.summary));
# Add PC1 and PC2 if there were more than 2 samples in the dataset
if (ncol(dataset) > 2) QC.summary <- cbind(QC.summary, as.data.frame(dataset.pca$rotation[,c("PC1","PC2")]));
# Add the control gene expression
QC.summary <- cbind(QC.summary, t(data.frame(exprs(dataset)[control.genes,], row.names=names(control.genes), check.names=FALSE)));
# Write to a tab-delimited text file
write.table(QC.summary, output.filename, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA);

save.image('300613_QCed_Bronc_RMA_Complete.RData')

rin <- pData(dataset)$RIN
batch <- pData(dataset)$BATCH
cel <- row.names(pData(dataset))

# BRONC PCA BY BATCH & DATASET
pcashape <- c(rep(1,length(cel)))
pcacolor <- c(rep('red',length(cel)))
pcacolor[which(batch==13)]<-"orange"
pcacolor[which(batch==14)]<-"green"
pcacolor[which(batch==15)]<-"blue"
pcacolor[which(batch==16)]<-"purple"
output.filename <- file.path(output.path, sprintf("%s_PCA_by_BATCH.pdf", output.prefix));
pdf(output.filename);
	# Note: 10% extra space is alloted in order to give more room for the sample names
	plot(
		dataset.pca$rotation[,c("PC1","PC2")], pch=pcashape, cex=2, col=pcacolor, bg=pcacolor,
		xlim=1.1 * range(dataset.pca$rotation[,"PC1"]), ylim=1.1 * range(dataset.pca$rotation[,"PC2"]),
		xlab=sprintf("PC1 (%2.0f%%)", percent.variance["PC1"]), ylab=sprintf("PC2 (%2.0f%%)", percent.variance["PC2"]),
		main="Principal Component Analysis (PCA)\nacross all genes in all samples"
	);
	# Write the sample names over the points
#	text(x=dataset.pca$rotation[,"PC1"], y=dataset.pca$rotation[,"PC2"], cex=0.7, col="black", labels=sampleNames(dataset));
	legend("topleft", c("Batch12","Batch13","Batch14","Batch15","Batch16"), pch=c(rep(1,5)), col=c('red','orange','green','blue','purple'), cex=0.8)
	dev.off();

# BRONC PCA BY RLE & NUSE
pcashape <- c(rep(1,length(cel)))
pcashape[which(!match(cel,row.names(QC.summary)[which(QC.summary$RLE>0.15)])==FALSE)] <- 2
pcacolor <- c(rep('green',length(cel)))
pcacolor[which(!match(cel,row.names(QC.summary)[which(QC.summary$NUSE>1.05)])==FALSE)] <- "red"
#pcacolor[which(!match(cel,row.names(QC.summary)[which(QC.summary$RLE>0.1)]))] <- "red"
output.filename <- file.path(output.path, sprintf("%s_PCA_by_RLE_0.15_NUSE_1.05.pdf", output.prefix));
pdf(output.filename);
	# Note: 10% extra space is alloted in order to give more room for the sample names
	plot(
		dataset.pca$rotation[,c("PC1","PC2")], pch=pcashape, cex=2, col=pcacolor, bg=pcacolor,
		xlim=1.1 * range(dataset.pca$rotation[,"PC1"]), ylim=1.1 * range(dataset.pca$rotation[,"PC2"]),
		xlab=sprintf("PC1 (%2.0f%%)", percent.variance["PC1"]), ylab=sprintf("PC2 (%2.0f%%)", percent.variance["PC2"]),
		main="Principal Component Analysis (PCA)\nacross all genes in all samples"
	);
	# Write the sample names over the points
#	text(x=dataset.pca$rotation[,"PC1"], y=dataset.pca$rotation[,"PC2"], cex=0.7, col="black", labels=sampleNames(dataset));
	legend("topleft", c("NUSE > 1.05","NUSE < 1.05","RLE > 0.15","RLE < 0.15"), lty=c(1,1,NA,NA), pch=c(NA,NA,2,1), col=c('red','green','black','black'), cex=1)
#	legend("bottomright", c("Dataset 1", "Dataset 2"), pch=c(1,2), col="black", cex=0.8)
	dev.off();

#BRONC PCA BY RIN
#pcashape <- c(rep(1,length(cel))
pcacolor <- c(rep('green',length(cel)))
pcacolor[which(rin<7)] <- "orange"
pcacolor[which(rin<3)] <- "red"
output.filename <- file.path(output.path, sprintf("%s_PCA_by_RIN.pdf", output.prefix));
pdf(output.filename);
	# Note: 10% extra space is alloted in order to give more room for the sample names
	plot(
		dataset.pca$rotation[,c("PC1","PC2")], pch=1, cex=2, col=pcacolor, bg=pcacolor,
		xlim=1.1 * range(dataset.pca$rotation[,"PC1"]), ylim=1.1 * range(dataset.pca$rotation[,"PC2"]),
		xlab=sprintf("PC1 (%2.0f%%)", percent.variance["PC1"]), ylab=sprintf("PC2 (%2.0f%%)", percent.variance["PC2"]),
		main="Principal Component Analysis (PCA)\nacross all genes in all samples"
	);
	# Write the sample names over the points
#	text(x=dataset.pca$rotation[,"PC1"], y=dataset.pca$rotation[,"PC2"], cex=0.7, col="black", labels=sampleNames(dataset));
	legend("topleft", c("RIN <=3", "RIN >3, <7", "RIN >=7"), pch=c(rep(1,3)), col=c('red','orange','green'), cex=0.8)
	dev.off();
