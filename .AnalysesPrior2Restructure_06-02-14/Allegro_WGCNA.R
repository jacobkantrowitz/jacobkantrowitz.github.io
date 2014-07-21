# Allegro_WGCNA
# Jake Kantrowitz
# 11/04/13
# Script to perform WGCNA on Allegro Data

# Set the working directory to the location of the data
setwd("/Users/jacobkantrowitz/Dropbox/BUSM/Spira_Lab/Data")

# Load the necessary libraries and set required options
library(WGCNA)
options(stringsAsFactors = FALSE)

# Read in the data
print("Reading in data...")
data = readRDS("Bronc_708_First_CEL.rds")
# Extract the phenotype and expression data
phenData = pData(data)
exprData = exprs(data)
# Transpose expression data matrix for working with WGCNA functions
exprData = t(exprData)

nGenes = dim(exprData)[1]

# Appropriatley coerce categorical phenotypes into factors
# Column 1 - BATCH
phenData$BATCH <- factor(phenData$BATCH)
# Column 3 - CANCER
phenData$CANCER <- factor(phenData$CANCER)
# Column 4 - SMK
phenData$SMK <- factor(phenData$SMK)
# Column 5 - COPD
phenData$COPD <- factor(phenData$COPD)
# Column 6 - GENDER
phenData$GENDER <- factor(phenData$GENDER)

# 1 checking data for excessive missing values and identification of outlier microarray samples
gsg = goodSamplesGenes(exprData, verbose = 3)

# If gsg$allOK is FALSE, all genes have not passed the cuts/filtering. This means we need to
# remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exprData = exprData[gsg$goodSamples, gsg$goodGenes]
}

# Next we cluster the samples to see if there are any obvious outliers. We use the function
# 'flashClust' that provides faster hierarchical clutering than the standard function 'hclust'
sampleTree = flashClust(dist(exprData), method="average")

sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)



# There may be some outliers - can remove  by hand or use an automatic approach
# Choose a height cut that will remove the offending sample, say 70 in this case (add a line there)
# and use a branch cut at that height

# add the line above which to cut
cutHeight = 70
abline(h=cutHeight, col="Red")

# Determine the cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
#table(clust)

# clust 1 contains the samples we want to keep
keepSamples = (clust==1)
exprData = exprData[keepSamples,]
exprData = data.frame(exprData)
nGenes = ncol(exprData)
nSamples = nrow(exprData)

# Matching Clinical Trait Data
# We now look at the trait data and match the samples for which they were measured to the
# expression samples

# Pull the clinical traits into a new data frame
traitData = phenData[,3:7]

##########################################################################################
##########################################################################################
##########################################################################################

# tutorial for WGCNA # 1

# Matching Clinical Trait Data
# We now look at the trait data and match the samples for which they were measured to the
# expression samples

traitData = read.csv("ClinicalTraits.csv")
#dim(traitData)
#names(traitData)

# remove columns that hold information we do not need
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36)]
#dim(allTraits)
#names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$Mice)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()

# We now have the epxression data in the variable datExpr and the corresponding clinical traits in the
# variable datTraits. Before we continue with network construction and module detection, we visualize
# how the clinical traits relate to the sample dendogram

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method="average")
# Convert traits to a color representation - white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendogram and the colors underneath
plotDendroAndColors(sampleTree2, traitColors, groupLabels=names(datTraits), main = "Sample dendogram and heatmap")

# The last step is to save the relevant expression and trait data for use in the next steps of the tutorials
save(datExpr, datTraits, file="FemaleLiver-01-dataInput.RData")

################################################################################
################################################################################
# Tutorial Next Step - Network analysis of liver expression data in female mice
# Automatic network construction and module detection



