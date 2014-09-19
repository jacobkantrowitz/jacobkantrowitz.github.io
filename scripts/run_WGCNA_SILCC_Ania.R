# load data; ania is loading matrices
#lnames <- load("~/Meta_Analysis/PCAN/data/input_for_WGCNA_nonevers_no2002.rda")

# set the working directory to the directory for experiments started 2014-08-13
setwd("/restricted/projectnb/pulmarray/LinGA_protected/Allegro/COPD_Cancer/experiments/2014-08-13")

# load the data using a script that removes patients with:
#  Cancer = NA
#  SMK    = 3
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/scripts/AllegroSetup.R")
holdEset <- removeBioReps(eset)
eset <- holdEset

#### START ANIA'S CODE

.libPaths(c("/unprotected/projects/cbmhive/R_packages/R-3.0.0//", .libPaths()))
print(.libPaths)
.libPaths(c("/usr3/graduate/kantro/R_library", .libPaths()))
source("/restricted/projectnb/pulmarray/LinGA_protected/Allegro/COPD_Cancer/scripts/Ania_WGCNA_wrapper_permutations.R")
library(WGCNA)#, lib.loc="~/R_packages/R-3.0.0/")
library(heatmap3)#, lib.loc="/home/aniat/R_packages/R-2.15.1/")


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# args[1] = nCores (number of cores to use)
nCores <- as.numeric(args[1])
# args[2] = nPerm (counter for permutation)
nPerm <- as.numeric(args[2])

# standard from WGCNA - finds the most biologically relevant modules, per WGCNA paper
# gender module will likely be smaller than this, 5 genes likely associated with gender
minModSize = 30

## Estimate soft threshold parameter to achieve scale-free topology
pdf("Allegro_WGCNA_SoftThreshold_plots.pdf")
#where does this field exp.aa.c.combat come from?
power = plotSoftThreshold(datExpr=exprs(eset), RsquaredCut=0.8)
dev.off()

cat(paste0("power = ", power))

# exp.aa.c.combat in ania's case is a numeric matrix, gene (rows) x sample (column)

## Perform WGCNA 
# for Airway try smaller modules to detect sex module
# coexpress = wgcna(data, softPower=power, minModuleSize=10, saveDissTOM=paste0(path, study, "_",power, "_", Rsq, "_TOM.Rda"), plotModuleHeatmaps=paste0(path, study, "_module_heatmaps.pdf"))
# save(coexpress, file=paste0(path, study, "_wgcna_results_", power, "_", Rsq, "minModuleSize10.rda"))
# fro Airway, try using non-ComBatted data
# coexpress = wgcna(data, softPower=power, minModuleSize=20, saveDissTOM=paste0(path, study, "_",power, "_", Rsq, "_TOM_noComBat.Rda"), plotModuleHeatmaps=paste0(path, study, "_module_heatmaps_noComBat_allGenes.pdf"))
# save(coexpress, file=paste0(path, study, "_wgcna_results_", power, "_", Rsq, "_noComBat_allGenes.rda"))
coexpress = wgcna(exprs(eset), softPower=power, minModuleSize = minModSize, 
                  saveAdjacency="Adjacency.Rda", 
                  saveDissTOM="TOM.Rda", 
                  plotModuleHeatmaps="Module_Heatmaps.pdf",
                  nCores=nCores)
save(coexpress, file="WGCNA_results.Rda")
gc()  
