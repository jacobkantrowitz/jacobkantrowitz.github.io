library(WGCNA)

lnames <- load("~/Meta_Analysis/PCAN/data/processed/air_all_grt_lam_lgv_ann_exp.rda")  # "all.ann" "all.exp" "air.ann" "air.exp" "lam.ann" "lam.exp" "grt.ann" "grt.exp" "lgv.ann" "lgv.exp"
lnames <- load("~/Meta_Analysis/PCAN/data/input_for_WGCNA_nonevers_no2002.rda")  # "exp.aa.c.combat"  "ann.aa.c"         "exp.aa.n.combat"  "ann.aa.n"         "exp.llg.d.combat" "ann.llg.d"        "exp.llg.n.combat" "ann.llg.n"  

air.ann <- air.ann[ air.ann$smoking_status!="never", ]
air.ann <- drop.levels(air.ann)
air.ann <- air.ann[ air.ann$filename %in% c(as.character(ann.aa.c$filename), as.character(ann.aa.n$filename)), ]

path.pcan <- "~/Meta_Analysis/PCAN/"

min.mod <- 9
# path.results <- "~/Meta_Analysis/PCAN/results/WGCNA/min9/"; min.mod <- 9
# path.results <- "~/Meta_Analysis/PCAN/results/WGCNA/min5/"; min.mod <- 5
# path.results <- "~/Meta_Analysis/PCAN/results/WGCNA/"; min.mod <- 20
# path.results <- paste0("~/Meta_Analysis/PCAN/results/WGCNA/140313_nonevers_min9_perm1_Rsq0.8/")
path.results <- paste0("~/Meta_Analysis/PCAN/results/WGCNA/140317_nonevers_min30_perm1_Rsq0.8/")







setwd(path.pcan)

### Cancer network
# combine annotations
traits <- intersect(names(air.ann), names(all.ann))
# traits <- traits[ traits!="filename"]

### NORMAL (Cancer network)

# ann.aa.n <- rbind(all.ann[all.ann$filename %in% colnames(exp.aa.n.combat), traits], 
#                   air.ann[air.ann$filename %in% colnames(exp.aa.n.combat), traits])
ann.aa.n.orig <- ann.aa.n
ann.aa.n <- ann.aa.n[, traits]
ann.aa.n$sex <- as.numeric(ann.aa.n$sex)
ann.aa.n$race <- as.numeric(ann.aa.n$race)
ann.aa.n$smoking_status <- as.numeric(ann.aa.n$smoking_status)
ann.aa.n$copd_status <- as.numeric(ann.aa.n$copd_status)
ann.aa.n$cancer_status <- as.numeric(ann.aa.n$cancer_status)

setwd(path.results)
res <- load("air_all_normal_combat_wgcna_results_1.rda")
res <- coexpress

datExpr <- t(exp.aa.n.combat)
datTraits <- ann.aa.n[ match(rownames(datExpr), ann.aa.n$filename), colnames(ann.aa.n)!="filename"]

title <- "Module-trait relationships \n (Cancer dataset, Normal network)"
filename <- paste0(path.results, "figures/WGCNA_aa_normal_module-trait_heatmap_min", min.mod, ".pdf")

res.aa.n <- res
datExpr.aa.n <- datExpr
datTraits.aa.n <- datTraits
title.aa.n <- title
filename.aa.n <- filename

### CANCER (Cancer network)
# ann.aa.c <- rbind(all.ann[all.ann$filename %in% colnames(exp.aa.c.combat), traits], 
#                   air.ann[air.ann$filename %in% colnames(exp.aa.c.combat), traits])
ann.aa.c.orig <- ann.aa.c
ann.aa.c <- ann.aa.c[, traits]
ann.aa.c$sex <- as.numeric(ann.aa.c$sex)
ann.aa.c$race <- as.numeric(ann.aa.c$race)
ann.aa.c$smoking_status <- as.numeric(ann.aa.c$smoking_status)
ann.aa.c$copd_status <- as.numeric(ann.aa.c$copd_status)
ann.aa.c$cancer_status <- as.numeric(ann.aa.c$cancer_status)

setwd(path.results)
res <- load("air_all_cancer_combat_wgcna_results_1.rda")
res <- coexpress

datExpr <- t(exp.aa.c.combat)
datTraits <- ann.aa.c[ match(rownames(datExpr), ann.aa.c$filename), colnames(ann.aa.c)!="filename"]

title <- "Module-trait relationships \n (Cancer dataset, Cancer network)"
filename <- paste0(path.results, "figures/WGCNA_aa_cancer_module-trait_heatmap_min", min.mod, ".pdf")

res.aa.c <- res
datExpr.aa.c <- datExpr
datTraits.aa.c <- datTraits
title.aa.c <- title
filename.aa.c <- filename

### Premalignancy network
# lnames <- load("~/Meta_Analysis/REAGENT/data/processed/all_air_lam_grt_lgv_ann_exp.rda")
names(lam.ann)[ names(lam.ann)=="copd_status_der"] <- "copd_status"
# combine annotations
traits <- Intersect(names(lam.ann), names(lgv.ann), names(grt.ann))
# traits <- traits[ traits!="filename"]

####
lgv.ann$cancer_status <- as.factor(gsub("no cancer", "no", lgv.ann$cancer_status))
lgv.ann$cancer_status <- as.factor(gsub("cancer", "yes", lgv.ann$cancer_status))

lgv.ann$copd_status <- as.factor(gsub("no COPD", "no", lgv.ann$copd_status))
lgv.ann$copd_status <- as.factor(gsub("COPD", "yes", lgv.ann$copd_status))

####


### NORMAL (Premalignancy network)
# ann.llg.n <- rbind(lam.ann[lam.ann$filename %in% colnames(exp.llg.n.combat), traits], 
#                    grt.ann[grt.ann$filename %in% colnames(exp.llg.n.combat), traits], 
#                   lgv.ann[lgv.ann$filename %in% colnames(exp.llg.n.combat), traits])
ann.llg.n.orig <- ann.llg.n
ann.llg.n <- ann.llg.n[, traits]
ann.llg.n$sex <- as.numeric(ann.llg.n$sex)
ann.llg.n$race <- as.numeric(ann.llg.n$race)
ann.llg.n$smoking_status <- as.numeric(ann.llg.n$smoking_status)
ann.llg.n$copd_status <- as.numeric(ann.llg.n$copd_status)
ann.llg.n$cancer_status <- as.numeric(ann.llg.n$cancer_status)
ann.llg.n$dysplasia_status <- as.numeric(ann.llg.n$dysplasia_status)

setwd(path.results)
res <- load("lam_grt_lgv_normal_combat_wgcna_results_1.rda")
res <- coexpress

datExpr <- t(exp.llg.n.combat); dim(datExpr)
datTraits <- ann.llg.n[ match(rownames(datExpr), ann.llg.n$filename), colnames(ann.llg.n)!="filename"]; dim(datTraits)

title <- "Module-trait relationships \n (Premalignancy dataset, Normal network)"
filename <- paste0(path.results, "figures/WGCNA_llg_normal_module-trait_heatmap_min", min.mod, ".pdf")

res.llg.n <- res
datExpr.llg.n <- datExpr
datTraits.llg.n <- datTraits
title.llg.n <- title
filename.llg.n <- filename

### DYSPLASIA (Premalignancy network)
# ann.llg.d <- rbind(lam.ann[lam.ann$filename %in% colnames(exp.llg.d.combat), traits], 
#                    grt.ann[grt.ann$filename %in% colnames(exp.llg.d.combat), traits], 
#                    lgv.ann[lgv.ann$filename %in% colnames(exp.llg.d.combat), traits])
ann.llg.d.orig <- ann.llg.d
ann.llg.d <- ann.llg.d[, traits]
ann.llg.d$sex <- as.numeric(ann.llg.d$sex)
ann.llg.d$race <- as.numeric(ann.llg.d$race)
ann.llg.d$smoking_status <- as.numeric(ann.llg.d$smoking_status)
ann.llg.d$copd_status <- as.numeric(ann.llg.d$copd_status)
ann.llg.d$cancer_status <- as.numeric(ann.llg.d$cancer_status)
# ann.llg.d$DYSGR_CAT <- as.numeric(ann.llg.d$DYSGR_CAT)

setwd(path.results)
res <- load("lam_grt_lgv_dysplasia_combat_wgcna_results_1.rda")
res <- coexpress

datExpr <- t(exp.llg.d.combat); dim(datExpr)
datTraits <- ann.llg.d[ match(rownames(datExpr), ann.llg.d$filename), colnames(ann.llg.d)!="filename"]; dim(datTraits)

title <- "Module-trait relationships \n (Premalignancy dataset, Dysplasia network)"
filename <- paste0(path.results, "figures/WGCNA_llg_dysplasia_module-trait_heatmap_min", min.mod, ".pdf")

res.llg.d <- res
datExpr.llg.d <- datExpr
datTraits.llg.d <- datTraits
title.llg.d <- title
filename.llg.d <- filename

### UNIVERSAL FROM THIS POINT ON
datExpr <- datExpr.aa.n
datTraits <- datTraits.aa.n
res <- res.aa.n
title <- title.aa.n
filename <- filename.aa.n

datExpr <- datExpr.aa.c
datTraits <- datTraits.aa.c
res <- res.aa.c
title <- title.aa.c
filename <- filename.aa.c

datExpr <- datExpr.llg.n
datTraits <- datTraits.llg.n
res <- res.llg.n
title <- title.llg.n
filename <- filename.llg.n

datExpr <- datExpr.llg.d
datTraits <- datTraits.llg.d
res <- res.llg.d
title <- title.llg.d
filename <- filename.llg.d

setwd(path.pcan)
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0 <- res$MEs
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "pairwise.complete.obs", method="pearson");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 
moduleTraitCor.aa.n <- moduleTraitCor
moduleTraitCor.aa.c <- moduleTraitCor
moduleTraitCor.llg.n <- moduleTraitCor
moduleTraitCor.llg.d <- moduleTraitCor

## Save data
# Min module size 5
# save(datExpr.aa.c, datExpr.aa.n, datExpr.llg.d, datExpr.llg.n, datTraits.aa.c, datTraits.aa.n, datTraits.llg.d, datTraits.llg.n, file="~/Meta_Analysis/PCAN/data/processed/datTraitsExprs_aa_llg_n_c_d_minModuleSize5.rda")
# save(moduleTraitCor.aa.c, moduleTraitCor.aa.n, moduleTraitCor.llg.d, moduleTraitCor.llg.n, file="~/Meta_Analysis/PCAN/data/processed/moduleTraitCors_aa_llg_n_c_d_minModuleSize5.rda")
# save(res.aa.c, res.aa.n, res.llg.d, res.llg.n, file="~/Meta_Analysis/PCAN/data/processed/coexpresss_aa_llg_n_c_d_minModuleSize5.rda")

# Min module size 9
# save(datExpr.aa.c, datExpr.aa.n, datExpr.llg.d, datExpr.llg.n, datTraits.aa.c, datTraits.aa.n, datTraits.llg.d, datTraits.llg.n, file="~/Meta_Analysis/PCAN/data/processed/datTraitsExprs_aa_llg_n_c_d_minModuleSize9.rda")
# save(moduleTraitCor.aa.c, moduleTraitCor.aa.n, moduleTraitCor.llg.d, moduleTraitCor.llg.n, file="~/Meta_Analysis/PCAN/data/processed/moduleTraitCors_aa_llg_n_c_d_minModuleSize9.rda")
# save(res.aa.c, res.aa.n, res.llg.d, res.llg.n, file="~/Meta_Analysis/PCAN/data/processed/coexpresss_aa_llg_n_c_d_minModuleSize9.rda")

# Min module size 20
# save(datExpr.aa.c, datExpr.aa.n, datExpr.llg.d, datExpr.llg.n, datTraits.aa.c, datTraits.aa.n, datTraits.llg.d, datTraits.llg.n, file="~/Meta_Analysis/PCAN/data/processed/datTraitsExprs_aa_llg_n_c_d_minModuleSize20.rda")
# save(moduleTraitCor.aa.c, moduleTraitCor.aa.n, moduleTraitCor.llg.d, moduleTraitCor.llg.n, file="~/Meta_Analysis/PCAN/data/processed/moduleTraitCors_aa_llg_n_c_d_minModuleSize20.rda")
# save(res.aa.c, res.aa.n, res.llg.d, res.llg.n, file="~/Meta_Analysis/PCAN/data/processed/coexpresss_aa_llg_n_c_d_minModuleSize20.rda")

# # Min module size 30
# save(datExpr.aa.c, datExpr.aa.n, datExpr.llg.d, datExpr.llg.n, datTraits.aa.c, datTraits.aa.n, datTraits.llg.d, datTraits.llg.n, file="~/Meta_Analysis/PCAN/data/processed/datTraitsExprs_aa_llg_n_c_d_minModuleSize30.rda")
# save(moduleTraitCor.aa.c, moduleTraitCor.aa.n, moduleTraitCor.llg.d, moduleTraitCor.llg.n, file="~/Meta_Analysis/PCAN/data/processed/moduleTraitCors_aa_llg_n_c_d_minModuleSize30.rda")
# save(res.aa.c, res.aa.n, res.llg.d, res.llg.n, file="~/Meta_Analysis/PCAN/data/processed/coexpresss_aa_llg_n_c_d_minModuleSize30.rda")

# # Min module size 09 140313
save(datExpr.aa.c, datExpr.aa.n, datExpr.llg.d, datExpr.llg.n, datTraits.aa.c, datTraits.aa.n, datTraits.llg.d, datTraits.llg.n, file="~/Meta_Analysis/PCAN/data/processed/140317_datTraitsExprs_aa_llg_n_c_d_minModuleSize30.rda")
save(moduleTraitCor.aa.c, moduleTraitCor.aa.n, moduleTraitCor.llg.d, moduleTraitCor.llg.n, file="~/Meta_Analysis/PCAN/data/processed/140317_moduleTraitCors_aa_llg_n_c_d_minModuleSize30.rda")
save(res.aa.c, res.aa.n, res.llg.d, res.llg.n, file="~/Meta_Analysis/PCAN/data/processed/140317_coexpresss_aa_llg_n_c_d_minModuleSize30.rda")



pdf(file=filename, height=10)
#   sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = title)

dev.off()


### Calculate numerical overlaps with sex and smoking gene sets
# use check_smoking_overlap_v2.R

## SEX
# view modules with highest sex correlation
head(moduleTraitCor[order(abs(moduleTraitCor[, "sex"]), decreasing=TRUE), ])

# view modules that contain signature sex genes
x <- res$modules[ res$modules$Gene %in% sex.sig.unlisted, ]; x
sort(table(res$modules$Color[ res$modules$Gene %in% sex.sig.unlisted]))

# view size of module that includes sex genes
sort(table(res$modules$Color[ res$modules$Color %in% x$Color]))

## SMOKING
# view modules with highest sex correlation
head(moduleTraitCor[order(abs(moduleTraitCor[, "smoking_status"]), decreasing=TRUE), ])

# view modules that contain signature sex genes
smk.sig.unlisted <- paste0(smk.ens, "_at")
smk.sig.unlisted <- smk.sig.unlisted[ !smk.sig.unlisted=="NA_at"]
x <- res$modules[ res$modules$Gene %in% smk.sig.unlisted, ]; x
y<-sort(table(res$modules$Color[ res$modules$Gene %in% smk.sig.unlisted]))

# view size of module that includes sex genes
z<-sort(table(res$modules$Color[ res$modules$Color %in% x$Color]))
zz <- merge(melt(y), melt(z), by="Var.1", all=TRUE)
zzz <- zz[ order(zz$value.x, decreasing=TRUE), ]
