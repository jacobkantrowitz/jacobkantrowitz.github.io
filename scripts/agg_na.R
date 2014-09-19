#author: Ana Pavel

# ls *pred*.txt | cut -d "_" -f 4 | cut  -d "." -f 1 > list.txt
# paste -s -d "," list.txt > list2.txt

.libPaths(c(sprintf("/unprotected/projects/cbmhive/R_packages/R-%s", getRversion()), .libPaths()));

# prefix <- "p1_mirna_training_95/140125_"
prefix <- "~/Meta_Analysis/lungevity/lungevity-biomarker/results/140130"

iterations <- c(100,10,11,12,13,14,15,16,17,18,19,1,20,21,22,23,24,25,26,27,28,29,2,30,31,32,33,34,35,36,37,38,39,3,40,41,42,43,44,45,46,47,48,49,4,50,51,52,53,54,55,56,57,58,59,5,60,61,62,63,64,65,66,67,68,69,6,70,71,72,73,74,75,76,77,78,79,7,80,81,82,83,84,85,86,87,88,89,8,90,91,92,93,94,95,96,97,98,99,9)
N = length(iterations)

filename <- paste(prefix, "_predictions_", iterations[1], ".txt", sep="")
data.agg <- as.matrix(read.table(filename, check.names=F, header=F))
data.cnt <- 0 + !is.na(data.agg)
data.agg[is.na(data.agg)] <- 0
for(i in 2:N) {
    filename <- paste(prefix, "_predictions_", iterations[i], ".txt", sep="")
    data.perf <- as.matrix(read.table(filename, check.names=F, header=F))
    data.perf.cnt <- 0 + !is.na(data.perf)
    data.perf[is.na(data.perf)] <- 0
    
    data.agg <- data.agg + data.perf
    data.cnt <- data.cnt + data.perf.cnt
}

data.agg <- data.agg/data.cnt

# colnames(data.agg) <- c("xval.iteration",
# "dataset", "sample.filter", "norm", "transform", "feature.filter", "phenotype", "feat.select", "num.feat", "model.predict",
# "train.ACC", "train.SENS", "train.SPEC", "train.PPV", "train.NPV", "train.MCC", "train.AUC", "train.MAQC2",
# "test.ACC", "test.SENS", "test.SPEC", "test.PPV", "test.NPV", "test.MCC", "test.AUC", "test.MAQC2")
# 
# write.table(data.agg, paste(prefix, "p1_mirna_training_95.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)

# m <- read.table(file.out); #m<-m[-1,]
m <- as.data.frame(data.agg)
colnames(m) <- c("iteration", "sample_filter", "transformation", "data_type", "gene_filter", 
                 "phenotype", "feature_selection", "number_features", "prediction_mode", 
                 "ACC.tr", "SENS.tr", "SPEC.tr", "PPV.tr", "NPV.tr", "MCC.tr", "AUC.tr", "MAQC2.tr", 
                 "ACC.ts", "SENS.ts", "SPEC.ts", "PPV.ts", "NPV.ts", "MCC.ts", "AUC.ts", "MAQC2.ts", 
                 "tr.0", "tr.1", "ts.0", "ts.1")

m[, 2:9] <- apply(m[, 2:9], 2, as.character)

for (row in 1:nrow(m)) {
  for (col in 2:9) {
    m[row, col] <- param.list[[col-1]][as.numeric(m[row, col])]
  }
}

m[, c(2,3,4,6,7,9)] <- lapply(m[, c(2,3,4,6,7,9)], as.factor)
m[, c(5,8)] <- lapply(m[, c(5,8)], as.numeric)
summary(m)

write.table(m, paste0(prefix, "_predictions.agg.readable.txt"), sep = "\t", quote = F, row.names = F)





print(warnings())
