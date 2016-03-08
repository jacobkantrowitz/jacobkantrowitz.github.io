setwd("/protected/projects/pulmarray/Allegro/COPD_Cancer/experiments/Expmt_0045_Hallmarks_in_COPD/")
source("/protected/projects/pulmarray/Allegro/COPD_Cancer/scripts/AllegroSetup.R")

data <- readRDS("dataForGSVACorEmpiric.RDS")

target <- data$eset
target <- target[!is.na(fData(target)$Symbol), ]
featureNames(target) <- fData(target)$Symbol
temp <- geneset.cor.test(data$genesets$hypoxia, data$genesets$atf4up, target, num.cors=1000)
save(temp, file="empiricCorrelationsHypoxiaATF41000-10.RData")

#load("empiricCorrelationsHypoxiaATF41000-01.RData")
#load("empiricCorrelationsHypoxiaATF41000-02.RData")
load("empiricCorrelationsHypoxiaATF41000-03.RData")
temp3 <- temp
load("empiricCorrelationsHypoxiaATF41000-04.RData")
temp4 <- temp
load("empiricCorrelationsHypoxiaATF41000-05.RData")
temp5 <- temp
load("empiricCorrelationsHypoxiaATF41000-06.RData")
temp6 <- temp
#load("empiricCorrelationsHypoxiaATF41000-07.RData")
load("empiricCorrelationsHypoxiaATF41000-08.RData")
temp8 <- temp
load("empiricCorrelationsHypoxiaATF41000-09.RData")
temp9 <- temp
load("empiricCorrelationsHypoxiaATF41000-10.RData")
temp10 <- temp
load("empiricCorrelationsHypoxiaATF41000.RData")
temp1 <- temp
