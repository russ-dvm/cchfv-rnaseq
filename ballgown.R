library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)

pheno_data <- read.table("~/cchf/ref_files/cchv.key", h = T)

dataDir <- "~/cchf/6_reestimate/"
samplePattern <- "Lib"


bg <- ballgown(dataDir = dataDir, samplePattern = samplePattern, pData = pheno_data)

bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = T)
fpkm <- texpr(bg_filt, meas="FPKM")
fpkm <- log2(fpkm+1)
# fpkm <- as.data.frame(fpkm)

ballgown::geneNames(bg_filt)[grep("MX1", ballgown::geneNames(bg_filt))]
plot(fpkm[147929,] ~ pheno_data$status, border=c(1,2), las=2)
plot(fpkm[147912,] ~ pheno_data$status, border=c(1,2), las=2)
fpkm[147912,]

plotTranscripts(ballgown::geneIDs(bg)[147912], bg, sample=c("Lib_10", "Lib_11", "Lib_12"))

ggplot(fpkm, aes(x=))
t(fpkm[147932,])
t(as.matrix(fpkm[147912,]))


ggplot(fpkm, aes(x=log(FPKM.Lib_10+1), y = log(FPKM.Lib_2+1))) + geom_point() + theme_bw()
cor(fpkm$FPKM.Lib_1, fpkm$FPKM.Lib_2)
cor(fpkm$FPKM.Lib_1, fpkm$FPKM.Lib_3)
cor(fpkm$FPKM.Lib_3, fpkm$FPKM.Lib_2)

cor(fpkm$FPKM.Lib_4, fpkm$FPKM.Lib_5)
cor(fpkm$FPKM.Lib_4, fpkm$FPKM.Lib_6)
cor(fpkm$FPKM.Lib_6, fpkm$FPKM.Lib_5)

cor(fpkm$FPKM.Lib_7, fpkm$FPKM.Lib_8)
cor(fpkm$FPKM.Lib_7, fpkm$FPKM.Lib_9)
cor(fpkm$FPKM.Lib_9, fpkm$FPKM.Lib_8)

cor(fpkm$FPKM.Lib_10, fpkm$FPKM.Lib_11)
cor(fpkm$FPKM.Lib_10, fpkm$FPKM.Lib_12)
cor(fpkm$FPKM.Lib_12, fpkm$FPKM.Lib_11)

cor(fpkm$FPKM.Lib_13, fpkm$FPKM.Lib_14)
cor(fpkm$FPKM.Lib_13, fpkm$FPKM.Lib_15)
cor(fpkm$FPKM.Lib_15, fpkm$FPKM.Lib_14)

cor(fpkm$FPKM.Lib_16, fpkm$FPKM.Lib_17)
cor(fpkm$FPKM.Lib_16, fpkm$FPKM.Lib_18)
cor(fpkm$FPKM.Lib_18, fpkm$FPKM.Lib_17)

cor(fpkm$FPKM.Lib_19, fpkm$FPKM.Lib_20)
cor(fpkm$FPKM.Lib_19, fpkm$FPKM.Lib_21)
cor(fpkm$FPKM.Lib_21, fpkm$FPKM.Lib_20)

cor(fpkm$FPKM.Lib_22, fpkm$FPKM.Lib_23)
cor(fpkm$FPKM.Lib_22, fpkm$FPKM.Lib_24)
cor(fpkm$FPKM.Lib_24, fpkm$FPKM.Lib_23)

cor(fpkm$FPKM.Lib_25, fpkm$FPKM.Lib_26)
cor(fpkm$FPKM.Lib_25, fpkm$FPKM.Lib_27)
cor(fpkm$FPKM.Lib_27, fpkm$FPKM.Lib_26)

cor(fpkm$FPKM.Lib_28, fpkm$FPKM.Lib_29)
cor(fpkm$FPKM.Lib_28, fpkm$FPKM.Lib_30)
cor(fpkm$FPKM.Lib_30, fpkm$FPKM.Lib_29)

cor(fpkm$FPKM.Lib_31, fpkm$FPKM.Lib_32)
cor(fpkm$FPKM.Lib_31, fpkm$FPKM.Lib_33)
cor(fpkm$FPKM.Lib_33, fpkm$FPKM.Lib_32)

cor(fpkm$FPKM.Lib_34, fpkm$FPKM.Lib_35)
cor(fpkm$FPKM.Lib_34, fpkm$FPKM.Lib_36)
cor(fpkm$FPKM.Lib_36, fpkm$FPKM.Lib_35)





ballgown::transcriptNames(bg_filt)


