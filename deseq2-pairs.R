library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(genefilter)

huh7 <- as.matrix(read.csv("~/cchf/tmp/hepg2/gene_count_matrix.csv", row.names="gene_id"))
colData_huh7 <- read.table("~/cchf/tmp/hepg2/cchv.key", row.names=1, sep="\t", h=T)
huh7 <- as.matrix(read.csv("~/cchf/tmp/huh7/gene_count_matrix.csv", row.names="gene_id"))
colData_huh7 <- read.table("~/cchf/tmp/huh7/cchv.key", row.names=1, sep="\t", h=T)
colData_huh7$status <- relevel(colData_huh7$status, ref="uninfected")
colData_huh7$timepoint <- as.factor(colData_huh7$timepoint)
colData_huh7$group <- as.factor(paste(colData_huh7$status, colData_huh7$timepoint, sep="_"))
colData_huh7 <- colData_huh7[-1]
str(colData_huh7)
dds_huh7 <- DESeqDataSetFromMatrix(countData = huh7,
                                   colData = colData_huh7,
                                   design = ~ group)
dds_huh7 <- dds_huh7[ rowSums(counts(dds_huh7)) > 1, ]
dds <- DESeq(dds_huh7, betaPrior = F)

resultsNames(dds)

thr <- results(dds, contrast = c("group", "infected_3", "uninfected_3"))
thrShrunk <- lfcShrink(dds, contrast = c("group", "infected_3", "uninfected_3"), res = thr)

sev <- results(dds, contrast = c("group", "infected_7", "uninfected_7"))
sevShrunk <- lfcShrink(dds, contrast = c("group", "infected_7", "uninfected_7"), res = sev)
head(sev[order(sev$padj),], 10)

plotMA(thr, ylim = c(-5,5))
plotMA(thrShrunk, ylim = c(-5, 5))

top30 <- head(thr[order(thr$padj),], 30)
top30
top30shrunk <- head(thr[order(thrShrunk$padj),], 100)
top30shrunk

plotCounts(dds, gene="MSTRG.5463", intgroup = c("group"))

top30rows <- head(order(rowVars(assay(top30shrunk)), decreasing = T), 30)

rld <- rlog(dds_huh7, blind = F)
top30genes <- head(order(rowVars(assay(rld)), decreasing = T), 100)
mat <- assay(rld)[rownames(top30shrunk), c("Lib_25", "Lib_26", "Lib_27", "Lib_28", "Lib_29", "Lib_30")]
mat
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("group")])
pheatmap(mat, main="hepg2")
