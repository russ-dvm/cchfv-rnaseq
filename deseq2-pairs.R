library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(genefilter)


#### DATA IMPORT ####
##ensembl ids
ens <- read.table("~/cchf/7_deseq/ensembl_ids.txt", h=F, sep = "\t", stringsAsFactors = F)
colnames(ens) <-  c("MSTRG", "gene_id", "transcript_id", "gene_name")

## Import & refactor HEPG2 data
hepg2 <- as.matrix(read.csv("~/cchf/7_deseq/hep/gene_count_matrix.csv", row.names="gene_id"))
colData_hepg2 <- read.table("~/cchf/7_deseq/hep/cchv.key", row.names=1, sep="\t", h=T)
colData_hepg2$timepoint <- as.factor(colData_hepg2$timepoint)

## Create groups and set ref level
colData_hepg2$group <- as.factor(paste(colData_hepg2$status, colData_hepg2$timepoint, sep="_"))
colData_hepg2$group <- relevel(colData_hepg2$group, ref = "uninfected_1")

## Import & refactor HUH7 data
huh7 <- as.matrix(read.csv("~/cchf/7_deseq/huh/gene_count_matrix.csv", row.names="gene_id"))
colData_huh7 <- read.table("~/cchf/7_deseq/huh/cchv.key", row.names=1, sep="\t", h=T)
colData_huh7$timepoint <- as.factor(colData_huh7$timepoint)

## Create groups and set ref level
colData_huh7$group <- as.factor(paste(colData_huh7$status, colData_huh7$timepoint, sep="_"))
colData_huh7$group <- relevel(colData_huh7$group, ref = "uninfected_1")


#### Create the DeSEQ objects ####

# HEPG2
dds_hepg2 <- DESeqDataSetFromMatrix(countData = hepg2,
                                    colData = colData_hepg2,
                                    design = ~ group)
dds_hepg2 <- dds_hepg2[ rowSums(counts(dds_hepg2)) > 1, ]


# HUH7
dds_huh7 <- DESeqDataSetFromMatrix(countData = huh7,
                                   colData = colData_huh7,
                                   design = ~ group)
dds_huh7 <- dds_huh7[ rowSums(counts(dds_huh7)) > 1, ]


#### Run DESeq ####
# HEPG2
dds_hepg2 <- DESeq(dds_hepg2, betaPrior = F)
resultsNames(dds_hepg2)
# HUH7
dds_huh7 <- DESeq(dds_huh7, betaPrior = F)
resultsNames(dds_huh7)

#### DIRECT PAIRWISE COMPARISONS ####
## Set parameters
## Played a bit with the lfc - placing at 2 (ie requiring a 2^2 = 4 fold change for significance reduces the number of results drastically - up to 0 in some cases.)
## Defaults: lfc <- 0, alpha <- 0.05
lfc <- 0.585
alpha <- 0.05


## HEPG2
oneHep <- results(dds_hepg2, contrast = c("group", "infected_1", "uninfected_1"), alpha = alpha, lfcThreshold = lfc)
thrHep <- results(dds_hepg2, contrast = c("group", "infected_3", "uninfected_3"), alpha = alpha, lfcThreshold = lfc)
sevHep <- results(dds_hepg2, contrast = c("group", "infected_7", "uninfected_7"), alpha = alpha, lfcThreshold = lfc)

## Shrink HEPG2
oneHepShrunk <- lfcShrink(dds_hepg2, contrast = c("group", "infected_1", "uninfected_1"), res = oneHep)
thrHepShrunk <- lfcShrink(dds_hepg2, contrast = c("group", "infected_3", "uninfected_3"), res = thrHep)
sevHepShrunk <- lfcShrink(dds_hepg2, contrast = c("group", "infected_7", "uninfected_7"), res = sevHep)

## HUH7
oneHuh <- results(dds_huh7, contrast = c("group", "infected_1", "uninfected_1"), alpha = alpha, lfcThreshold = lfc)
thrHuh <- results(dds_huh7, contrast = c("group", "infected_3", "uninfected_3"), alpha = alpha, lfcThreshold = lfc)
sevHuh <- results(dds_huh7, contrast = c("group", "infected_7", "uninfected_7"), alpha = alpha, lfcThreshold = lfc)
## Shrink HUH7
oneHuhShrunk <- lfcShrink(dds_huh7, contrast = c("group", "infected_1", "uninfected_1"), res = oneHuh)
thrHuhShrunk <- lfcShrink(dds_huh7, contrast = c("group", "infected_3", "uninfected_3"), res = thrHuh)
sevHuhShrunk <- lfcShrink(dds_huh7, contrast = c("group", "infected_7", "uninfected_7"), res = sevHuh)


#### SIGNIFICANT GENES ONLY ####
sigHepOne <- subset(oneHep, oneHep$padj < 0.05)
sigHepThr <- subset(thrHep, thrHep$padj < 0.05)
sigHepSev <- subset(sevHep, sevHep$padj < 0.05)

sigHuhOne <- subset(oneHuh, oneHuh$padj < 0.05)
sigHuhThr <- subset(thrHuh, thrHuh$padj < 0.05)
sigHuhSev <- subset(sevHuh, sevHuh$padj < 0.05)

sigList <- list("HepG2 Day 1" = sigHepOne, "HepG2 Day 3" = sigHepThr, "HepG2 Day 7" = sigHepSev, "Huh7 Day 1" = sigHuhOne, "Huh7 Day 3" = sigHuhThr, "Huh7 Day 7" = sigHuhSev)
sigListDF <- lapply(sigList, data.frame)


#### MA PLOTS ####
## Arrange plots
par(mfrow=c(3,2))
## HEPG2
plotMA(oneHep, ylim = c(-5,5), main = "HEPG2, Day 1, raw")
plotMA(oneHepShrunk, ylim = c(-5, 5), main = "HEPG2, Day 1, Shrunk")
plotMA(thrHep, ylim = c(-5,5), main = "HEPG2, Day 3, raw")
plotMA(thrHepShrunk, ylim = c(-5, 5), main = "HEPG2, Day 3, Shrunk")
plotMA(sevHep, ylim = c(-5,5), main = "HEPG2, Day 7, raw")
plotMA(sevHepShrunk, ylim = c(-5, 5), main = "HEPG2, Day 7, Shrunk")

##HUH7
plotMA(oneHuh, ylim = c(-5,5), main = "HUH7, Day 1, raw")
plotMA(oneHuhShrunk, ylim = c(-5, 5), main = "HUH7, Day 1, Shrunk")
plotMA(thrHuh, ylim = c(-5,5), main = "HUH7, Day 3, raw")
plotMA(thrHuhShrunk, ylim = c(-5, 5), main = "HUH7, Day 3, Shrunk")
plotMA(sevHuh, ylim = c(-5,5), main = "HUH7, Day 7, raw")
plotMA(sevHuhShrunk, ylim = c(-5, 5), main = "HUH7, Day 7, Shrunk")


#### PCA  PLOTS ####
vstHep <- vst(dds_hepg2)
vstHuh <- vst(dds_huh7)

library(tidyr)
pcaHep <- plotPCA(vstHep, intgroup = c("group"), returnData = T)
varHep <- round(100 * attr(pcaHep, "percentVar"))
pcaHep <- separate(pcaHep, group, c("status", "timepoint"), sep = "_")

pcaHuh <- plotPCA(vstHuh, intgroup = c("group"), returnData = T)
varHuh <- round(100 * attr(pcaHuh, "percentVar"))
pcaHuh <- separate(pcaHuh, group, c("status", "timepoint"), sep = "_")

a <- ggplot(pcaHep, aes(x = PC1, y = PC2, colour = timepoint, shape = status)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", varHep[1], "% variance")) + 
  ylab(paste0("PC2: ", varHep[2], "% variance")) + 
  coord_fixed() + 
  ggtitle("HEPG2")
b <- ggplot(pcaHuh, aes(x = PC1, y = PC2, colour = timepoint, shape = status)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", varHuh[1], "% variance")) + 
  ylab(paste0("PC2: ", varHuh[2], "% variance")) + 
  coord_fixed() +
  ggtitle("HUH7")

library(gridExtra)
grid.arrange(a, b)

#### GENE CLUSTERING ####
###BY MSTRG
## HEPG2
topHepOne <- head(oneHep[order(oneHep$padj),], 50)
topHepThr <- head(thrHep[order(thrHep$padj),], 50)
topHepSev <- head(sevHep[order(sevHep$padj),], 50)

matHepOne <- assay(vstHep)[rownames(topHepOne), c("Lib_19", "Lib_20", "Lib_21", "Lib_22", "Lib_23", "Lib_24")]
matHepThr <- assay(vstHep)[rownames(topHepThr), c("Lib_25", "Lib_26", "Lib_27", "Lib_28", "Lib_29", "Lib_30")]
matHepSev <- assay(vstHep)[rownames(topHepSev), c("Lib_31", "Lib_32", "Lib_33", "Lib_34", "Lib_35", "Lib_36")]

matHepOne <- matHepOne - rowMeans(matHepOne)
matHepThr <- matHepThr - rowMeans(matHepThr)
matHepSev <- matHepSev - rowMeans(matHepSev)

anno_hep <- as.data.frame(colData(vstHep)[,c("status", "timepoint")])

pheatmap(matHepOne, annotation_col = anno_hep, main = "HepG2")
pheatmap(matHepThr, annotation_col = anno_hep, main = "HepG2")
pheatmap(matHepSev, annotation_col = anno_hep, main = "HepG2")

## HUH7
topHuhOne <- head(oneHuh[order(oneHuh$padj),], 50)
topHuhThr <- head(thrHuh[order(thrHuh$padj),], 50)
topHuhSev <- head(sevHuh[order(sevHuh$padj),], 50)

matHuhOne <- assay(vstHuh)[rownames(topHuhOne), c("Lib_1", "Lib_2", "Lib_3", "Lib_4", "Lib_5", "Lib_6")]
matHuhThr <- assay(vstHuh)[rownames(topHuhThr), c("Lib_7", "Lib_8", "Lib_9", "Lib_10", "Lib_11", "Lib_12")]
matHuhSev <- assay(vstHuh)[rownames(topHuhSev), c("Lib_13", "Lib_14", "Lib_15", "Lib_16", "Lib_17", "Lib_18")]

matHuhOne <- matHuhOne - rowMeans(matHuhOne)
matHuhThr <- matHuhThr - rowMeans(matHuhThr)
matHuhSev <- matHuhSev - rowMeans(matHuhSev)

anno_huh <- as.data.frame(colData(vstHuh)[,c("status", "timepoint")])

pheatmap(matHuhOne, annotation_col = anno_huh, main = "Huh7")
pheatmap(matHuhThr, annotation_col = anno_huh, main = "Huh7")
pheatmap(matHuhSev, annotation_col = anno_huh, main = "Huh7")

### BY GENE NAME or ENS ID
matList <- list("Hep One" = matHepOne, "Hep Three" = matHepThr, "Hep Sev" = matHepSev, "Huh One" = matHuhOne, "Huh Three" = matHuhThr, "Huh Seven" = matHuhSev)

genify <- function(x, type){
  x <- merge(x, ens, by.x = 0, by.y = "MSTRG")
  rownames(x) <- x[,type]
  x <- x[,c(-1,-8,-9,-10)]
}

matListId <- lapply(matList, genify, type = "gene_id")
matListName <- lapply(matList, genify, type = "gene_name")

for (x in grep("Hep", names(matListId))){
  pheatmap(matListId[[x]], annotation_col = anno_hep, main = names(matListId)[x])
  pheatmap(matListName[[x]], annotation_col = anno_hep, main = names(matListName)[x])
}

for (x in grep("Huh", names(matListId))){
  pheatmap(matListId[[x]], annotation_col = anno_huh, main = names(matListId)[x])
  pheatmap(matListName[[x]], annotation_col = anno_huh, main = names(matListName)[x])
}


### PLOT COUNTS ####
## Clean up the plots
## HEPG2
listHepOne <- as.list(rownames(topHepOne)[c(1:12)])
listHepThr <- as.list(rownames(topHepThr)[c(1:12)])
listHepSev  <- as.list(rownames(topHepSev)[c(1:12)])

dataHepOne <- lapply(listHepOne, plotCounts2, dds = dds_hepg2, intgroup = c("status", "timepoint", "group", "replicate"), returnData = T)
dataHepThr <- lapply(listHepThr, plotCounts2, dds = dds_hepg2, intgroup = c("status", "timepoint", "group", "replicate"), returnData = T)
dataHepSev <- lapply(listHepSev, plotCounts2, dds = dds_hepg2, intgroup = c("status", "timepoint", "group", "replicate"), returnData = T)

library(plyr)
dfHepOne <- ldply(dataHepOne, data.frame)
dfHepThr <- ldply(dataHepThr, data.frame)
dfHepSev <- ldply(dataHepSev, data.frame)

ggplot(dfHepOne, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap("gene", scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 1, HepG2")

ggplot(subset(dfHepOne, timepoint == 1), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 1, HepG2")

ggplot(dfHepThr, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 3, HepG2")

ggplot(subset(dfHepThr, timepoint == 3), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = replicate, shape = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 3, HepG2")

ggplot(dfHepSev, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 7, HepG2")

ggplot(subset(dfHepSev, timepoint == 7), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 7, HepG2")

## HUH7
listHuhOne <- as.list(rownames(topHuhOne)[c(1:12)])
listHuhThr <- as.list(rownames(topHuhThr)[c(1:12)])
listHuhSev  <- as.list(rownames(topHuhSev)[c(1:12)])

dataHuhOne <- lapply(listHuhOne, plotCounts2, dds = dds_huh7, intgroup = c("status", "timepoint", "group"), returnData = T)
dataHuhThr <- lapply(listHuhThr, plotCounts2, dds = dds_huh7, intgroup = c("status", "timepoint", "group"), returnData = T)
dataHuhSev <- lapply(listHuhSev, plotCounts2, dds = dds_huh7, intgroup = c("status", "timepoint", "group"), returnData = T)

dfHuhOne <- ldply(dataHuhOne, data.frame)
dfHuhThr <- ldply(dataHuhThr, data.frame)
dfHuhSev <- ldply(dataHuhSev, data.frame)

ggplot(dfHuhOne, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 1, Huh7")

ggplot(subset(dfHuhOne, timepoint == 1), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 1, Huh7")

ggplot(dfHuhThr, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 3, Huh7")

ggplot(subset(dfHuhThr, timepoint == 3), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 3, Huh7")

ggplot(dfHuhSev, aes(y = count, x = group)) + geom_point(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 7, Huh7")

ggplot(subset(dfHuhSev, timepoint == 7), aes(y = count, x = group)) + geom_jitter(size = 2, aes(colour = status)) + theme_bw() + facet_wrap(~gene, scales = "free_y") + theme(axis.text.x = element_text(angle = 60, hjust =1)) + ggtitle("Top 12 DE genes, Day 7, Huh7")

#### VOLCANO PLOT ####
a <- subset(thrHuh, !is.na(padj))
volc <- data.frame(logFC = a$log2FoldChange, negLogPval = -log10(a$padj))
maxval <- max(volc$logFC)
minval <- abs(min(volc$logFC))
lim <- round(max(c(minval, maxval)))
pval <- alpha
lfc <- 0.586
volc$sig <- (volc$logFC >= lfc | volc$logFC <= -lfc) & volc$negLogPval >= -log10(pval)

ggplot(volc, aes(x = logFC, y = negLogPval)) + 
  geom_point(aes(colour = sig)) +
  geom_hline(yintercept = -log10(pval), linetype = "dashed", size = 0.5, colour = "red") +
  geom_vline(xintercept = -lfc, linetype = "dashed", size = 0.5, colour = "blue") +
  geom_vline(xintercept = lfc, linetype = "dashed", size = 0.5, colour = "blue") +
  scale_x_continuous(limits = c(-lim, lim)) +
  theme_bw() +
  theme(text = element_text(size = 10)) + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("grey50", "red")) +
  theme(legend.position = "none") +
  xlab(expression(log[2]~fold~change)) +
  ylab(expression(-log[10](adjusted~p~value))) +
  annotate(geom = "text", label = nrow(subset(volc, sig == T & logFC >= lfc)), y = max(volc$negLogPval), x = lim/2) +
  annotate(geom = "text", label = nrow(subset(volc, sig == T & logFC <= -lfc)), y = max(volc$negLogPval), x = -lim/2) +
  annotate(geom = "text", label = nrow(subset(volc, sig == F)) - sum(is.na(volc$negLogPval)), y = max(volc$negLogPval), x = 0) +
  annotate(geom = "text", label = paste("p = ", pval), colour = "red", x = -lim, y = pval, hjust = 0, vjust = -1) +
  annotate(geom = "text", label = lfc, colour = "blue", x = lfc + 0.1, y = 0, hjust = "left", vjust = "top") +
  annotate(geom = "text", label = -lfc, colour = "blue", x = -lfc - 0.1, y = 0, hjust = "right", vjust = "top")





### ANNOTATE AND EXPORT ####
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

allResults <- list("oneHuh" = oneHuh, "thrHuh" = thrHuh, "sevHuh" = sevHuh, "oneHep" = oneHep, "thrHep" = thrHep, "sevHep" = sevHep)
allResultsDF <- lapply(allResults, data.frame)
allResultsGene <- lapply(allResultsDF, merge, by.x = 0, by.y = "MSTRG", y = ens[,c(1,2)])

sigListGene <- lapply(sigListDF, merge, y = ens[,c(1,2)], by.x = 0, by.y = "MSTRG")

stringify <- function(x){toString(x)}

sigListGeneAnno <- lapply(sigListGene, 
       function(x){
         x$symbol <- mapIds(org.Hs.eg.db,
                             key = x$gene_id,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first");
         x$entrez_id <- mapIds(org.Hs.eg.db,
                               key = x$gene_id,
                               column = "ENTREZID",
                               keytype ="ENSEMBL",
                               multiVals = "first");
         x$GO <- mapIds(org.Hs.eg.db,
                        key = x$gene_id,
                        column = "GO",
                        keytype = "ENSEMBL",
                        multiVals = stringify);
         x$EVIDENCE <- mapIds(org.Hs.eg.db,
                              key = x$gene_id,
                              column = "EVIDENCE",
                              keytype = "ENSEMBL",
                              multiVals = stringify);
        return(x)
       })



write_out <- function(x, y){
  z <- gsub(" ", "", x)
  loc <- paste("~/Dropbox/temp/", z, ".txt")
  write.table(y[[x]], file = loc, row.names = F, quote = F, sep = "\t")
}

sapply(names(sigListGeneAnno), write_out, y = sigListGeneAnno)

plotCounts2 <- function (dds, gene, intgroup = "condition", normalized = TRUE, 
                         transform = TRUE, main, xlab = "group", returnData = FALSE, 
                         replaced = FALSE, pc, ...) 
{
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(colData(dds)))) 
    stop("all variables in 'intgroup' must be columns of colData")
  stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
                                                             "factor"))))
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene, 
                                                                    ]
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
                              levels(colData(dds)[[intgroup[2]]]), function(x, 
                                                                            y) paste(x, y, sep = " : "))))
    droplevels(factor(apply(as.data.frame(colData(dds)[, 
                                                       intgroup, drop = FALSE]), 1, paste, collapse = " : "), 
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
                 1, paste, collapse = " : "))
  }
  data <- data.frame(count = cnts + pc, group = as.integer(group))
  logxy <- if (transform) 
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    return(data.frame(count = data$count, colData(dds)[intgroup], gene = main))
  plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count, 
       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
       xlab = xlab, ylab = ylab, main = main, ...)
  axis(1, at = seq_along(levels(group)), levels(group))
}

