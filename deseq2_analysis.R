library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)

#### Data import from prepDE.py ####
countData <- as.matrix(read.csv("~/cchf/6a_reestimate_for_deseq2/gene_count_matrix.csv", row.names="gene_id"))
colData <- read.table("~/cchf/ref_files/cchv.key", row.names=1, sep="\t", h=T)
colData$status <- relevel(colData$status, ref="uninfected")
colData$timepoint <- as.factor(colData$timepoint)
str(colData)


## "design = ~ <variable_1> + <variable_2>" == we want to test for the effect of variable_2 controlling for the effect of variable_1
## From the vignette: "As condition is the variable of interest, we put it at the end of the formula"
## Status:timepoint == status-specific differences over time
## Add cell_line:timepoint? cell_line specific differences over time?
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ status + timepoint + cell_line + status:timepoint)


## Many rows have a total count of <= 1. Remove them.
dds <- dds[ rowSums(counts(dds)) > 1, ]

#### Some prelim metrics ####
## Transform the data using the rlog. Blind set to false to ensure that the variables being investigated, status vs cell_line, will not contribute to variance-mean trend of the experiment.
rld <- rlog(dds, blind = F)
vsd <- vst(dds, blind = F)
# rld_2 <- rld
#### Find the euclidean mean-sample distances ####
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$status, rld$cell_line, rld$replicate, rld$timepoint, sep = " - " )
colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- viridis_pal()(1005)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#### PCA ####
## Modify the function slightly to incorporate an extra variable - in this case, deafults to cell line,
## so that we can pass that variable to ggplot2's shape aesthetic. Makes for a nicer plot - fewer colours.
plotPCA2<-function (object, intgroup = "condition", extragroup = "status", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  extragroup.df <- as.data.frame(colData(object)[, extragroup, drop =F])
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, name = colnames(object), extra = extragroup.df)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", shape = extragroup)) +
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    coord_fixed()
}

a<-plotPCA2(rld, intgroup = c("cell_line", "timepoint"), returnData = T)

ggplot(data = a, aes_string(x = "PC1", y = "PC2", color = "group", shape = "status")) +
  geom_point(size = 3) + 
  theme_bw() +
  xlab(paste0("PC1: ", round(attr(a, "percentVar")[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(attr(a, "percentVar")[2] * 100), "% variance")) + 
  coord_fixed()
plotPCA(vsd, intgroup = c("cell_line", "timepoint"))

#### Differential Expression ####
dd1 <- DESeq(dds, test = "LRT", reduced = ~ status + timepoint + cell_line)

## Alpha value == adjusted p value below the given FDR cutoff. Default = 0.1
## lfcthreshold = raises the log2 fold change threshold, so that we test for genes that show more substnaitial changes due to treatment. Equivalent to saying show significant effects for genes that change by a factor of 2^lfcthreshold. e.g. lfcthreshold = 1 looks for genes that either double or halve their expression.
res1 <- results(dd1, alpha = 0.1)
res1$symbol <- mcols(dd1)$symbol

p1 <- plotCounts(dd1, which.min(res1$padj), intgroup = c("timepoint", "status", "cell_line"), returnData = T)

ggplot(p1, aes(x = timepoint, y = count, colour = cell_line, group = status, shape = status)) + geom_point() + geom_smooth(se = F, method = "loess") + scale_y_log10() + ggtitle(row.names(dd1[which.min(res1$padj),]))

betas <- coef(dd1)
colnames(betas)

topGenes <- head(order(res1$padj), 10)
mat <- betas[topGenes, -c(1,2)]
thr <- 10
mat[mat < -thr] <- -thr
mat[mat > thr] <= thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col = F)
