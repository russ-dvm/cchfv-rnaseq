library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(genefilter)


#### Data import from prepDE.py ####
# countData <- as.matrix(read.csv("~/cchf/tmp/gene_count_matrix_c0.csv", row.names="gene_id"))
# countData1 <- as.matrix(read.csv("~/cchf/tmp/gene_count_matrix_c1.csv", row.names="gene_id"))
# countData10 <- as.matrix(read.csv("~/cchf/tmp/gene_count_matrix_c10.csv", row.names="gene_id"))

#hepg2
hepg2 <- as.matrix(read.csv("~/cchf/tmp/hepg2/gene_count_matrix.csv", row.names="gene_id"))
colData_hepg2 <- read.table("~/cchf/tmp/hepg2/cchv.key", row.names=1, sep="\t", h=T)
colData_hepg2$status <- relevel(colData_hepg2$status, ref="uninfected")
colData_hepg2$timepoint <- as.factor(colData_hepg2$timepoint)
str(colData_hepg2)

## "design = ~ <variable_1> + <variable_2>" == we want to test for the effect of variable_2 controlling for the effect of variable_1
## From the vignette: "As condition is the variable of interest, we put it at the end of the formula"
## Status:timepoint == status-specific differences over time
## Add cell_line:timepoint? cell_line specific differences over time?
dds_hepg2 <- DESeqDataSetFromMatrix(countData = hepg2,
                              colData = colData_hepg2,
                              design = ~ status + timepoint + status:timepoint)

## Many rows have a total count of <= 1. Remove them.
dds_hepg2 <- dds_hepg2[ rowSums(counts(dds_hepg2)) > 6, ]

#### Some prelim metrics ####
## Transform the data using the rlog. Blind set to false to ensure that the variables being investigated, status vs cell_line, will not contribute to variance-mean trend of the experiment.
rld_hepg2 <- rlog(dds_hepg2, blind = F)

# vsd <- vst(dds, blind = F)
# rld_2 <- rld
#### Find the euclidean mean-sample distances ####
sampleDists <- dist(t(assay(rld_hepg2)))

sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste( rld$status, rld$cell_line, rld$replicate, rld$timepoint, sep = " - " )
colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- viridis_pal()(1005)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main="")

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

a<-plotPCA2(rld_hepg2, intgroup = c("timepoint"), returnData = T)

ggplot(a, aes(x = PC1, y = PC2, color = timepoint, shape = status)) +
  geom_point(size = 3) + 
  theme_bw() +
  xlab(paste0("PC1: ", round(attr(a, "percentVar")[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(attr(a, "percentVar")[2] * 100), "% variance")) + 
  coord_fixed() +
  ggtitle("HepG2")


#### Differential Expression ####
dds_hepg2 <- DESeq(dds_hepg2, test = "LRT", reduced = ~ status + timepoint)

## Alpha value == adjusted p value below the given FDR cutoff. Default = 0.1
## lfcthreshold = raises the log2 fold change threshold, so that we test for genes that show more substnaitial changes due to treatment. Equivalent to saying show significant effects for genes that change by a factor of 2^lfcthreshold. e.g. lfcthreshold = 1 looks for genes that either double or halve their expression. Default = 0
res_hepg2 <- results(dds_hepg2, alpha = 0.05, lfcthreshold = 0)
res_hepg2
summary(res_hepg2)

## MA Plot
res1_hepg2 <- results(dds_hepg2, lfcThreshold=2)
res2_hepg2 <- lfcShrink(dds_hepg2, contrast=c("status","infected","uninfected"), res=res1_hepg2)
plotMA(res2_hepg2, ylim = c(-2, 2), main="HepG2")

## Significant only
sig_hepg2 <- subset(res_hepg2, padj < 0.05)
table(sig_hepg2$padj < 0.05)

## Down and upregulated
# downreg <- head(resSig[ order(resSig$log2FoldChange), ])
upreg_hepg2 <- head(res_hepg2[ order(res_hepg2$padj, decreasing = TRUE), ], 30)


##Heatmap
upreg_top <- head(upreg_hepg2, 20)
topVarGenes <- head(order(rowVars(assay(rld_hepg2)), decreasing = TRUE), 50)
mat  <- assay(rld_hepg2)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld_hepg2)[, c("timepoint","status")])
pheatmap(mat, annotation_col = anno, main="HepG2 time")
