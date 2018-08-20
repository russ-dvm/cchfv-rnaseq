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

#huh7
huh7 <- as.matrix(read.csv("~/cchf/tmp/huh7/gene_count_matrix.csv", row.names="gene_id"))

colData_huh7 <- read.table("~/cchf/tmp/huh7/cchv.key", row.names=1, sep="\t", h=T)
colData_huh7$status <- relevel(colData_huh7$status, ref="uninfected")
colData_huh7$timepoint <- as.factor(colData_huh7$timepoint)
# colData_huh7$group <- as.factor(paste(colData_huh7$status, colData_huh7$timepoint, sep="_"))
colData_huh7 <- colData_huh7[-1]
str(colData_huh7)


## "design = ~ <variable_1> + <variable_2>" == test for the effect of variable_2 controlling for the effect of variable_1
## From the vignette: "As condition is the variable of interest, we put it at the end of the formula"
## I think we need an interaction term for this analysis
## Interaction terms are used in regression to analyze how two independent variablese might interact to change the response of the Y variable. 
## In other words: time affects the expression of gene Y, and infection changes the expression of gene Y.
## But, is the effect of infection on the expression of gene Y influenced by time? 
## We need to look at the interaction between time and infection.
## Remember: this is just a fancy Y = a + bX equation ;)

dds_huh7 <- DESeqDataSetFromMatrix(countData = huh7,
                                   colData = colData_huh7,
                                   design = ~ status + timepoint + status:timepoint)

## Many rows have a total count of <= 1. Remove them.
dds_huh7 <- dds_huh7[ rowSums(counts(dds_huh7)) > 1, ]

#### Some prelim metrics ####
## Transform the data using the rlog. Blind set to false to ensure that the variables being investigated, status vs cell_line, will not contribute to variance-mean trend of the experiment.
rld_huh7 <- rlog(dds_huh7, blind = F)

# vsd <- vst(dds, blind = F)
# rld_2 <- rld
#### Find the euclidean mean-sample distances ####
sampleDists <- dist(t(assay(rld_huh7)))

sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste( rld_huh7$status, rld_huh7$cell_line, rld_huh7$replicate, rld_huh7$timepoint, sep = " - " )
colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- viridis_pal()(1005)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main="")

#### PCA ####
a<-plotPCA2(rld_huh7, intgroup = c("timepoint"), returnData = T)

ggplot(data = a, aes_string(x = "PC1", y = "PC2", color = "timepoint", shape = "status")) +
  geom_point(size = 3) + 
  theme_bw() +
  xlab(paste0("PC1: ", round(attr(a, "percentVar")[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(attr(a, "percentVar")[2] * 100), "% variance")) + 
  coord_fixed() +
  ggtitle("Huh7")


#### Differential Expression ####
## Want to use likelihood ratio test for the multifactorial design (time and infection status)
## Requires a "reduced" parameter. Genes with a small p value at this point are those which at one or more time points (after time 0) showed a status-specific effect. This EXCLUDES genes that move up/down over time in both infected/non-infected.
## Thus the ratio is our default design (specified above) vs the null hypothesis that 
dds_huh7 <- DESeq(dds_huh7, test = "LRT", reduced = ~ status + timepoint)
dds_huh7 <- DESeq(dds_huh7, betaPrior = F)
resultsNames(dds_huh7)

## resultsNames provides information on different tests that were conducted. In this case, it looks like this:
## [1] "Intercept"                     "status_infected_vs_uninfected" "timepoint_3_vs_1"              "timepoint_7_vs_1"            [5] "statusinfected.timepoint3"     "statusinfected.timepoint7" 

## This model has interaction terms, which can be identified as the terms that contain two variable names: in this case, status and time. The term "statusinfected.timepoint3" is referring to a test for if the infected vs. uninfected cells have a gene expression fold change that is different at timepoint 3 than at timepoint 1. As another example, "timepoint_3_vs_1" is a test for whether there is a different log-fold change in the uninfected cells at Day 3 vs Day 1. 


## Alpha value == adjusted p value below the given FDR cutoff. Default = 0.1
## lfcthreshold = raises the log2 fold change threshold, so that we test for genes that show more substnaitial changes due to treatment. Equivalent to saying show significant effects for genes that change by a factor of 2^lfcthreshold. e.g. lfcthreshold = 1 looks for genes that either double or halve their expression.
res_huh7 <- results(dds_huh7, alpha = 0.05)

## Ordered by pvalue
ordered_res <- res_huh7[order(res_huh7$padj),]

## Counts plot
top_gene <- plotCounts2(dds_huh7, "MSTRG.19218", intgroup = c("status", "timepoint"), returnData = T)
top_gene <- rbind(top_gene, plotCounts2(dds_huh7, "MSTRG.19218", intgroup = c("status", "timepoint"), returnData = T ))

ggplot(top_gene, aes(x = as.numeric(timepoint), y = count, colour = status, group = status)) + geom_point() + geom_smooth(se=F, method = "loess") + theme_classic() + facet_grid(.~gene, scales = "free_y") + scale_y_log10() 

## Heatmap test
betas <- coef(dds_huh7)
colnames(betas)
top <- head(order(res_huh7$padj), 50)
mat <- betas[top, -c(1)]
thr <- 5
mat[mat < -thr] <- -thr
mat[mat > -thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col = F)

#### Functions ####
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

## Modified to return gene name when returnData = T
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
