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

