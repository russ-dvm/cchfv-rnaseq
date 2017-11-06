library(tidyverse)
library(ggpubr)


hep <- read.csv("~/Dropbox/temp/gene_count_matrix.csv", h = T)
huh <- read.csv("~/Dropbox/temp/huh_gene_counts.csv", h = T)
a <- cor.test(hep$Lib_22, hep$Lib_23, method = "pearson")
shapiro.test(hep$Lib_22)

for (i in 2:ncol(hep)){
  
  print(cor(hep[,i], hep[,i+1], method = "spearman"))
  print(colnames(hep)[i])
  
}

for (i in 2:ncol(huh)){
  
  print(cor(huh[,i], huh[,i+1], method = "spearman"))
  print(colnames(huh)[i])
  
}
