library(tidyverse)
library(ggpubr)
library(plyr)


hep <- read.csv("~/Dropbox/temp/gene_count_matrix.csv", h = T)
huh <- read.csv("~/Dropbox/temp/huh_gene_counts.csv", h = T)


cor_fnc <- function(x, offset){ 
  cor_list <- list()
    for (i in 1:(ncol(x)-1)){
    j <- i
    i <- i + offset

  if (i %% 3 == 0){
    liba <- paste("^Lib_", i-2, "$", sep="")
    libb <- paste("^Lib_", i-1, "$", sep="")
    libc <- paste("^Lib_", i, "$", sep="")
    
    colLiba <- grep(liba, colnames(x))
    colLibb <- grep(libb, colnames(x))
    colLibc <- grep(libc, colnames(x))
    
    comp1 <- paste(liba, libb, sep="-")
    comp1 <- gsub("\\^", "", comp1)
    comp1 <- gsub("\\$", "", comp1)
    cor_ab <- cor(x[,colLiba], x[,colLibb], method = "spearman")
    
    comp2 <- paste(liba, libc, sep="-")
    comp2 <- gsub("\\^", "", comp2)
    comp2 <- gsub("\\$", "", comp2)
    cor_ac <- cor(x[,colLiba], x[,colLibc], method = "spearman")

    comp3 <- paste(libb, libc, sep="-")
    comp3 <- gsub("\\^", "", comp3)
    comp3 <- gsub("\\$", "", comp3)
    cor_bc <- cor(x[,colLibb], x[,colLibc], method = "spearman")

    cor_list[[j-2]] <- c(comp1, cor_ab)
    cor_list[[j-1]] <- c(comp2, cor_ac)
    cor_list[[j]] <- c(comp3, cor_bc)
  }
  }
  return(cor_list)
}

hepList <- cor_fnc(hep, 18)
hepDf <- ldply(hepList)
hepDfRound <- hepDf 
hepDfRound$V2 <- signif(as.numeric(hepDfRound$V2), digits = 3)

huhList <- cor_fnc(huh, 0)
huhDf <- ldply(huhList)
huhDfRound <- huhDf
huhDfRound$V2 <- signif(as.numeric(huhDfRound$V2), digits = 3)
