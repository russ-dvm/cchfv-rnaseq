library(tidyverse)

## KEY
key <- read.table("~/cchf/ref_files/cchv_key.txt", h = T, sep = "\t")


## READS ALIGNED TO THE HUMAN GENOME
aligned <- read.table("~/cchf/3_aligned/complete_summary.txt", h = F, sep = "\t")
colnames(aligned) <- c("Aligner", "Library", "Metric", "Value")


aligned <- merge(aligned, key[,c(1,2,3,5)], by.x = "Library", by.y = "library")


ggplot(subset(aligned, aligned$Metric == "rate"), aes(x = cell_line, y = Value, fill = status)) + 
  geom_boxplot() 

tot <- subset(aligned, aligned$Metric == "total")
sum(subset(tot, tot$cell_line == "HepG2" & tot$status == "infected")$Value)

rate <- subset(aligned, aligned$Metric == "rate")
rate
summary(rate)

## READS ALIGNED TO THE CCHFV GENOME
source("~/cchf/cchfv-rnaseq/virus.R")
p1

sum(subset(virus, cell_line == "HepG2" & timepoint == 1 & status == "infected")$Reads)
