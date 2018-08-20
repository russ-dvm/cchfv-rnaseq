library(tidyverse)

virus <- read.table("~/cchf/8_viral_alignment/summaries_RK/summary.txt", h = T)
key <- read.table("~/cchf/ref_files/cchv.key", h = T)

virus <- merge(virus, key[,c(1,2,3,5)], by.x = "Library", by.y = "library")

p1 <- ggplot(virus, aes(x = timepoint, y = Reads)) + 
  geom_point(aes(shape = cell_line, colour = status), size = 3) + 
  theme_bw() + 
  scale_shape_manual(values = c(1,4)) +
  xlab("Day") +
  scale_x_continuous(breaks = c(1:7)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black")) +
  theme(text = element_text(size = 12))

p1


## A more detailed look at read depth/coverage of the viral genome
## The coverage files were generated with bedootls
coverage <- read.table("~/cchf/8_viral_alignment/bams_RK/all_cov.txt", sep = "\t", h = F)
colnames(coverage) <- c("chrom", "pos", "depth", "library")

coverage %>% group_by(library) %>%  group_by(chrom) # %>% summarise(mean = mean(depth), n = n()))
b <- coverage %>% group_by(chrom, library)
avgs <- b %>%  summarise(avg_cov = mean(depth))
avgs <- as.data.frame(avgs)
avgs$avg_cov <- as.integer(avgs$avg_cov)
avg_info <- merge(avgs, key[,c(1,2,3,5)], by.x = "library", by.y = "library") 
avg_info$seg <- NA
avg_info[avg_info$chrom == "KY484034",]$seg <- "L"
avg_info[avg_info$chrom == "KY484035",]$seg <- "M"
avg_info[avg_info$chrom == "KY484036",]$seg <- "S"


ggplot(avg_info, aes(x = seg, y = avg_cov, fill = status)) + 
  geom_boxplot() +
  ylab("Average per base depth of coverage") +
  xlab("CCHFV segment") +
  theme_classic() +
  theme(text = element_text(size = 12)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black"))


p3 <- ggplot(subset(avg_info, chrom == "KY484034"), aes(x = chrom, y = avg_cov, fill = status)) + geom_boxplot()

