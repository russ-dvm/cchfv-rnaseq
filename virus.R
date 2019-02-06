library(tidyverse)

virus <- read.table("~/cchf/8_viral_alignment/summaries_RK/summary.txt", h = T)
key <- read.table("~/cchf/ref_files/cchv.key", h = T)

virus <- merge(virus, key[,c(1,2,3,5)], by.x = "Library", by.y = "library")
colnames(virus)[3] <- "Cell line"
colnames(virus)[4] <- "Status"

##omit day 7 as pe RK
virus_no7 <- subset(virus, timepoint != 7)

p1 <- ggplot(virus_no7, aes(x = timepoint, y = Reads)) + 
  geom_point(aes(shape = `Cell line`, colour = Status), size = 3) +
  theme_classic() + 
  scale_shape_manual(values = c(1,4)) +
  xlab("Day") +
  ylab("Number of reads") +
  scale_x_continuous(breaks = c(1:7)) +
  # theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black")) +
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
avg_info$timepoint <- paste("Day", avg_info$timepoint)
colnames(avg_info)[5] <- "Status"
##omit day 7 as per RK
avg_info_no7 <- subset(avg_info, timepoint != "Day 7")

ggplot(avg_info_no7, aes(x = seg, y = avg_cov, fill = Status)) + 
  geom_boxplot() +
  ylab("Average depth of coverage") +
  xlab("CCHFV segment") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  # theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black")) +
  facet_grid(cell_line ~ timepoint, scales = "free_y")

p3 <- ggplot(subset(avg_info, chrom == "KY484034"), aes(x = chrom, y = avg_cov, fill = Status)) + geom_boxplot()

##Export table, easier to summarize using a pivot table in excel
write.table(file = "~/Dropbox/temp/avg_info_no7.txt", avg_info_no7, row.names = F, quote = F, sep = "\t")