library(ggplot2)

virus <- read.table("~/cchf/8_viral_alignment/summaries_RK/summary.txt", h = T)
key <- read.table("~/cchf/ref_files/cchv.key", h = T)

virus <- merge(virus, key[,c(1,2,3,5)], by.x = "Library", by.y = "library")

ggplot(virus, aes(x = timepoint, y = Reads)) + 
  geom_point(aes(shape = cell_line, colour = status), size = 3) + 
  theme_bw() + 
  scale_shape_manual(values = c(1,4)) +
  ggtitle("Number of reads aligned to the CCHFV genome in two cell lines") +
  xlab("Day") +
  scale_x_continuous(breaks = c(1:7))
