library(tidyverse)
library(ggpubr)
library(gridExtra)

lib1 <- read.table("~/cchf/7_deseq/huh/Lib_1.counts.txt", h = T)
lib2 <- read.table("~/cchf/7_deseq/huh/Lib_2.counts.txt", h = T)
lib3 <- read.table("~/cchf/7_deseq/huh/Lib_3.counts.txt", h = T)
lib4 <- read.table("~/cchf/7_deseq/huh/Lib_4.counts.txt", h = T)
lib5 <- read.table("~/cchf/7_deseq/huh/Lib_5.counts.txt", h = T)
lib6 <- read.table("~/cchf/7_deseq/huh/Lib_6.counts.txt", h = T)
lib7 <- read.table("~/cchf/7_deseq/huh/Lib_7.counts.txt", h = T)
lib8 <- read.table("~/cchf/7_deseq/huh/Lib_8.counts.txt", h = T)
lib9 <- read.table("~/cchf/7_deseq/huh/Lib_9.counts.txt", h = T)
lib10 <- read.table("~/cchf/7_deseq/huh/Lib_10.counts.txt", h = T)
lib11 <- read.table("~/cchf/7_deseq/huh/Lib_11.counts.txt", h = T)
lib12 <- read.table("~/cchf/7_deseq/huh/Lib_12.counts.txt", h = T)
lib13 <- read.table("~/cchf/7_deseq/huh/Lib_13.counts.txt", h = T)
lib14 <- read.table("~/cchf/7_deseq/huh/Lib_14.counts.txt", h = T)
lib15 <- read.table("~/cchf/7_deseq/huh/Lib_15.counts.txt", h = T)
lib16 <- read.table("~/cchf/7_deseq/huh/Lib_16.counts.txt", h = T)
lib17 <- read.table("~/cchf/7_deseq/huh/Lib_17.counts.txt", h = T)
lib18 <- read.table("~/cchf/7_deseq/huh/Lib_18.counts.txt", h = T)

huhUn1 <- inner_join(lib1, lib2, by = "transcript", suffix = c(".lib1", ".lib2"))
huhUn1 <- inner_join(huhUn1, lib3, by = "transcript")

huhInf1 <- inner_join(lib4, lib5, by = "transcript", suffix = c(".lib4", ".lib5"))
huhInf1 <- inner_join(huhInf1, lib6, by = "transcript")

huhUn3 <- inner_join(lib7, lib8, by = "transcript", suffix = c(".lib7", ".lib8"))
huhUn3 <- inner_join(huhUn3, lib9, by = "transcript")

huhInf3 <- inner_join(lib10, lib11, by = "transcript", suffix = c(".lib10", ".lib11"))
huhInf3 <- inner_join(huhInf3, lib12, by = "transcript")

huhUn7 <- inner_join(lib13, lib14, by = "transcript", suffix = c(".lib13", ".lib14"))
huhUn7 <- inner_join(huhUn7, lib15, by = "transcript")

huhInf7 <- inner_join(lib16, lib17, by = "transcript", suffix = c(".lib16", ".lib17"))
huhInf7 <- inner_join(huhInf7, lib18, by = "transcript")

#### HUH DAY 1 UNINFECTED ####
huhUn1_1v2 <- ggplot(huhUn1, aes(x = FPKM.lib1, y = FPKM.lib2)) + geom_point() + stat_cor(method = "spearman") + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib1") + ylab("Lib2")
huhUn1_1v3 <- ggplot(huhUn1, aes(x = FPKM.lib1, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib1") + ylab("Lib3")
huhUn1_2v3 <- ggplot(huhUn1, aes(x = FPKM.lib2, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib2") + ylab("Lib3")

grid.arrange(huhUn1_1v2, huhUn1_1v3, huhUn1_2v3, nrow = 1, ncol = 3)

#### HUH DAY 3 UNINFECTED ####
huhUn3_7v8 <- ggplot(huhUn3, aes(x = FPKM.lib7, y = FPKM.lib8)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib7") + ylab("Lib8")
huhUn3_7v9 <- ggplot(huhUn3, aes(x = FPKM.lib7, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib7") + ylab("Lib9")
huhUn3_8v9 <- ggplot(huhUn3, aes(x = FPKM.lib8, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib8") + ylab("Lib9")

grid.arrange(huhUn3_7v8, huhUn3_7v9, huhUn3_8v9, nrow = 1, ncol = 3)

#### HUH DAY 7 UNINFECTED ####
huhUn7_13v14 <- ggplot(huhUn7, aes(x = FPKM.lib13, y = FPKM.lib14)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib13") + ylab("Lib14")
huhUn7_13v15 <- ggplot(huhUn7, aes(x = FPKM.lib13, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib13") + ylab("Lib15")
huhUn7_14v15 <- ggplot(huhUn7, aes(x = FPKM.lib14, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib14") + ylab("Lib15")

grid.arrange(huhUn7_13v14, huhUn7_13v15, huhUn7_14v15, nrow = 1, ncol = 3)

#### HUH DAY 1 INFECTED ####
huhInf1_4v5 <- ggplot(huhInf1, aes(x = FPKM.lib4, y = FPKM.lib5)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib4") + ylab("Lib5")
huhInf1_4v6 <- ggplot(huhInf1, aes(x = FPKM.lib4, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib4") + ylab("Lib6")
huhInf1_5v6 <- ggplot(huhInf1, aes(x = FPKM.lib5, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib5") + ylab("Lib6")

grid.arrange(huhInf1_4v5, huhInf1_4v6, huhInf1_5v6, nrow = 1, ncol = 3)

#### HUH DAY 3 INFECTED ####
huhInf3_10v11 <- ggplot(huhInf3, aes(x = FPKM.lib10, y = FPKM.lib11)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib10") + ylab("Lib11")
huhInf3_10v12 <- ggplot(huhInf3, aes(x = FPKM.lib10, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib10") + ylab("Lib12")
huhInf3_11v12 <- ggplot(huhInf3, aes(x = FPKM.lib11, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib11") + ylab("Lib12")

grid.arrange(huhInf3_10v11, huhInf3_10v12, huhInf3_11v12, nrow = 1, ncol = 3)

#### HUH DAY 7 INFECTED ####
huhInf7_16v17 <- ggplot(huhInf7, aes(x = FPKM.lib16, y = FPKM.lib17)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib16") + ylab("Lib17")
huhInf7_16v18 <- ggplot(huhInf7, aes(x = FPKM.lib16, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib16") + ylab("Lib18")
huhInf7_17v18 <- ggplot(huhInf7, aes(x = FPKM.lib17, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib17") + ylab("Lib18")

grid.arrange(huhInf7_16v17, huhInf7_16v18, huhInf7_17v18, nrow = 1, ncol = 3)


lib19 <- read.table("~/cchf/7_deseq/hep/Lib_19.counts.txt", h = T)
lib20 <- read.table("~/cchf/7_deseq/hep/Lib_20.counts.txt", h = T)
lib21 <- read.table("~/cchf/7_deseq/hep/Lib_21.counts.txt", h = T)
lib22 <- read.table("~/cchf/7_deseq/hep/Lib_22.counts.txt", h = T)
lib23 <- read.table("~/cchf/7_deseq/hep/Lib_23.counts.txt", h = T)
lib24 <- read.table("~/cchf/7_deseq/hep/Lib_24.counts.txt", h = T)
lib25 <- read.table("~/cchf/7_deseq/hep/Lib_25.counts.txt", h = T)
lib26 <- read.table("~/cchf/7_deseq/hep/Lib_26.counts.txt", h = T)
lib27 <- read.table("~/cchf/7_deseq/hep/Lib_27.counts.txt", h = T)
lib28 <- read.table("~/cchf/7_deseq/hep/Lib_28.counts.txt", h = T)
lib29 <- read.table("~/cchf/7_deseq/hep/Lib_29.counts.txt", h = T)
lib30 <- read.table("~/cchf/7_deseq/hep/Lib_30.counts.txt", h = T)
lib31 <- read.table("~/cchf/7_deseq/hep/Lib_31.counts.txt", h = T)
lib32 <- read.table("~/cchf/7_deseq/hep/Lib_32.counts.txt", h = T)
lib33 <- read.table("~/cchf/7_deseq/hep/Lib_33.counts.txt", h = T)
lib34 <- read.table("~/cchf/7_deseq/hep/Lib_34.counts.txt", h = T)
lib35 <- read.table("~/cchf/7_deseq/hep/Lib_35.counts.txt", h = T)
lib36 <- read.table("~/cchf/7_deseq/hep/Lib_36.counts.txt", h = T)

hepUn1 <- inner_join(lib19, lib20, by = "transcript", suffix = c(".lib19", ".lib20"))
hepUn1 <- inner_join(hepUn1, lib21, by = "transcript")

hepInf1 <- inner_join(lib22, lib23, by = "transcript", suffix = c(".lib22", ".lib23"))
hepInf1 <- inner_join(hepInf1, lib24, by = "transcript")

hepUn3 <- inner_join(lib25, lib26, by = "transcript", suffix = c(".lib25", ".lib26"))
hepUn3 <- inner_join(hepUn3, lib27, by = "transcript")

hepInf3 <- inner_join(lib28, lib29, by = "transcript", suffix = c(".lib28", ".lib29"))
hepInf3 <- inner_join(hepInf3, lib30, by = "transcript")

hepUn7 <- inner_join(lib31, lib32, by = "transcript", suffix = c(".lib31", ".lib32"))
hepUn7 <- inner_join(hepUn7, lib33, by = "transcript")

hepInf7 <- inner_join(lib34, lib35, by = "transcript", suffix = c(".lib34", ".lib35"))
hepInf7 <- inner_join(hepInf7, lib36, by = "transcript")

#### hep DAY 1 UNINFECTED ####
hepUn1_1v2 <- ggplot(hepUn1, aes(x = FPKM.lib19, y = FPKM.lib20)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib19") + ylab("Lib20")
hepUn1_1v3 <- ggplot(hepUn1, aes(x = FPKM.lib19, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib19") + ylab("Lib21")
hepUn1_2v3 <- ggplot(hepUn1, aes(x = FPKM.lib20, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib20") + ylab("Lib21")

grid.arrange(hepUn1_1v2, hepUn1_1v3, hepUn1_2v3, nrow = 1, ncol = 3)

#### hep DAY 3 UNINFECTED ####
hepUn3_7v8 <- ggplot(hepUn3, aes(x = FPKM.lib25, y = FPKM.lib26)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib25") + ylab("Lib26")
hepUn3_7v9 <- ggplot(hepUn3, aes(x = FPKM.lib25, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib25") + ylab("Lib27")
hepUn3_8v9 <- ggplot(hepUn3, aes(x = FPKM.lib26, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib26") + ylab("Lib27")

grid.arrange(hepUn3_7v8, hepUn3_7v9, hepUn3_8v9, nrow = 1, ncol = 3)

#### hep DAY 7 UNINFECTED ####
hepUn7_13v14 <- ggplot(hepUn7, aes(x = FPKM.lib31, y = FPKM.lib32)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib31") + ylab("Lib32")
hepUn7_13v15 <- ggplot(hepUn7, aes(x = FPKM.lib31, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib31") + ylab("Lib33")
hepUn7_14v15 <- ggplot(hepUn7, aes(x = FPKM.lib32, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib32") + ylab("Lib33")

grid.arrange(hepUn7_13v14, hepUn7_13v15, hepUn7_14v15, nrow = 1, ncol = 3)

#### hep DAY 1 INFECTED ####
hepInf1_4v5 <- ggplot(hepInf1, aes(x = FPKM.lib22, y = FPKM.lib23)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib22") + ylab("Lib23")
hepInf1_4v6 <- ggplot(hepInf1, aes(x = FPKM.lib22, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib22") + ylab("Lib24")
hepInf1_5v6 <- ggplot(hepInf1, aes(x = FPKM.lib23, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib23") + ylab("Lib24")

grid.arrange(hepInf1_4v5, hepInf1_4v6, hepInf1_5v6, nrow = 1, ncol = 3)

#### hep DAY 3 INFECTED ####
hepInf3_10v11 <- ggplot(hepInf3, aes(x = FPKM.lib28, y = FPKM.lib29)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib28") + ylab("Lib29")
hepInf3_10v12 <- ggplot(hepInf3, aes(x = FPKM.lib28, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib28") + ylab("Lib30")
hepInf3_11v12 <- ggplot(hepInf3, aes(x = FPKM.lib29, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib29") + ylab("Lib30")

grid.arrange(hepInf3_10v11, hepInf3_10v12, hepInf3_11v12, nrow = 1, ncol = 3)

#### hep DAY 7 INFECTED ####
hepInf7_16v17 <- ggplot(hepInf7, aes(x = FPKM.lib34, y = FPKM.lib35)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib34") + ylab("Lib35")
hepInf7_16v18 <- ggplot(hepInf7, aes(x = FPKM.lib34, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib34") + ylab("Lib36")
hepInf7_17v18 <- ggplot(hepInf7, aes(x = FPKM.lib35, y = FPKM)) + geom_point() + stat_cor() + geom_smooth(method = "lm", se = F) + theme_bw() + theme(panel.grid = element_blank()) + xlab("Lib35") + ylab("Lib36")

grid.arrange(hepInf7_16v17, hepInf7_16v18, hepInf7_17v18, nrow = 1, ncol = 3)
