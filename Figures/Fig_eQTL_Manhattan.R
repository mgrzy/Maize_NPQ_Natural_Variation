library(data.table)
library(tidyverse)
library(patchwork)
library(bigmemory)
library(ggsci)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

dat <- read.csv("../Data/Maize_NPQ_Natural_Variation/data/work/exp_eqt_rnav4.csv")
geno <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.desc")
map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = F)

a <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/eQTL_rnav4.csv.gz", data.table = F)

nCHR <- length(unique(a$CHROM))
a$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(a$CHROM))){
  nbp[i] <- max(a[a$CHROM == i,]$POS)
  a[a$CHROM == i,"BPcum"] <- a[a$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

b <- a[,c(1,2,3,4,5,8,11,14,17,20,23,26,27)]

dat$rs3_177339629 <- geno[which(map$SNP=="rs3_177339629"),]
dat$rs3_177339629[dat$rs3_177339629==0] <- "C"
dat$rs3_177339629[dat$rs3_177339629==2] <- "T"

axis.set <- b %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

PsbS.data <- b %>%
  filter(CHROM==3 & POS > 175652936 & POS < 178038945)
PsbS.LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_177340601.ld")
PsbS.LD <- PsbS.LD[,c(6,7)]
colnames(PsbS.LD) <- c("SNP", "R2")
PsbS.data <- merge(PsbS.data[,c(1:6)], PsbS.LD, by="SNP")
PsbS.data <- PsbS.data[order(PsbS.data$R2, decreasing = F),]
rm(PsbS.LD)

g.psbs <- b %>%
  filter(-log10(Zm00001d042697.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042697.MLM))) +
  geom_segment(aes(x = 728674539, xend=728674539, y=2, yend=25.5)) + 
  annotate("text", x=728674539, y=25, hjust=1.2, label=expression(italic(PSBS))) + 
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none") 

g.psbs.zoom <- ggplot(PsbS.data, aes(x = POS, y = -log10(Zm00001d042697.MLM), colour=R2)) +
  geom_vline(xintercept = 177338945) +
  geom_point(size=4) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD againt\n3:177,340,601") + 
  labs(x = "Chromosome 3",  y = expression(-log[10](italic(p)))) + 
  scale_x_continuous(labels = paste0(c("176", "177", "178"), " MB"),
                     breaks = c(176000000, 177000000,178000000), limits = c(175652936, 178038945)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.2, 0.8)) 
  
g.psbs + inset_element(g.psbs.zoom, left = 0.4, bottom = 0.5, right = 1, top = 1)

g.psbs.box <- ggplot(dat, aes(rs3_177339629, Zm00001d042697, fill=rs3_177339629)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(PSBS)~(FPKM))) + 
  xlab("3:177,339,629") + 
  scale_fill_d3()

g.psbs + inset_element(g.psbs.zoom, left = 0.4, bottom = 0.5, right = 1, top = 1) + g.psbs.box + plot_layout(widths = c(4,1))


g.psi3 <- b %>%
  filter(-log10(Zm00001d042669.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042669.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.irm1 <- b %>%
  filter(-log10(Zm00001d005657.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d005657.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.trx1 <- b %>%
  filter(-log10(Zm00001d042017.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042017.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.acht3 <- b %>%
  filter(-log10(Zm00001d022518.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d022518.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.pmi1 <- b %>%
  filter(-log10(Zm00001d029761.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d029761.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.oep37 <- b %>%
  filter(-log10(Zm00001d017171.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d017171.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.psbs + g.psi3 + g.irm1 + g.trx1 + g.acht3 + g.pmi1 + g.oep37

ggsave(filename = "../Data/Maize_NPQ_Natural_Variation/figures/eQTL.png", device = "png", width = 12, height = 6)

