library(data.table)
library(tidyverse)
library(patchwork)
library(bigmemory)
library(ggsci)
library(jcolors)

theme_set(theme_classic(base_size = 16))
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

dat$rs3_177339629 <- geno[which(map$SNP=="rs3_177340601"),]
dat$rs3_177339629[dat$rs3_177339629==0] <- "A"
dat$rs3_177339629[dat$rs3_177339629==2] <- "G"

dat$rs2_182382835 <- geno[which(map$SNP=="rs2_182382835"),]
dat$rs2_182382835[dat$rs2_182382835==0] <- "C"
dat$rs2_182382835[dat$rs2_182382835==2] <- "G"

dat$rs3_147531265 <- geno[which(map$SNP=="rs3_147531265"),]
dat$rs3_147531265[dat$rs3_147531265==0] <- "T"
dat$rs3_147531265[dat$rs3_147531265==2] <- "G"

dat$rs1_86064177 <- geno[which(map$SNP=="rs1_86064177"),]
dat$rs1_86064177[dat$rs1_86064177==0] <- "G"
dat$rs1_86064177[dat$rs1_86064177==2] <- "A"

dat$rs5_187859720 <- geno[which(map$SNP=="rs5_187859720"),]
dat$rs5_187859720[dat$rs5_187859720==0] <- "G"
dat$rs5_187859720[dat$rs5_187859720==2] <- "A"

rs2_182382835
axis.set <- b %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

############
### PsbS ###
############
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
  geom_segment(aes(x = 728674539, xend=728674539, y=2, yend=25.5), size=0.5) + 
  annotate("text", x=728674539, y=25, hjust=1.2, label=expression(italic(PSBS)), size=10) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + # points for RMIP values
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
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.3, 0.8), 
        legend.direction="horizontal", legend.text = element_text(angle = 45, hjust = 0.5)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))
  
g.psbs.box <- ggplot(dat, aes(rs3_177339629, Zm00001d042697, fill=rs3_177339629)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(PSBS)~(FPKM))) + 
  xlab("3:177,340,601") + 
  scale_fill_d3()

g.psbs + inset_element(g.psbs.zoom, left = 0.4, bottom = 0.4, right = 1, top = 1) + g.psbs.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/PsbS_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

############
### PSI3 ###
############

g.psi3 <- b %>%
  filter(-log10(Zm00001d042669.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042669.MLM))) +
  geom_segment(aes(x = 727490936, xend=728674539, y=2, yend=10), size=0.5) + 
  annotate("text", x=727490936, y=10, hjust=1.2, label=expression(italic(PSI3)), size=10) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/PSI3_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

############
### IRM1 ###
############
IRM1.data <- b %>%
  filter(CHROM==2 & POS > 181902891 & POS < 182859003)
IRM1.LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs2_182382835.ld")
IRM1.LD <- IRM1.LD[,c(6,7)]
colnames(IRM1.LD) <- c("SNP", "R2")
IRM1.data <- merge(IRM1.data[,c(1:5,8)], IRM1.LD, by="SNP")
IRM1.data <- IRM1.data[order(IRM1.data$R2, decreasing = F),]
rm(IRM1.LD)

g.irm1 <- b %>%
  filter(-log10(Zm00001d005657.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d005657.MLM))) +
  geom_segment(aes(x = 489304280, xend=489304280, y=2, yend=45), size=0.5) + 
  annotate("text", x=489304280, y=42, hjust=1.2, label=expression(italic(IRM1)), size=10) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.irm1.zoom <- ggplot(IRM1.data, aes(x = POS, y = -log10(Zm00001d005657.MLM), colour=R2)) +
  geom_vline(xintercept = 182382405) +
  geom_point(size=4) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD againt\n2:182,382,835") + 
  labs(x = "Chromosome 3",  y = expression(-log[10](italic(p)))) + 
  scale_x_continuous(labels = paste0(c("182,382,800", "182,383,400"), " kB"),
                   breaks = c(182382800, 182383400)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.75, 0.8), 
        legend.direction="horizontal", legend.text = element_text(angle = 45, hjust = 0.5)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))

g.irm1.box <- ggplot(dat, aes(rs2_182382835, Zm00001d005657, fill=rs2_182382835)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(IRM1)~(FPKM))) + 
  xlab("2:182,382,835") + 
  scale_fill_d3()

g.irm1 + inset_element(g.irm1.zoom, left = 0.3, bottom = 0.4, right = 1, top = 1) + g.irm1.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/IRM1_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

############
### TRX1 ###
############
TRX1.data <- b %>%
  filter(CHROM==3 & POS > 147030941 & POS < 148299517)
TRX1.LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_147531265_2.ld")
TRX1.LD <- TRX1.LD[,c(6,7)]
colnames(TRX1.LD) <- c("SNP", "R2")
TRX1.data <- merge(TRX1.data[,c(1:5,9)], TRX1.LD, by="SNP")
TRX1.data <- TRX1.data[order(TRX1.data$R2, decreasing = F),]
rm(TRX1.LD)

g.trx1 <- b %>%
  filter(-log10(Zm00001d042017.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042017.MLM))) +
  geom_segment(aes(x = 698866766, xend=698866766, y=2, yend=25), size=0.5) + 
  annotate("text", x=698866766, y=24, hjust=1.2, label=expression(italic(TRX1 - Y1)), size=10) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) +
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.trx1.zoom <- ggplot(TRX1.data, aes(x = POS, y = -log10(Zm00001d042017.MLM), colour=R2)) +
  geom_vline(xintercept = 147530941) +
  geom_point(size=4) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD againt\n3:147,531,265") + 
  labs(x = "Chromosome 3",  y = expression(-log[10](italic(p)))) + 
  scale_x_continuous(labels = paste0(c("147,300", "147,600", "147,900"), " kB"),
                     breaks = c(147300000, 147600000,147900000)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.22, 0.8), 
        legend.direction="horizontal", legend.text = element_text(angle = 45, hjust = 0.5)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))

g.trx1.box <- ggplot(dat, aes(rs3_147531265, Zm00001d042017, fill=rs3_147531265)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(TRX-Y1)~(FPKM))) + 
  xlab("3:147,531,265") + 
  scale_fill_d3()

g.trx1 + inset_element(g.trx1.zoom, left = 0.375, bottom = 0.35, right = 1, top = 1) + g.trx1.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/TRX1_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

#############
### ACHT3 ###
#############

g.acht3 <- b %>%
  filter(-log10(Zm00001d022518.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d022518.MLM))) +
  geom_segment(aes(x = 1610524859, xend=1610524859, y=2, yend=12), size=0.5) + 
  annotate("text", x=1610524859, y=12, hjust=1, label=expression(italic(ACHT3)), size=10) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/ACHT3_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

############
### PMI1 ###
############

PMI1.data <- b %>%
  filter(CHROM==1 & POS > 85914177 & POS < 86164177)
PMI1.LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs1_86064177.ld")
PMI1.LD <- PMI1.LD[,c(6,7)]
colnames(PMI1.LD) <- c("SNP", "R2")
PMI1.data <- merge(PMI1.data[,c(1:5,11)], PMI1.LD, by="SNP")
PMI1.data <- PMI1.data[order(PMI1.data$R2, decreasing = F),]
rm(PMI1.LD)

g.pmi1 <- b %>%
  filter(-log10(Zm00001d029761.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d029761.MLM))) +
  geom_segment(aes(x = 86064194, xend=86064194, y=2, yend=16), size=0.5) + 
  annotate("text", x=86064194, y=16, hjust=-0.25, label=expression(italic(PMI1)), size=10) +
  geom_point(aes(colour=as.character(CHROM)), size=3) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.pmi1.zoom <- ggplot(PMI1.data, aes(x = POS, y = -log10(Zm00001d029761.MLM), colour=R2)) +
  geom_vline(xintercept = 86064277) +
  geom_point(size=4) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD againt\n1:86,064,177") + 
  labs(x = "Chromosome 3",  y = expression(-log[10](italic(p)))) + 
  scale_x_continuous(labels = paste0(c("860,00", "860,50", "861,000"), " kB"),
                     breaks = c(86000000, 86050000,86100000)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.75, 0.8), 
        legend.direction="horizontal", legend.text = element_text(angle = 45, hjust = 0.5)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))

g.pmi1.box <- ggplot(dat, aes(rs1_86064177, Zm00001d029761, fill=rs1_86064177)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(PMI1)~(FPKM))) + 
  xlab("1:86,064,177") + 
  scale_fill_d3()

g.pmi1 + inset_element(g.pmi1.zoom, left = 0.3, bottom = 0.35, right = 1, top = 1) + g.pmi1.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/PMI1_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)


#############
### OEP37 ###
#############
OEP37.data <- b %>%
  filter(CHROM==5 & POS > 187812124 & POS < 187860118)
OEP37.LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs5_187859720.ld")
OEP37.LD <- OEP37.LD[,c(6,7)]
colnames(OEP37.LD) <- c("SNP", "R2")
OEP37.data <- merge(OEP37.data[,c(1:5,12)], OEP37.LD, by="SNP")
OEP37.data <- OEP37.data[order(OEP37.data$R2, decreasing = F),]
rm(OEP37.LD)

g.oep37 <- b %>%
  filter(-log10(Zm00001d017171.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d017171.MLM))) +
  geom_segment(aes(x = 1221622419, xend=1221622419, y=2, yend=45), size=0.5) + 
  annotate("text", x=1221622419, y=45, hjust=-0.25, label=expression(italic(OEP37)), size=10) +
  geom_point(aes(colour=as.character(CHROM)), size=3) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g.oep37.zoom <- ggplot(OEP37.data, aes(x = POS, y = -log10(Zm00001d017171.MLM), colour=R2)) +
  geom_vline(xintercept = 187855558) +
  geom_point(size=4) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD againt\n5:187,859,720") + 
  labs(x = "Chromosome 3",  y = expression(-log[10](italic(p)))) + 
  scale_x_continuous(labels = paste0(c("187,820", "187,840", "187,860"), " kB"),
                     breaks = c(187820000, 187840000,187860000)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position = c(0.55, 0.8), 
        legend.direction="horizontal", legend.text = element_text(angle = 45, hjust = 0.5)) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0),
         size = guide_legend(title.position="top", title.hjust = 0))

g.oep37.box <- ggplot(dat, aes(rs5_187859720, Zm00001d017171, fill=rs5_187859720)) + 
  geom_boxplot() + 
  theme(legend.position = "none") + 
  ylab(expression(italic(OEP37)~(FPKM))) + 
  xlab("5:187,859,720") + 
  scale_fill_d3()

g.oep37 + inset_element(g.oep37.zoom, left = 0.01, bottom = 0.35, right = 0.5, top = 1) + g.oep37.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/eQTL/OEP37_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)
