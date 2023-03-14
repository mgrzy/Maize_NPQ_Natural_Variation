library(data.table)
library(tidyverse)
library(patchwork)
library(jcolors)

# Set them to classic and axis text to black rather then grey
theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

# Load data

gff <- fread("../../BigData/Maize/v4/GFF/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3.gz",skip = 5, fill = T)
gff2 <- fread("../../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz",skip = 5, fill = T)

###

RMIP_2020 <- fread("../Data/Maize_NPQ_Natural_Variation/data/figures/Fig1c_RMIP.csv")
RMIP_2021 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/RMIPclean2021.csv")

###

RMIP_2020 <- RMIP_2020[RMIP>1,]
RMIP_2020$trait <- factor(RMIP_2020$trait, levels = c("Induction", "Steady-state", "Relaxation", "Relaxed-stage", 
                                                      "PSII", "PSIIrate"))
RMIP_2021$trait <- factor(RMIP_2021$trait, levels = c("Induction", "Steady-state", "Relaxation", "Relaxed-stage", 
                                                      "PSII", "PSIIrate"))

greeks <- list("Induction", "Steady-state", "Relaxation", "Relaxed-stage", bquote(phi[PSII]), bquote(phi[PSII_2]))

### PSI-PsbS
psi.psbs <- list()
psi.psbs$rmip2020 <- RMIP_2020[CHROM==3  & POS > 175652936 & POS < 178038945,]
psi.psbs$rmip2021 <- RMIP_2021[CHROM==3  & POS > 175652936 & POS < 178038945,]
psi.psbs$gen <- gff[V1=="Chr3" & V4 > 175652936 & V5 < 178038945]
psi.psbs$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_176152936.ld")
psi.psbs$psbs <- psi.psbs$gen[grep("Zm00001d042697",psi.psbs$gen$V9)]
psi.psbs$psi <- psi.psbs$gen[grep("Zm00001d042669",psi.psbs$gen$V9)]
psi.psbs$psbs_v5 <- gff2[grep("Zm00001eb146510", gff2$V9)]
psi.psbs$psi_v5 <- gff2[grep("Zm00001eb146270", gff2$V9)]

g.psi.rmip2020 <- ggplot(psi.psbs$rmip2020, aes(POS, RMIP)) + 
  geom_vline(xintercept = c(176152936, 177340601), linetype=2, colour="red") + 
  geom_point(size=4, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  xlim(175652936, 178038945) + ylab("RMIP\n2020") + 
  scale_color_manual(values = c("yellow", "orange", "darkgreen", "green",
                                        "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) # colour values for each trait

g.psi.rmip2021 <- ggplot(psi.psbs$rmip2021, aes(POS, RMIP)) + 
  geom_vline(xintercept = c(176152936, 177340601), linetype=2, colour="red") + 
  geom_point(size=4, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = c(0.9,0.7)) +
  xlim(175652936, 178038945) + ylab("RMIP\n2021") + 
  scale_color_manual(values = c("yellow", "orange", "darkgreen", "green",
                                        "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) # colour values for each trait

g.psi.LD <- ggplot(psi.psbs$LD, aes(BP_B, R2)) + 
  geom_vline(xintercept = c(176152936, 177340601), linetype=2, colour="red") + 
  geom_point(size=4, alpha=0.9, colour="dodgerblue4") + 
  geom_point(data=psi.psbs$LD[psi.psbs$LD$BP_B==177340601,], aes(BP_B, R2), colour="red", size=4, alpha=0.9) + 
  geom_point(data=psi.psbs$LD[psi.psbs$LD$BP_B==176152936,], aes(BP_B, R2), colour="red", size=4, alpha=0.9) + 
  ylab(expression(LD~(r^2))) + 
  xlab("Chromosome 3") + 
  theme() + 
  scale_x_continuous(labels = paste0(c("176,000", "176,500", "177,000", "177,500"), " kB"),
                     breaks = c(176000000, 176500000, 177000000,177500000), limits = c(175652936, 178038945))

g.psipsbs.gene <- ggplot() + 
  geom_point(data=psi.psbs$gen[psi.psbs$gen$V7=="+",], aes(x=V4, y=0.5), shape="\u25BA", size=5) +
  geom_point(data=psi.psbs$gen[psi.psbs$gen$V7=="-",], aes(x=V4, y=-0.5), shape="\u25C4", size=5) +
  geom_point(data=psi.psbs$psbs, aes(x=V4, y=-0.5), shape="\u25C4", size=5, colour="red") +
  geom_point(data=psi.psbs$psi, aes(x=V4, y=0.5), shape="\u25BA", size=5, colour="red") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ylim(-2, 2) + 
  xlim(c(175652936, 178038945))

g.psbs <- ggplot() + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.5, "cm"), ends = "first")) + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=179565678, y=-0.5, label=expression(Zm00001d042697~"("*italic(PSBS)*")")) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.psi <- ggplot() + 
  geom_segment(data=psi.psbs$psi_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.5, "cm"))) + 
  geom_segment(data=psi.psbs$psi_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=178410159, y=-0.5, label=expression(Zm00001d042669~"("*italic(PSI3)*")")) + 
  geom_segment(aes(y=0.5, yend=0, x=178405643, xend=178405643), size=1, colour="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)


g.12 <- (g.psi.rmip2020 / g.psi.rmip2021 / g.psi.LD / g.psipsbs.gene) / (g.psi + g.psbs) + 
  plot_layout(heights = c(2,2,2,2,1.5))

ggsave(plot = g.12, filename = "../Data/Maize_NPQ_Natural_Variation/figures/Fig2/PSI3_PsbS.svg", device = "svg", 
       units = "mm", width = 75, height = 40, scale = 4)
#########################################################################
#########################################################################

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
  geom_point(data=PsbS.data[PsbS.data$POS==176152936,], aes(POS, -log10(Zm00001d042697.MLM)), colour="red", size=4, alpha=0.9) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD relative\n3:177,340,601") + 
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
  ylab(expression(italic(PSBS)~Expression~(FPKM))) + 
  xlab("3:177,340,601") + 
  scale_fill_d3()

g.psbs + inset_element(g.psbs.zoom, left = 0.4, bottom = 0.4, right = 1, top = 1) + g.psbs.box + plot_layout(widths = c(4,1))

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/Fig2/PsbS_eQTL.svg", units = "mm", width = 75, height = 40, scale = 4)

##############################################################
##############################################################
library(tidyverse)
library(bigmemory)
library(data.table)
library(patchwork)
library(ggsci)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

geno <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.desc")
map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map")

genes <- read.csv("../npq/data/meta/GenesToMark2.csv")

blup2020 <- read.csv("../npq/data/work/blups/blup2020gwasRNA.csv")
blup2021 <- read.csv("../npq/data/work/blups/blup2021gwasRNA.csv")

RMIP <- read.csv("../npq/results/data/RMIP2020.csv")

blup2020$Year <- "2020"
blup2021$Year <- "2021"

res <- rbind(blup2020, blup2021)

npq <- read.csv("../npq/data/work/NpqRaw2020.csv")
colnames(npq)[9:27] <- c(0.266333, 0.41333, 0.746667, 1.08, 2.08, 3.08, 4.08, 5.08 , 6.08, 7.08 , 8.08, 
                         9.08, 10.08, 10.26333, 10.41333, 10.91333, 11.91333, 15.08, 20.08)

npq <- npq %>% 
  select(Taxa,9:27) %>%
  pivot_longer(2:20) %>%
  group_by(Taxa, name) %>%
  summarise(NPQ=mean(value))

colnames(npq)[2] <- "Time"
npq$Time <- as.numeric(npq$Time)

qy2020 <- read_csv("../npq/data/work/Qy_Raw2020.csv")
colnames(qy2020)[9:27] <- c(0.266333, 0.41333, 0.746667, 1.08, 2.08, 3.08, 4.08, 5.08 , 6.08, 7.08 , 8.08, 9.08, 
                            10.08, 10.26333, 10.41333, 10.91333, 11.91333, 15.08, 20.08)

qy2020 <- qy2020 %>% 
  select(Taxa,9:27) %>%
  pivot_longer(2:20) %>%
  group_by(Taxa, name) %>%
  summarise(QY=mean(value))

colnames(qy2020)[2] <- "Time"
qy2020$Time <- as.numeric(qy2020$Time)

res$rs3_176152936 <- rep(geno[which(map$SNP=="rs3_176152936"),], 2)
res$rs3_176152936[res$rs3_176152936==0] <- "C"
res$rs3_176152936[res$rs3_176152936==2] <- "T"

res$rs2_182382835 <- rep(geno[which(map$SNP=="rs2_182382835"),], 2)
res$rs2_182382835[res$rs2_182382835==0] <- "C"
res$rs2_182382835[res$rs2_182382835==2] <- "G"

res$rs3_147531265 <- rep(geno[which(map$SNP=="rs3_147531265"),], 2)
res$rs3_147531265[res$rs3_147531265==0] <- "T"
res$rs3_147531265[res$rs3_147531265==2] <- "G"

res$rs7_179281748 <- rep(geno[which(map$SNP=="rs7_179281748"),], 2)
res$rs7_179281748[res$rs7_179281748==0] <- "C"
res$rs7_179281748[res$rs7_179281748==2] <- "T"

res$rs5_187859720 <- rep(geno[which(map$SNP=="rs5_187859720"),], 2)
res$rs5_187859720[res$rs5_187859720==0] <- "G"
res$rs5_187859720[res$rs5_187859720==2] <- "A"

res$rs1_86064177 <- rep(geno[which(map$SNP=="rs1_86064177"),], 2)
res$rs1_86064177[res$rs1_86064177==0] <- "G"
res$rs1_86064177[res$rs1_86064177==2] <- "A"

res$rs3_177340601 <- rep(geno[which(map$SNP=="rs3_177340601"),], 2)
res$rs3_177340601[res$rs3_177340601==0] <- "A"
res$rs3_177340601[res$rs3_177340601==2] <- "G"


x <- merge(npq, res[,c(1,35:42)], by="Taxa")
x2 <- merge(qy2020, res[,c(1,35:42)], by="Taxa")

###################
### PsbS allele ###
###################

gNPQ0 <- ggplot() + 
  stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_177340601)), fun=mean, geom="point", size=2.5, alpha=0.6) + 
  stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_177340601)), fun.data="mean_se", geom="linerange", size=1.1, width=0.25) + 
  #stat_summary(data=x, aes(Time, NPQ, fill=paste(rs3_176152936)), fun.data="mean_sdl", geom="ribbon", size=1.1, alpha=0.6) + 
  #stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_176152936)), fun=mean, geom="line", size=1.1, alpha=0.6) + 
  geom_segment(aes(y=0, yend=0, x=10.25, xend=20.5), size=5) + # Mark dark part 
  geom_segment(aes(y=0, yend=0, x=0, xend=10.25), size=5, colour="gold") +
  xlab("Time (minutes)") + 
  ylab("NPQ") + 
  scale_color_d3(name="3:177340601") + 
  scale_fill_d3() + 
  theme(legend.position = c(0.9,0.9), legend.key.size = unit(1, "cm"))
# guides(shape = guide_legend(override.aes = list(size = 10)))

gO1 <- ggplot(res, aes(Year, NPQmaxme, fill=rs3_177340601)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ[max]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(NPQmaxme, na.rm=T)+0.1, yend=max(NPQmaxme, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(NPQmaxme, na.rm=T)+0.1, yend=max(NPQmaxme, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$NPQmaxme, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$NPQmaxme, na.rm=T)+0.1, label="*", size=10) + 
  ylim(min(res$NPQmaxme, na.rm=T)-0.2,max(res$NPQmaxme, na.rm=T)+0.4)

gO2 <- ggplot(res, aes(Year, cf3npq, fill=rs3_177340601)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~start[darkH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(cf3npq, na.rm=T)+0.1, yend=max(cf3npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(cf3npq, na.rm=T)+0.1, yend=max(cf3npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$cf3npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$cf3npq, na.rm=T)+0.1, label="*", size=10) +
  ylim(min(res$cf3npq, na.rm=T)-0.2,max(res$cf3npq, na.rm=T)+0.4)


gO3 <- ggplot(res, aes(Year, bf3npq, fill=rs3_177340601)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~amplitude[darkH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(bf3npq, na.rm=T)+0.1, yend=max(bf3npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(bf3npq, na.rm=T)+0.1, yend=max(bf3npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$bf3npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$bf3npq, na.rm=T)+0.1, label="**", size=10) + 
  ylim(min(res$bf3npq, na.rm=T)-0.2,max(res$bf3npq, na.rm=T)+0.4)


gO4 <- ggplot(res, aes(Year, af4npq, fill=rs3_177340601)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~amplitude[darkE]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(af4npq, na.rm=T)+0.1, yend=max(af4npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(af4npq, na.rm=T)+0.1, yend=max(af4npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$af4npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$af4npq, na.rm=T)+0.1, label="**", size=10)+ 
  ylim(min(res$af4npq, na.rm=T)-0.2,max(res$af4npq, na.rm=T)+0.4)


gNPQ0 | ((gO1 | gO2) / (gO3 | gO4))  + plot_annotation(tag_levels = "a")

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/Fig2/PsbS_effect.svg", units = "mm", width = 75, height = 40, scale = 4)

###################
### PSI3 allele ###
###################

gNPQ1 <- ggplot() + 
  stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_176152936)), fun=mean, geom="point", size=2.5, alpha=0.6) + 
  stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_176152936)), fun.data="mean_se", geom="linerange", size=1.1, width=0.25) + 
  #stat_summary(data=x, aes(Time, NPQ, fill=paste(rs3_176152936)), fun.data="mean_sdl", geom="ribbon", size=1.1, alpha=0.6) + 
  #stat_summary(data=x, aes(Time, NPQ, colour=paste(rs3_176152936)), fun=mean, geom="line", size=1.1, alpha=0.6) + 
  geom_segment(aes(y=0, yend=0, x=10.25, xend=20.5), size=5) + # Mark dark part 
  geom_segment(aes(y=0, yend=0, x=0, xend=10.25), size=5, colour="gold") +
  xlab("Time (minutes)") + 
  ylab("NPQ") + 
  scale_color_d3(name="3:176,152,936") + 
  scale_fill_d3() + 
  theme(legend.position = c(0.9,0.9), legend.key.size = unit(1, "cm"))

gA1 <- ggplot(res, aes(Year, bf1npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~asymptote[lightH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(bf1npq, na.rm=T)+0.2, yend=max(bf1npq, na.rm=T)+0.2)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(bf1npq, na.rm=T)+0.2, yend=max(bf1npq, na.rm=T)+0.2)) + 
  annotate("text", x=1, y=max(res$bf1npq, na.rm=T)+0.2, label="***", size=10) + 
  annotate("text", x=2, y=max(res$bf1npq, na.rm=T)+0.2, label="***", size=10) + 
  xlab("")

gA2 <- ggplot(res, aes(Year, af2npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~asymptote[lightE]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(af2npq, na.rm=T)+0.2, yend=max(af2npq, na.rm=T)+0.2)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(af2npq, na.rm=T)+0.2, yend=max(af2npq, na.rm=T)+0.2)) + 
  annotate("text", x=1, y=max(res$af2npq, na.rm=T)+0.2, label="***", size=10) + 
  annotate("text", x=2, y=max(res$af2npq, na.rm=T)+0.2, label="***", size=10) + 
  xlab("")

gA3 <- ggplot(res, aes(Year, NPQmaxme, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ[max]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(NPQmaxme, na.rm=T)+0.1, yend=max(NPQmaxme, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(NPQmaxme, na.rm=T)+0.1, yend=max(NPQmaxme, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$NPQmaxme, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$NPQmaxme, na.rm=T)+0.1, label="***", size=10) + 
  xlab("") +
  ylim(min(res$NPQmaxme, na.rm=T)-0.2,max(res$NPQmaxme, na.rm=T)+0.4)


gA4 <- ggplot(res, aes(Year, cf3npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~start[darkH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(cf3npq, na.rm=T)+0.1, yend=max(cf3npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(cf3npq, na.rm=T)+0.1, yend=max(cf3npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$cf3npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$cf3npq, na.rm=T)+0.1, label="***", size=10) + 
  xlab("") + 
  ylim(min(res$cf3npq, na.rm=T)-0.2,max(res$cf3npq, na.rm=T)+0.4)

gA5 <- ggplot(res, aes(Year, af3npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~slope[darkH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(af3npq, na.rm=T)+0.1, yend=max(af3npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(af3npq, na.rm=T)+0.1, yend=max(af3npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$af3npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$af3npq, na.rm=T)+0.1, label="***", size=10) + 
  ylim(min(res$af3npq, na.rm=T)-0.2,max(res$af3npq, na.rm=T)+0.4)


gA6 <- ggplot(res, aes(Year, npqslop2me, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~slope[darkL]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(npqslop2me, na.rm=T)+0.2, yend=max(npqslop2me, na.rm=T)+0.2)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(npqslop2me, na.rm=T)+0.2, yend=max(npqslop2me, na.rm=T)+0.2)) + 
  annotate("text", x=1, y=max(res$npqslop2me, na.rm=T)+0.2, label="***", size=10) + 
  annotate("text", x=2, y=max(res$npqslop2me, na.rm=T)+0.2, label="***", size=10)

gA7 <- ggplot(res, aes(Year, bf3npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~amplitude[darkH]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(bf3npq, na.rm=T)+0.1, yend=max(bf3npq, na.rm=T)+0.1)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(bf3npq, na.rm=T)+0.1, yend=max(bf3npq, na.rm=T)+0.1)) + 
  annotate("text", x=1, y=max(res$bf3npq, na.rm=T)+0.1, label="***", size=10) + 
  annotate("text", x=2, y=max(res$bf3npq, na.rm=T)+0.1, label="***", size=10) + 
  ylim(min(res$bf3npq, na.rm=T)-0.2,max(res$bf3npq, na.rm=T)+0.4)


gA8 <- ggplot(res, aes(Year, af4npq, fill=rs3_176152936)) + 
  geom_boxplot() +
  theme(legend.position = 'none') + 
  scale_fill_d3() + 
  ylab(expression(italic(NPQ~amplitude[darkE]))) + 
  geom_segment(aes(x=0.75, xend=1.25, y=max(af4npq, na.rm=T)+0.2, yend=max(af4npq, na.rm=T)+0.2)) + 
  geom_segment(aes(x=1.75, xend=2.25, y=max(af4npq, na.rm=T)+0.2, yend=max(af4npq, na.rm=T)+0.2)) + 
  annotate("text", x=1, y=max(res$af4npq, na.rm=T)+0.2, label="***", size=10) + 
  annotate("text", x=2, y=max(res$af4npq, na.rm=T)+0.2, label="***", size=10)

(gNPQ1 + plot_spacer()) / (gA1 | gA2 | gA3 | gA4) / (gA5 | gA6 | gA7 | gA8) + plot_annotation(tag_levels = "a")

gNPQ1 | ((gA3 | gA4) / (gA5 | gA7))  + plot_annotation(tag_levels = "a")


ggsave("../Data/Maize_NPQ_Natural_Variation/figures/Fig2/PSI3_effect.svg", units = "mm", width = 75, height = 40, scale = 4)

