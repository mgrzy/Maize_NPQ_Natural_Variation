library(data.table)
library(tidyverse)
library(patchwork)

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
  geom_point(size=4, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  xlim(175652936, 178038945) + ylab("RMIP\n2020") + 
  scale_color_manual(values = c("yellow", "orange", "darkgreen", "green",
                                        "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) # colour values for each trait

g.psi.rmip2021 <- ggplot(psi.psbs$rmip2021, aes(POS, RMIP)) + 
  geom_point(size=4, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = c(0.9,0.7)) +
  xlim(175652936, 178038945) + ylab("RMIP\n2021") + 
  scale_color_manual(values = c("yellow", "orange", "darkgreen", "green",
                                        "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) # colour values for each trait

g.psi.LD <- ggplot(psi.psbs$LD, aes(BP_B, R2)) + 
  geom_point(size=4, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  xlab("Chromosome 3") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  xlim(c(175652936, 178038945))

g.psipsbs.gene <- ggplot() + 
  geom_point(data=psi.psbs$gen[psi.psbs$gen$V7=="+",], aes(x=V4, y=0.5), shape="\u25BA", size=6) +
  geom_point(data=psi.psbs$gen[psi.psbs$gen$V7=="-",], aes(x=V4, y=-0.5), shape="\u25C4", size=6) +
  geom_point(data=psi.psbs$psbs, aes(x=V4, y=-0.5), shape="\u25C4", size=6, colour="red") +
  geom_point(data=psi.psbs$psi, aes(x=V4, y=0.5), shape="\u25BA", size=6, colour="red") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-2, 2) + 
  scale_x_continuous(labels = paste0(c("176,000", "176,500", "177,000", "177,500"), " kB"),
                     breaks = c(176000000, 176500000, 177000000,177500000), limits = c(175652936, 178038945)) + 
  xlab("Chromosome 3")

g.psbs <- ggplot() + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.5, "cm"), ends = "first")) + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  geom_segment(aes(y=0.5, yend=-0.18, x=179566467, xend=179566467), size=1, colour="red") + 
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
  geom_segment(aes(y=0.5, yend=-0.18, x=178405582, xend=178405582), size=0.5, colour="red") + 
  geom_segment(aes(y=0.5, yend=-0.18, x=178405643, xend=178405643), size=0.5, colour="red") + 
  geom_segment(aes(y=0.5, yend=-0.18, x=178407540, xend=178407540), size=0.5, colour="red") + 
  geom_segment(aes(y=0.5, yend=-0.18, x=178407565, xend=178407565), size=0.5, colour="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)


g.12 <- (g.psi.rmip2020 / g.psi.rmip2021 / g.psi.LD / g.psipsbs.gene) / (g.psi + g.psbs) + 
  plot_layout(heights = c(2,2,2,2,1.5))

ggsave(plot = g.12, filename = "../Data/Maize_NPQ_Natural_Variation/figures/PSI3_PsbS.svg", device = "svg", 
       units = "mm", width = 75, height = 40, scale = 4)

### IRM1
irm1 <- list()
irm1$rmip2020 <- RMIP_2020[CHROM==2  & POS > 182332835 & POS < 182432835,]
irm1$gen <- gff[V1=="Chr2" & V4 > 182332835 & V5 < 182432835]
irm1$irm1_v5 <- gff2[grep("Zm00001eb098640", gff2$V9)]

irm1$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs2_182382835.ld")

g.irm1.rmip2020 <- ggplot(irm1$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, colour="yellow") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  xlim(182382305, 182384068) + 
  ylab("RMIP\n2020")

g.irm1.LD <- ggplot(irm1$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
  #      axis.ticks.x = element_blank()) + 
  scale_x_continuous(labels = paste0(c("182,382,800", "182,383,200", "182,383,600")),
                   breaks = c(182382800, 182383200, 182383600), limits = c(182382305, 182384068)) + 
  xlab("Chromosome 2")

g.irm1.gene <- ggplot() + 
  geom_segment(data=irm1$gen[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=irm1$gen[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
         axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

g.irm1 <- ggplot() + 
  geom_segment(data=irm1$irm1_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.5, "cm"), ends = "last")) + 
  geom_segment(data=irm1$irm1_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=181981838, y=-0.75, label=expression(Zm00001d005657~"("*italic(IRM1)*")")) + 
  geom_segment(aes(y=0.5, yend=0, x=181981544, xend=181981544), size=2, colour="red") + 
  geom_segment(aes(y=0.5, yend=0, x=181981600, xend=181981600), size=2, colour="red") + 
  geom_segment(aes(y=0.5, yend=0, x=181981652, xend=181981652), size=2, colour="red") + 
  geom_segment(aes(y=0.5, yend=0, x=181982035, xend=181982035), size=2, colour="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.3 <- g.irm1.rmip2020 / g.irm1.LD / g.irm1 + plot_layout(heights = c(4,4,1))

ggsave(plot = g.3, filename = "../Data/Maize_NPQ_Natural_Variation/figures/IRM1.svg", device = "svg", width = 12, height = 6)

### TRX
trx <- list()
trx$rmip2020 <- RMIP_2020[CHROM==3  & POS > 147481265 & POS < 147581265,]
trx$gen <- gff[V1=="Chr3" & V4 > 147481265 & V5 < 147581265]
trx$trx_v5 <- gff2[grep("Zm00001eb140270", gff2$V9)]
trx$trx <- gff[grep("Zm00001d042017", gff$V9)]

trx$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_147531265.ld")

g.trx.rmip2020 <- ggplot(trx$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.position = c(0.9, 0.8)) + 
  xlim(147481265, 147581265) + 
  scale_color_manual(values = c("darkgreen", "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) + 
  ylab("RMIP\n2020")


g.trx.LD <- ggplot(trx$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  xlab("Chromosome 2") + 
  xlim(147481265, 147581265) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())

trx$gen[grep("Zm00001d042016", trx$gen$V9),]$V8 <- -1
trx$gen[grep("Zm00001d042017", trx$gen$V9),]$V8 <- -1
trx$gen[grep("Zm00001d042018", trx$gen$V9),]$V8 <- 1
trx$gen[grep("Zm00001d042019", trx$gen$V9),]$V8 <- -1

trx$gen$V8 <- as.numeric(trx$gen$V8)
trx$gen <- trx$gen[!which(is.na(trx$gen$V8)),]

g.trx.gene <- ggplot() + 
  geom_segment(data=trx$gen[V3=="gene",], aes(x=V4, xend=V5, y=V8, yend=V8)) + 
  geom_segment(data=trx$gen[V3=="CDS",], aes(x=V4, xend=V5, y=V8, yend=V8), size=5) + 
  geom_segment(data=trx$trx[V3=="gene",], aes(x=V4, xend=V5, y=-1, yend=-1), colour="red") + 
  geom_segment(data=trx$trx[V3=="CDS",], aes(x=V4, xend=V5, y=-1, yend=-1), size=5, colour="red") + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("147,500", "147,525", "147,550", "147,575"), " kB"), 
                     breaks = c(147500000, 147525000, 147550000, 147575000), limits = c(147481265, 147581265)) +
  ylim(-2,2) +
  xlab("Chromosome 3")

g.trx <- ggplot() + 
  geom_segment(data=trx$trx_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.5, "cm"), ends = "first")) + 
  geom_segment(data=trx$trx_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=148314515, y=-0.75, label=expression(Zm00001d042017~"("*italic(TRX-Y1)*")")) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1) 

g.4 <- g.trx.rmip2020 / g.trx.LD / g.trx.gene / g.trx + plot_layout(heights = c(4,4,1,1))

ggsave(plot = g.4, filename = "../Data/Maize_NPQ_Natural_Variation/figures/TRX.svg", device = "svg", width = 12, height = 6)

### ACHT3
acht3 <- list()
acht3$rmip2020 <- RMIP_2020[CHROM==7  & POS > 179231748 & POS < 179331748,]
acht3$gen <- gff[V1=="Chr7" & V4 > 179231748 & V5 < 179331748]
acht3$acht3_v5 <- gff2[grep("Zm00001eb330920", gff2$V9)]
acht3$acht3 <- gff[grep("Zm00001d022518", gff$V9)]

acht3$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs7_179281748.ld")

g.acht3.rmip2020 <- ggplot(acht3$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.position = c(0.9, 0.8)) + 
  xlim(179271748, 179291748) + 
  scale_color_manual(values = c("darkgreen", "green"),
                     name="Trait related to", na.translate = F, label=greeks) + 
  ylab("RMIP\n2020")


g.acht3.LD <- ggplot(acht3$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  xlim(c(179271748, 179291748)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())
 
g.acht3.gene <- ggplot() + 
  geom_segment(data=acht3$gen[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1)) + 
  geom_segment(data=acht3$gen[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5) + 
  geom_segment(data=acht3$acht3[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1), colour="red") + 
  geom_segment(data=acht3$acht3[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5, colour="red") + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("179,275", "179,280", "179,285", "179,290"), " kB"),
                     breaks = c(179275000, 179280000, 179285000, 179900000), limits = c(179271748, 179291748)) + 
  xlab("Chromosome 7")
  
g.acht3.acht <- ggplot() + 
  geom_segment(data=acht3$acht3_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=acht3$acht3_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=182742282, y=-0.8, label=expression(Zm00001d022518~"("*italic(ACHT3)*")")) + 
  geom_segment(aes(y=1, yend=0, x=182742647, xend=182742647), size=2, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=182742793, xend=182742793), size=2, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=182742814, xend=182742814), size=2, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=182742836, xend=182742836), size=2, colour="red") + 
  geom_point(aes(y=0.5, x=182742510), shape=25, fill="red", size=4) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.5 <- g.acht3.rmip2020 / g.acht3.LD / g.acht3.gene / g.acht3.acht + plot_layout(heights = c(4,4, 1, 1))

ggsave(plot = g.5, filename = "../Data/Maize_NPQ_Natural_Variation/figures/ACHT3.svg", device = "svg", width = 12, height = 6)

### PMI1
PMI1 <- list()
PMI1$rmip2020 <- RMIP_2020[CHROM==1  & POS > 86054177 & POS < 86074177,]
PMI1$gen <- gff[V1=="Chr1" & V4 > 86054177 & V5 < 86074177]
PMI1$PMI1_v5 <- gff2[grep("Zm00001eb022210", gff2$V9)]
PMI1$PMI1 <- gff[grep("Zm00001d029761", gff$V9)]

PMI1$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs1_86064177.ld")

g.PMI1.rmip2020 <- ggplot(PMI1$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.position = c(0.9, 0.8)) + 
  xlim(86054177, 86074177) + 
  scale_color_manual(values = c("grey50"),
                     name="Trait related to", na.translate = F, label=greeks) + 
  ylab("RMIP\n2020")

g.PMI1.LD <- ggplot(PMI1$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  scale_x_continuous(labels = paste0(c("860,64", "860,65", "860,62"), " kB"),
                                        breaks = c(86064000, 86065000, 86062000), limits = c(86062903, 86065425)) + 
  xlab("Chromosome 1")
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
  #      axis.ticks.x = element_blank())

g.PMI1.gene <- ggplot() + 
  geom_segment(data=PMI1$gen[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1)) + 
  geom_segment(data=PMI1$gen[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5) + 
  geom_segment(data=PMI1$PMI1[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1), colour="red") + 
  geom_segment(data=PMI1$PMI1[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5, colour="red") + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  #scale_x_continuous(labels = paste0(c("179,275", "179,280", "179,285", "179,290"), " kB"),
  #                   breaks = c(179275000, 179280000, 179285000, 179900000), limits = c(179271748, 179291748)) + 
  xlab("Chromosome 1")

mutPSI <- read.table("psi1.txt", header = F)

g.PMI1.pmi1 <- ggplot() + 
  geom_segment(data=PMI1$PMI1_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=PMI1$PMI1_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=85241864, y=-0.8, label=expression(Zm00001d029761~"("*italic(PMI1)*")")) + 
  geom_point(aes(y=0.5, x=85243339), shape=25, fill="red", size=4) + 
  geom_segment(data=mutPSI, aes(y=1, yend=0, x=V1, xend=V1), size=2, colour="red") + 
  #geom_segment(aes(y=1, yend=0, x=182742793, xend=182742793), size=2, colour="red") + 
  #geom_segment(aes(y=1, yend=0, x=182742814, xend=182742814), size=2, colour="red") + 
  #geom_segment(aes(y=1, yend=0, x=182742836, xend=182742836), size=2, colour="red") + 
  #geom_point(aes(y=0.5, x=182742510), shape=25, fill="red", size=4) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.6 <- g.PMI1.rmip2020 / g.PMI1.LD / g.PMI1.pmi1 + plot_layout(heights = c(4,4, 1, 1))

ggsave(plot = g.6, filename = "../Data/Maize_NPQ_Natural_Variation/figures/PMI1.svg", device = "svg", width = 12, height = 6)

### OEP37
OEP37 <- list()
OEP37$rmip2020 <- RMIP_2020[CHROM==5  & POS > 187809720 & POS < 187909720,]
OEP37$gen <- gff[V1=="Chr5" & V4 > 187809720 & V5 < 187909720]
OEP37$OEP37_v5 <- gff2[grep("Zm00001eb246750", gff2$V9)]
OEP37$OEP37 <- gff[grep("Zm00001d017171", gff$V9)]

OEP37$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs5_187859720.ld")

g.OEP37.rmip2020 <- ggplot(OEP37$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), legend.position = c(0.9, 0.8)) + 
  xlim(187809720, 187909720) + 
  scale_color_manual(values = c("grey50"),
                     name="Trait related to", na.translate = F, label=greeks) + 
  ylab("RMIP\n2020")

g.OEP37.LD <- ggplot(OEP37$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  xlim(187809720, 187909720) + 
  #scale_x_continuous(labels = paste0(c("860,64", "860,65", "860,62"), " kB"),
  #                   breaks = c(86064000, 86065000, 86062000), limits = c(86062903, 86065425)) + 
  xlab("Chromosome 5") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
      axis.ticks.x = element_blank())

g.OEP37.gene <- ggplot() + 
  geom_segment(data=OEP37$gen[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1)) + 
  geom_segment(data=OEP37$gen[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5) + 
  geom_segment(data=OEP37$OEP37[V3=="gene",], aes(x=V4, xend=V5, y=1, yend=1), colour="red") + 
  geom_segment(data=OEP37$OEP37[V3=="CDS",], aes(x=V4, xend=V5, y=1, yend=1), size=5, colour="red") + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("187,830", "187,860", "187,890"), " kB"),
                     breaks = c(187830000, 187860000, 187890000), limits = c(187809720, 187909720)) + 
  xlab("Chromosome 5") 

g.OEP37.OEP37 <- ggplot() + 
  geom_segment(data=OEP37$OEP37_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=OEP37$OEP37_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("text", x=186954049, y=-0.8, label=expression(Zm00001d017171~"("*italic(OEP37)*")")) + 
  #geom_point(aes(y=0.5, x=85243339), shape=25, fill="red", size=4) + 
  geom_segment(aes(y=1, yend=0, x=186953749, xend=186953749), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186952808, xend=186952808), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186952822, xend=186952822), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186952823, xend=186952823), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186955468, xend=186955468), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186955834, xend=186955834), size=1, colour="red") + 
  geom_segment(aes(y=1, yend=0, x=186955859, xend=186955859), size=1, colour="red") + 
  #geom_point(aes(y=0.5, x=182742510), shape=25, fill="red", size=4) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.7 <- g.OEP37.rmip2020 / g.OEP37.LD / g.OEP37.gene /g.OEP37.OEP37 + plot_layout(heights = c(4,4, 1, 1))

ggsave(plot = g.7, filename = "../Data/Maize_NPQ_Natural_Variation/figures/OEP37.svg", device = "svg", width = 12, height = 6)
