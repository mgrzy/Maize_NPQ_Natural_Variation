library(data.table)
library(tidyverse)
library(patchwork)

# Set them to classic and axis text to black rather then grey
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

# Load data

gff <- fread("../../BigData/Maize/v4/GFF/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3.gz",skip = 5, fill = T)
gff2 <- fread("../../BigData/Maize/v5/GFF/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz",skip = 5, fill = T)

###

RMIP_2020 <- fread("../Data/Maize_NPQ_Natural_Variation/data/figures/Fig1c_RMIP.csv")
RMIP_2021 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/RMIPclean2021.csv")

###

RMIP_2020 <- RMIP_2020[RMIP>1,]

### PSI-PsbS
psi.psbs <- list()
psi.psbs$rmip2020 <- RMIP_2020[CHROM==3  & POS > 175652936 & POS < 178038945,]
psi.psbs$rmip2021 <- RMIP_2021[CHROM==3  & POS > 175652936 & POS < 178038945,]
psi.psbs$gen <- gff[V1=="Chr3" & V4 > 175652936 & V5 < 178038945]
psi.psbs$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_176152936.ld")
psi.psbs$psbs <- psi.psbs$gen[grep("Zm00001d042697",psi.psbs$gen$V9)]
psi.psbs$psi <- psi.psbs$gen[grep("Zm00001d042669",psi.psbs$gen$V9)]
psi.psbs$psbs_v5 <- gff2[grep("Zm00001eb146510", gff2$V9)]

g.psi.rmip2020 <- ggplot(psi.psbs$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = c(0.8,0.8)) +
  xlim(175652936, 178038945)

g.psi.rmip2021 <- ggplot(psi.psbs$rmip2021, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = c(0.8,0.8)) +
  xlim(175652936, 178038945)

g.psi.LD <- ggplot(psi.psbs$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
  #      axis.ticks.x = element_blank()) + 
  #scale_x_continuous(labels = paste0(c("182,382,500", "182,383,000", "182,383,500"), " "),
  #                   breaks = c(182382500, 182383000, 182383500)) + 
  xlab("Chromosome 3") + 
  xlim(175652936, 178038945)

g.psipsbs.gene <- ggplot() + 
  geom_segment(data=psi.psbs$gen[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_segment(data=psi.psbs$psbs[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.25, "cm")), colour="red") + 
  geom_segment(data=psi.psbs$psi[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.25, "cm")), colour="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  xlim(175652936, 178038945)

g.psbs <- ggplot() + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=psi.psbs$psbs_v5[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  annotate("point", y=0.75, x=179566467, shape=25, size=5, fill="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  ylim(-1, 1)

g.psi <- ggplot() + 
  geom_segment(data=psi.psbs$psi[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=psi.psbs$psi[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

(g.psi.rmip2020 / g.psi.rmip2021 / g.psi.LD / g.psipsbs.gene) / (g.psi + g.psbs) + 
  plot_layout(heights = c(4,4,3,0.5,1))

###
irm1 <- list()
irm1$rmip2020 <- RMIP_2020[CHROM==2  & POS > 182332835 & POS < 182432835,]
irm1$gen <- gff[V1=="Chr2" & V4 > 182332835 & V5 < 182432835]
irm1$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs2_182382835.ld")

g.irm1.rmip2020 <- ggplot(irm1$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, colour="yellow") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  xlim(182382424, 182383814)

g.irm1.LD <- ggplot(irm1$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
  #      axis.ticks.x = element_blank()) + 
scale_x_continuous(labels = paste0(c("182,382,500", "182,383,000", "182,383,500"), " "),
                   breaks = c(182382500, 182383000, 182383500)) + 
  xlab("Chromosome 2")

g.irm1.gene <- ggplot() + 
  geom_segment(data=irm1$gen[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=irm1$gen[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
         axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

g.irm1.rmip2020 / g.irm1.LD / g.irm1.gene + plot_layout(heights = c(4,4,1))

###

acht3 <- list()
acht3$rmip2020 <- RMIP_2020[CHROM==7  & POS > 179231748 & POS < 179331748,]
acht3$gen <- gff[V1=="Chr7" & V4 > 179231748 & V5 < 179331748]
acht3$LD <- fread("../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs7_179281748.ld")

g.acht3.rmip2020 <- ggplot(acht3$rmip2020, aes(POS, RMIP)) + 
  geom_point(size=5, alpha=0.8, colour="yellow") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  xlim(179271748, 179291748)

g.acht3.LD <- ggplot(acht3$LD, aes(BP_B, R2)) + 
  geom_point(size=5, alpha=0.8, colour="dodgerblue4") + 
  ylab(expression(LD~(r^2))) + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
  #      axis.ticks.x = element_blank()) + 
  #scale_x_continuous(labels = paste0(c("182,382,500", "182,383,000", "182,383,500"), " "),
  #                   breaks = c(182382500, 182383000, 182383500)) + 
  xlab("Chromosome 2") + 
  xlim(179271748, 179291748)

g.acht3.gene <- ggplot() + 
  geom_segment(data=acht3$gen[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(0.25, "cm"))) + 
  geom_segment(data=acht3$acht3[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
             arrow = arrow(length = unit(0.25, "cm")), colour="red") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank()) + 
  xlim(179271748, 179291748)

g.acht3.acht <- ggplot() + 
  geom_segment(data=acht3$acht3[V3=="gene",], aes(x=V4, xend=V5, y=0, yend=0), 
               arrow = arrow(length = unit(1, "cm"))) + 
  geom_segment(data=acht3$acht3[V3=="CDS",], aes(x=V4, xend=V5, y=0, yend=0), size=5) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

g.acht3.rmip2020 / g.acht3.LD / g.acht3.gene / g.acht3.acht + plot_layout(heights = c(4,4, 1, 1))
