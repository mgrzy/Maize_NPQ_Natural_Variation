library(data.table)
library(tidyverse)
library(patchwork)
library(ggsci)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

a <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv503.csv.gz", data.table = F)

a <- a[a$chromosome %in% c(1:10),]
a$chromosome <- as.numeric(a$chromosome)
colnames(a)[2] <- "CHROM"
colnames(a)[4] <- "POS"


nCHR <- length(unique(a$CHROM))
a$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(a$CHROM))){
  nbp[i] <- max(a[a$CHROM == i,]$POS)
  a[a$CHROM == i,"BPcum"] <- a[a$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- a %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

g1 <- ggplot(a, aes(x = BPcum, y = -log10(npqslop1me))) +
  geom_hline(yintercept = 5.73, linetype=2) + 
  geom_segment(aes(x = 728648553, xend=728648553, y=2, yend=7.5), linewidth=0.5) + 
  annotate("text", x=728648553, y=7.5, hjust=1, label=expression(italic(PSBS)), size=5) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g2 <- ggplot(a, aes(x = BPcum, y = -log10(NPQmaxme))) +
  geom_hline(yintercept = 5.73, linetype=2) + 
  geom_segment(aes(x = 728648553, xend=728648553, y=2, yend=7.5), linewidth=0.5) + 
  annotate("text", x=728648553, y=7.5, hjust=1, label=expression(italic(PSBS)), size=5) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g3 <- ggplot(a, aes(x = BPcum, y = -log10(npqslop2me))) +
  geom_hline(yintercept = 5.73, linetype=2) + 
  geom_segment(aes(x = 728648553, xend=728648553, y=2, yend=7.5), linewidth=0.5) + 
  annotate("text", x=728648553, y=7.5, hjust=1, label=expression(italic(PSBS)), size=5) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g4 <- ggplot(a, aes(x = BPcum, y = -log10(NPQendme))) +
  geom_hline(yintercept = 5.73, linetype=2) + 
  geom_segment(aes(x = 728648553, xend=728648553, y=2, yend=7.5), linewidth=0.5) + 
  annotate("text", x=728648553, y=7.5, hjust=1, label=expression(italic(PSBS)), size=5) + 
  geom_point(aes(colour=as.character(CHROM)), size=3) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none")

g1 + g2 + g3 + g4 + plot_annotation(tag_levels = "a")

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/TWAS.png", units = "mm", width = 80, height = 40, scale = 4)
