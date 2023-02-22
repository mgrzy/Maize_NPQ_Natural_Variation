library(data.table)
library(tidyverse)
library(patchwork)
library(jcolors)

theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

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

axis.set <- b %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

g.psbs <- b %>%
  filter(-log10(Zm00001d042697.MLM)>2) %>%
  ggplot(aes(x = BPcum, y = -log10(Zm00001d042697.MLM))) +
  geom_point(aes(colour=as.character(CHROM))) + # points for RMIP values
  #geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  #scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = expression(-log[10](italic(p)))) + 
  theme(legend.position = "none") 

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
