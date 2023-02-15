library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)

# Set them to classic and axis text to black rather then grey
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

res <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv503.csv.gz", data.table = F)
res2 <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv942.csv.gz", data.table = F)
res3 <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2021/TWAS_WiDiv503_2021.csv.gz", data.table = F)
res4 <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2021/TWAS_WiDiv942_2021.csv.gz", data.table = F)

res$chromosome <- as.numeric(res$chromosome)
res <- res[res$chromosome %in% c(1:10),]

res2$chromosome <- as.numeric(res2$chromosome)
res2 <- res2[res2$chromosome %in% c(1:10),]

res3$chromosome <- as.numeric(res3$chromosome)
res3 <- res3[res3$chromosome %in% c(1:10),]

res4$chromosome <- as.numeric(res4$chromosome)
res4 <- res4[res4$chromosome %in% c(1:10),]

nCHR <- length(unique(res$chromosome))
res$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(res$chromosome))){
  nbp[i] <- max(res[res$chromosome == i,]$position_left)
  res[res$chromosome == i,"BPcum"] <- res[res$chromosome == i,"position_left"] + s
  s <- s + nbp[i]
}

nCHR <- length(unique(res2$chromosome))
res2$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(res2$chromosome))){
  nbp[i] <- max(res2[res2$chromosome == i,]$position_left)
  res2[res2$chromosome == i,"BPcum"] <- res2[res2$chromosome == i,"position_left"] + s
  s <- s + nbp[i]
}

nCHR <- length(unique(res$chromosome))
res3$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(res3$chromosome))){
  nbp[i] <- max(res3[res3$chromosome == i,]$position_left)
  res3[res3$chromosome == i,"BPcum"] <- res3[res3$chromosome == i,"position_left"] + s
  s <- s + nbp[i]
}

nCHR <- length(unique(res$chromosome))
res4$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(res4$chromosome))){
  nbp[i] <- max(res4[res4$chromosome == i,]$position_left)
  res4[res4$chromosome == i,"BPcum"] <- res4[res4$chromosome == i,"position_left"] + s
  s <- s + nbp[i]
}

axis.set <- res %>% 
  group_by(chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

trait <- colnames(res)[6:36]

for (i in 6:36){
  
g1 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res)) , linetype=2) +
  geom_point(data=res, aes(BPcum, -log10(res[,i]), colour=as.character(chromosome))) + 
  scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome")

g2 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res2)) , linetype=2) +
  geom_point(data=res2, aes(BPcum, -log10(res2[,i]), colour=as.character(chromosome))) + 
  scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome")

g3 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res3)) , linetype=2) +
  geom_point(data=res3, aes(BPcum, -log10(res3[,i]), colour=as.character(chromosome))) + 
  scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome")

g4 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res4)) , linetype=2) +
  geom_point(data=res4, aes(BPcum, -log10(res4[,i]), colour=as.character(chromosome))) + 
  scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome")

gg <- g1 + g2 + g3 + g4 + plot_annotation(tag_levels = "a")

ggsave(plot = gg, paste0("../Data/Maize_NPQ_Natural_Variation/figures/TWAS/", trait[i-5], ".png"), width = 10)

}

blup <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2020gwasRNA.csv")
blup2 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2021gwasRNA.csv")
ex <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Expression_MetaData_csv.gz")

a <- ex[,c("Taxa", "Study", "Zm00001d042697")]

a <- merge(a, blup[,c("Taxa", "NPQmaxme")], by="Taxa")
colnames(a)[4] <- "NPQmax_2020"
a <- merge(a, blup2[,c("Taxa", "NPQmaxme")], by="Taxa")
colnames(a)[5] <- "NPQmax_2021"

a <- a %>% pivot_longer(cols = 4:5)
a <- a[!a$Study=="WiDiv-942**",]
a$year <- str_split_fixed(a$name, "_", 2)[,2]


ggplot(a, aes(Zm00001d042697, value)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  facet_grid(Study~year) + 
  xlab(expression(italic(PsbS)~(FPKM))) + 
  ylab(expression(NPQ[italic(max)]))
