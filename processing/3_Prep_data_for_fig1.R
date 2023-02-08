library(tidyverse)
library(patchwork)
library(ggsci)
library(ggbeeswarm)
library(data.table)
library(readxl)


### Part one: Natural variaion in NPQ. 

npq <- read.csv("data/work/NpqRaw2020.csv")
traits <- read.csv("data/work/TraitsRaw2020.csv")

npq <- npq %>%
  filter(Taxa %in% c("PHGV6", "B73", "LH132", "PHP76"))

traits <- traits %>%
  filter(Taxa %in% c("PHGV6", "B73", "LH132", "PHP76"))

npq <- npq[,c(2,9:27)] # select columns

colnames(npq)[2:20] <- c(0.266333, 0.41333, 0.746667, 1.08, 2.08, 3.08, 4.08, 5.08 , 6.08, 7.08 , 8.08, 9.08, 
                         10.08, 10.26333, 10.41333, 10.91333, 11.91333, 15.08, 20.08) # change col names to minuets from NPQ protocol 

npq <- npq %>%
  pivot_longer(cols = 2:20) # transform to long format

npqSum <- npq %>%
  group_by(Taxa, name) %>%
  summarise(mean=mean(value), sd=sd(value)) # summarize NPQ data

traits <- traits %>%
  dplyr::select(c(2,9:38)) # select columns 

traits <- traits %>%
  pivot_longer(cols = 2:31) # transform to long format

traitsSum <- traits %>% 
  group_by(Taxa, name) %>%
  summarise(mean=mean(value), sd=sd(value)) # summarise trait data

traitsWide <- pivot_wider(traitsSum, 
                          id_cols = Taxa, names_from = name, values_from = mean)

datToPlotRaw <- traits %>%
  filter(Taxa %in% c("PHGV6", "B73", "LH132", "PHP76"),
         name %in% c("npqslop1me", "NPQmaxme", "npqslop2me", "NPQendme"))

write.csv(npq, "data/figures/Fig1a.csv", row.names = F)
write.csv(datToPlotRaw, "data/figures/Fig1b.csv", row.names = F)

#############################

map <- fread("../../BigData/WiDivGeno/RNAseqagpv4/rMVP/WiDiv752.geno.map", data.table = T)
RMIP <- fread("results/data/RMIP2020.csv", data.table = F)

RMIP <- RMIP[!RMIP$trait %in% c("qN_L1", "PC1", "PC2", "PC3"),]

RMIP$trait[RMIP$trait %in% c("af1npq", "af2npq", "npqslop1me")] <- "Induction"
RMIP$trait[RMIP$trait %in% c("bf1npq", "bf2npq", "NPQmaxme", "NPQ_L1")] <- "Steady-state"
RMIP$trait[RMIP$trait %in% c("af3npq", "af4npq", "npqslop2me")] <- "Relaxation"
RMIP$trait[RMIP$trait %in% c("bf3npq", "cf3npq", "bf4npq", "cf4npq", "NPQendme")] <- "Relaxed-stage"
RMIP$trait[RMIP$trait %in% c("phippsiiend", "af3phipsii", "bf3phipsii", "cf3phipsii", "af4phipsii", "bf4phipsii",
                             "cf4phipsii", "phipsiietos", "phtoNPQend", "QY_max", "Fv.Fm_L1", "QY_L1", "Rfd_L1")] <- "PSII"

RMIP <- merge(map[,1:3], RMIP, by="SNP", all.x = T)
RMIP <- as.data.frame(RMIP)

#RMIP$RMIP[RMIP$RMIP==1] <- NA
#RMIP$RMIP[RMIP$RMIP==2] <- NA

nCHR <- length(unique(RMIP$CHROM))
RMIP$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(RMIP$CHROM))){
  nbp[i] <- max(RMIP[RMIP$CHROM == i,]$POS)
  RMIP[RMIP$CHROM == i,"BPcum"] <- RMIP[RMIP$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

RMIP$trait <- factor(RMIP$trait, levels = c("Induction", "Steady-state", "Relaxation", "Relaxed-stage", 
                                            "PSII"))

greeks <- list("Induction", "Steady-state", "Relaxation", "Relaxed-stage", bquote(phi[PSII]))
RMIP2 <- na.omit(RMIP)
RMIP2$RMIP[RMIP2$RMIP<2] <- NA
RMIP2 <- na.omit(RMIP2)

write.csv(RMIP2, "data/figures/Fig1c.csv", row.names = F)
