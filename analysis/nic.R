library(tidyverse)
library(data.table)

map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = T)
RMIP <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/RMIP/RMIP2021.csv", data.table = F)

RMIP$trait[RMIP$trait %in% c("af1npq", "af2npq", "npqslop1me")] <- "Induction"
RMIP$trait[RMIP$trait %in% c("bf1npq", "bf2npq", "NPQmaxme", "NPQ_L1")] <- "Steady-state"
RMIP$trait[RMIP$trait %in% c("af3npq", "af4npq", "npqslop2me")] <- "Relaxation"
RMIP$trait[RMIP$trait %in% c("bf3npq", "cf3npq", "bf4npq", "cf4npq", "NPQendme")] <- "Relaxed-stage"
RMIP$trait[RMIP$trait %in% c("phippsiiend", "af3phipsii", "bf3phipsii", "cf3phipsii", "af4phipsii", "bf4phipsii",
                             "cf4phipsii", "phipsiietos", "phtoNPQend", "QY_max", "Fv.Fm_L1", "QY_L1", "Rfd_L1")] <- "PSII"

RMIP <- merge(map[,1:3], RMIP, by="SNP", all.x = T)
RMIP <- as.data.frame(RMIP)

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

write.csv(RMIP2, "../Data/Maize_NPQ_Natural_Variation/data/work/RMIPclean2021.csv", row.names = F)

####