library(data.table)
library(tidyverse)

ex <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Expression_MetaData_csv.gz")

genes <- c("Zm00001d042697", "Zm00001d042669" ,"Zm00001d005657", "Zm00001d042017", "Zm00001d022518", "Zm00001d029761", 
           "Zm00001d017171")

dat <- ex %>%
  select(c("Taxa", "Study", "PC1", "PC3", "PC3", "PC4", "PC5", genes))

tx.rna <- read.table("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.ind")
colnames(tx.rna) <- "Taxa"
tx.wgs <- read.table("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind")
colnames(tx.wgs) <- "Taxa"

dat.rna <- plyr::join(tx.rna, dat, by="Taxa")

dat.wgs <- plyr::join(tx.wgs, dat, by="Taxa")

write.csv(dat.rna, "../Data/Maize_NPQ_Natural_Variation/data/work/exp_eqt_rnav4.csv", row.names = F)
write.csv(dat.wgs, "../Data/Maize_NPQ_Natural_Variation/data/work/exp_eqt_wgsv5.csv", row.names = F)