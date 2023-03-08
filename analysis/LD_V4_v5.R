library(rMVP)
library(data.table)

dat <- read.table("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.ind", header = F)
geno <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.desc")
map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = F)

dat2 <- read.table("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind", header = F)
geno2 <- attach.big.matrix("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.desc")
map2 <- fread("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.map", data.table = F)

### PSI3

PSI3v4 <- data.frame(Taxa=dat$V1,
                   rs3_176152936=geno[which(map$SNP=="rs3_176152936"),])

PSI3v5 <- data.frame(Taxa=dat2$V1,
                     chr3_178405582=geno2[which(map2$SNP=="chr3_178405582"),],
                     chr3_178405643=geno2[which(map2$SNP=="chr3_178405643"),],
                     chr3_178407540=geno2[which(map2$SNP=="chr3_178407540"),],
                     chr3_178407565=geno2[which(map2$SNP=="chr3_178407565"),])

a <- merge(PSI3v4, PSI3v5, by="Taxa")
