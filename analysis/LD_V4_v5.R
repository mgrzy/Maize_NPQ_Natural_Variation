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

### IRM1
IRM1v4 <- data.frame(Taxa=dat$V1,
                     rs2_182382835=geno[which(map$SNP=="rs2_182382835"),])

IRM1v5 <- data.frame(Taxa=dat2$V1,
                     chr2_181981544=geno2[which(map2$SNP=="chr2_181981544"),],
                     chr2_181981600=geno2[which(map2$SNP=="chr2_181981600"),],
                     chr2_181981652=geno2[which(map2$SNP=="chr2_181981652"),],
                     chr2_181982035=geno2[which(map2$SNP=="chr2_181982035"),])

irm1 <- merge(IRM1v4, IRM1v5, by="Taxa")

###ACHT3

ACHT3v4 <- data.frame(Taxa=dat$V1,
                     rs3_176152936=geno[which(map$SNP=="rs7_179281748"),])

ACHTv5 <- data.frame(Taxa=dat2$V1,
                     chr7_182742647=geno2[which(map2$SNP=="chr7_182742647"),],
                     chr7_182742793=geno2[which(map2$SNP=="chr7_182742793"),],
                     chr7_182742814=geno2[which(map2$SNP=="chr7_182742814"),],
                     chr7_182742836=geno2[which(map2$SNP=="chr7_182742836"),])

acht <- merge(ACHT3v4, ACHTv5, by="Taxa")

###
