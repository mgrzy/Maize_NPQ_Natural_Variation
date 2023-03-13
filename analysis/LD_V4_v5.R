library(rMVP)
library(data.table)

dat <- read.table("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.ind", header = F)
geno <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.desc")
map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = F)

dat2 <- read.table("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind", header = F)
geno2 <- attach.big.matrix("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.desc")
map2 <- fread("../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.map", data.table = F)

### PSBS

PSBSv4 <- data.frame(Taxa=dat$V1,
                     rs3_177340601=geno[which(map$SNP=="rs3_177340601"),])

PSBSv5 <- data.frame(Taxa=dat2$V1,
                     chr3_179566467=geno2[which(map2$SNP=="chr3_179566467"),])

b <- merge(PSBSv4, PSBSv5, by="Taxa")

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
                     chr7_182742836=geno2[which(map2$SNP=="chr7_182742836"),],
                     chr7_182742510=geno2[which(map2$SNP=="chr7_182742510"),])

acht <- merge(ACHT3v4, ACHTv5, by="Taxa")

# OEP37

OEPv4 <-  data.frame(Taxa=dat$V1,
                     rs5_187859720=geno[which(map$SNP=="rs5_187859720"),])

OEPv5 <- data.frame(Taxa=dat2$V1,
                    chr5_186953749=geno2[which(map2$SNP=="chr5_186953749"),],
                    chr5_186952808=geno2[which(map2$SNP=="chr5_186952808"),],
                    chr5_186952822=geno2[which(map2$SNP=="chr5_186952822"),],
                    chr5_186952823=geno2[which(map2$SNP=="chr5_186952823"),],
                    chr5_186955468=geno2[which(map2$SNP=="chr5_186955468"),],
                    chr5_186955834=geno2[which(map2$SNP=="chr5_186955834"),],
                    chr5_186955859=geno2[which(map2$SNP=="chr5_186955859"),])

oep <- merge(OEPv4, OEPv5)

### PMI1

PMIv4 <-  data.frame(Taxa=dat$V1,
                     rs1_86064177=geno[which(map$SNP=="rs1_86064177"),])

PMIv5 <-  data.frame(Taxa=dat2$V1,
                     chr1_85243338=geno2[which(map2$SNP=="chr1_85243338"),],
                     chr1_85243339=geno2[which(map2$SNP=="chr1_85243339"),],
                     chr1_85241837=geno2[which(map2$SNP=="chr1_85241837"),],
                     chr1_85241957=geno2[which(map2$SNP=="chr1_85241957"),],
                     chr1_85241996=geno2[which(map2$SNP=="chr1_85241996"),],
                     chr1_85242086=geno2[which(map2$SNP=="chr1_85242086"),],
                     chr1_85242272=geno2[which(map2$SNP=="chr1_85242272"),],
                     chr1_85242275=geno2[which(map2$SNP=="chr1_85242275"),],
                     chr1_85242380=geno2[which(map2$SNP=="chr1_85242380"),],
                     chr1_85242512=geno2[which(map2$SNP=="chr1_85242512"),],
                     chr1_85242515=geno2[which(map2$SNP=="chr1_85242515"),],
                     chr1_85242704=geno2[which(map2$SNP=="chr1_85242704"),],
                     chr1_85242710=geno2[which(map2$SNP=="chr1_85242710"),],
                     chr1_85241837=geno2[which(map2$SNP=="chr1_85241837"),],
                     chr1_85242734=geno2[which(map2$SNP=="chr1_85242734"),],
                     chr1_85242830=geno2[which(map2$SNP=="chr1_85242830"),],
                     chr1_85242986=geno2[which(map2$SNP=="chr1_85242986"),],
                     chr1_85243046=geno2[which(map2$SNP=="chr1_85243046"),],
                     chr1_85243088=geno2[which(map2$SNP=="chr1_85243088"),],
                     chr1_85243097=geno2[which(map2$SNP=="chr1_85243097"),],
                     chr1_85241837=geno2[which(map2$SNP=="chr1_85241837"),],
                     chr1_85243154=geno2[which(map2$SNP=="chr1_85243154"),],
                     chr1_85243178=geno2[which(map2$SNP=="chr1_85243178"),],
                     chr1_85243241=geno2[which(map2$SNP=="chr1_85243241"),],
                     chr1_85243304=geno2[which(map2$SNP=="chr1_85243304"),])

pm <- merge(PMIv4, PMIv5)
