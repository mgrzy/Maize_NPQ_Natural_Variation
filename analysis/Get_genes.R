library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)

map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = F)

txdb <- makeTxDbFromGFF("../../BigData/Maize/v4/GFF/Zea_mays.AGPv4.36.gff3.gz")

###
RMIP <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/RMIP/RMIP2021.csv", data.table = F)
RMIP <- RMIP[RMIP$RMIP>4,]
RMIP <- merge(RMIP, map, by="SNP", all.y=F)

###
pos <- data.frame(chrom=RMIP$CHROM, start=RMIP$POS, end=RMIP$POS, trait=RMIP$trait, RMIP=RMIP$RMIP, pos=RMIP$POS)
pos$end <- pos$end+100000

###
db <- makeGRangesFromDataFrame(pos, seqnames.field = "chrom", start.field = "start", end.field = "end", keep.extra.columns = T)

gen <- transcripts(txdb)

res <- subsetByOverlaps(gen, db)
res <- mergeByOverlaps(res, db)
res <- as.data.frame(res)
res <- res[,c("res.seqnames", "res.start", "res.end", "res.tx_name", "trait", "RMIP", "pos")]
colnames(res)[4] <- c("Gen")
#res$Gen <- gsub("transcript:", "", res$Gen)
res$Gen <-  stringr::str_split_fixed(res$Gen, "_", 2)[,1]
#res <- res[!duplicated(res$Gen),]
rm(gen)

