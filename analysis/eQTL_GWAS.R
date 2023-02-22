library(rMVP)

dat <- read.csv("../Data/Maize_NPQ_Natural_Variation/data/work/exp_eqt_rnav4.csv")

geno <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.desc")
kin <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.kin.desc")
pc <- attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.pc.desc")[]
map <- fread("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.map", data.table = F)

co <- model.matrix(~dat$Study)
co <- co[,2:3]

pc <- cbind(pc, co)

res <- map
for(j in 7:13){
  imMVP <- MVP(
    phe=dat[, c(1, j)],
    geno=geno,
    map=map,
    K=kin,
    CV.MLM=pc,
    priority="speed",
    ncpus=10,
    vc.method="BRENT",
    method=c("MLM"),
    file.output=F
  )
  res <- cbind.data.frame(res, imMVP$mlm.results)
  gc()
}

fwrite(res, "../Data/Maize_NPQ_Natural_Variation/bigResults/eQTL_rnav4.csv.gz")
