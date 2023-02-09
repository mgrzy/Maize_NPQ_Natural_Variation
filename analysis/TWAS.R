library(data.table)
library(tidyverse)
library(sommer)

K <- as.matrix(bigmemory::attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.kin.desc"))
ind <- read.table("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.ind")

blup <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2020gwasRNA.csv")
ex <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Expression_MetaData_csv.gz")
gene <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Gene_MetaData_csv.gz")

ex.1 <- ex[ex$Study=="WiDiv-503", -c(2,3)]
ex.2 <- ex[ex$Study=="WiDiv-942", -c(2:3)]

ex.1.pc <- ex.1[,1:6]
ex.1 <- ex.1[,-c(2:6)]

ex.2.pc <- ex.1[,1:6]
ex.2 <- ex.1[,-c(1:6)]

pheno <- blup[,c(1,11)]

phenotype <- pheno
expressionMatrix <- ex.1
covariate <- ex.1.pc

commonName <- colnames(phenotype)[1]

X <- merge(phenotype, expressionMatrix, by=commonName)
X <- as.data.frame(X)

toRemove <- which(as.vector(is.na(X[,2])))

X <- X[-c(toRemove),]  
covariate <- covariate[-c(toRemove),]

ph <- as.numeric(X[,2])

genExp <- X[,-c(1,2)]

genToKeep <- apply(genExp, 2, FUN = function(x) {
  b <- table(x > 0)
  z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
  return(z)}
)

genToKeep <- which(genToKeep==TRUE)
genExp <- genExp[,genToKeep]
genExp <- as.data.frame(genExp)

covariate <- as.matrix(covariate[,-c(1)])

res <- apply(genExp, 2, function(x) 
  summary(lm(ph ~ covariate + x))$coefficients[7,4])

res <- data.frame(gene=names(res), p_val=res)

res <- merge(res, gene)

rownames(K) <- ind$V1
colnames(K) <- ind$V1

z <- data.frame(ph=ph, covariate, genExp)
z <- data.frame(Taxa=X$Taxa, z)

m <- mmer(ph ~  PC1 + PC2 + PC3,
          random=~vsr(Taxa,Gu=K),
          rcov=~units, nIters=3,
          data=z, verbose = FALSE)

resMix <- data.frame(gen=c(), pmlm=c())

for (i in 25001:ncol(z)){
  
  q <- z[,c(1:5,i)]
  
  g <- colnames(q)[6]
  
  colnames(q)[6] <- 'gene'
  
  m2 <- mmer(ph ~  PC1 + PC2 + PC3 + gene,
             random=~vsr(Taxa,Gu=K),
             rcov=~units, nIters=3,
             data=q, verbose = FALSE)
  a <- anova(m, m2)
  aa <- data.frame(gene=g, pmlm=a$PrChisq[2])
  resMix <- rbind.data.frame(resMix, aa)
}

write.csv(res, "GLM_TWAS_1.csv", row.names = F)
write.csv(resMix, "MLM_TWAS_6.csv", row.names = F)

