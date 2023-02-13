library(data.table)
library(tidyverse)

blup <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2020gwasRNA.csv")
blup2 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2021gwasRNA.csv")
ex <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Expression_MetaData_csv.gz")
gene <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Gene_MetaData_csv.gz")

ex.1 <- ex[ex$Study=="WiDiv-503", -c(2,3)]
ex.2 <- ex[ex$Study=="WiDiv-942", -c(2:3)]

ex.1.pc <- ex.1[,1:6]
ex.1 <- ex.1[,-c(2:6)]

ex.2.pc <- ex.2[,1:6]
ex.2 <- ex.2[,-c(2:6)]

pheno <- blup[,-c(2,3)]

X <- merge(pheno, ex.1, by="Taxa")
X <- as.data.frame(X)

toRemove <- which(as.vector(is.na(X[,2])))

X <- X[-c(toRemove),]  
ex.1.pc <- ex.1.pc[-c(toRemove),]

genExp <- X[,-c(1:33)]

genToKeep <- apply(genExp, 2, FUN = function(x) {
  b <- table(x > 0)
  z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
  return(z)}
)

genToKeep <- which(genToKeep==TRUE)
genExp <- genExp[,genToKeep]
genExp <- as.data.frame(genExp)

ex.1.pc <- as.matrix(ex.1.pc[,-c(1)])

res <- gene
res <- res[res$gene %in% colnames(genExp)]

for (i in 2:32){
  ph <- X[,i]

  a <- apply(genExp, 2, function(x) 
    summary(lm(ph ~ ex.1.pc + x))$coefficients[7,4])

  a <- data.frame(gene=names(a), p_val=a)
  colnames(a)[2] <- colnames(X)[i]

  res <- merge(res, a, by="gene")
  rm(a)
  print(paste(colnames(X)[i]))
}

######
###
######

X <- merge(pheno, ex.2, by="Taxa")
X <- as.data.frame(X)

toRemove <- which(as.vector(is.na(X[,2])))

X <- X[-c(toRemove),]  
ex.2.pc <- ex.2.pc[-c(toRemove),]

genExp <- X[,-c(1:33)]

genToKeep <- apply(genExp, 2, FUN = function(x) {
  b <- table(x > 0)
  z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
  return(z)}
)

genToKeep <- which(genToKeep==TRUE)
genExp <- genExp[,genToKeep]
genExp <- as.data.frame(genExp)

ex.2.pc <- as.matrix(ex.2.pc[,-c(1)])

res2 <- gene
res2 <- res2[res2$gene %in% colnames(genExp)]

for (i in 2:32){
  ph <- X[,i]
  
  a <- apply(genExp, 2, function(x) 
    summary(lm(ph ~ ex.2.pc + x))$coefficients[7,4])
  
  a <- data.frame(gene=names(a), p_val=a)
  colnames(a)[2] <- colnames(X)[i]
  
  res2 <- merge(res2, a, by="gene")
  rm(a)
  print(paste(colnames(X)[i]))
}

###

pheno <- blup2[,-c(2,3)]

X <- merge(pheno, ex.1, by="Taxa")
X <- as.data.frame(X)

toRemove <- which(as.vector(is.na(X[,2])))

X <- X[-c(toRemove),]  
ex.1.pc <- ex.1.pc[-c(toRemove),]

genExp <- X[,-c(1:33)]

genToKeep <- apply(genExp, 2, FUN = function(x) {
  b <- table(x > 0)
  z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
  return(z)}
)

genToKeep <- which(genToKeep==TRUE)
genExp <- genExp[,genToKeep]
genExp <- as.data.frame(genExp)

ex.1.pc <- as.matrix(ex.1.pc[,-c(1)])

res3 <- gene
res3 <- res3[res3$gene %in% colnames(genExp)]

for (i in 2:32){
  ph <- X[,i]
  
  a <- apply(genExp, 2, function(x) 
    summary(lm(ph ~ ex.1.pc + x))$coefficients[7,4])
  
  a <- data.frame(gene=names(a), p_val=a)
  colnames(a)[2] <- colnames(X)[i]
  
  res3 <- merge(res3, a, by="gene")
  rm(a)
  print(paste(colnames(X)[i]))
}

###

pheno <- blup2[,-c(2,3)]

X <- merge(pheno, ex.2, by="Taxa")
X <- as.data.frame(X)

toRemove <- which(as.vector(is.na(X[,2])))

X <- X[-c(toRemove),]  
ex.2.pc <- ex.2.pc[-c(toRemove),]

genExp <- X[,-c(1:33)]

genToKeep <- apply(genExp, 2, FUN = function(x) {
  b <- table(x > 0)
  z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
  return(z)}
)

genToKeep <- which(genToKeep==TRUE)
genExp <- genExp[,genToKeep]
genExp <- as.data.frame(genExp)

ex.2.pc <- as.matrix(ex.2.pc[,-c(1)])

res4 <- gene
res4 <- res4[res4$gene %in% colnames(genExp)]

for (i in 2:32){
  ph <- X[,i]
  
  a <- apply(genExp, 2, function(x) 
    summary(lm(ph ~ ex.2.pc + x))$coefficients[7,4])
  
  a <- data.frame(gene=names(a), p_val=a)
  colnames(a)[2] <- colnames(X)[i]
  
  res4 <- merge(res4, a, by="gene")
  rm(a)
  print(paste(colnames(X)[i]))
}

###

fwrite(res, "../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv503.csv.gz")
fwrite(res2, "../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv942.csv.gz")
fwrite(res3, "../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2021/TWAS_WiDiv503_2021.csv.gz")
fwrite(res4, "../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2021/TWAS_WiDiv942_2021.csv.gz")

