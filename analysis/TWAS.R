library(data.table)
library(tidyverse)

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
