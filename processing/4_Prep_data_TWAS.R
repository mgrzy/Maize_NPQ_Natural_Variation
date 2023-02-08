library(data.table)
library(tidyverse)

blup <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2020gwasRNA.csv")
ex <- fread("../../BigData/Maize/v4/WiDivExpression/Wisconsin942_Expression.csv")
ex.meta2 <- fread("../../BigData/Maize/MetaData/WiDiv_Meta_PCA.csv")

### Subset meta data from expression table
ex.meta <- ex[,c(1:5)]
ex.Taxa <- colnames(ex)[-c(1:5)]
ex.Taxa <- gsub(" ", "_", ex.Taxa)

ex <- data.table::transpose(ex[,-c(1:5)])

ex <- data.frame(Taxa=ex.Taxa, ex)
colnames(ex)[2:ncol(ex)] <- ex.meta$gene

x <- merge(ex.meta2, ex, by="Taxa")

fwrite(x, "../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Expression_MetaData_csv.gz")
fwrite(ex.meta, "../Data/Maize_NPQ_Natural_Variation/data/work/expression/WiDiv_Gene_MetaData_csv.gz")
