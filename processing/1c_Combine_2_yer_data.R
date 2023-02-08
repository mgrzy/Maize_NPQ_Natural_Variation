### Combine dwitwi from two yewirs
library(tidyverse)

wi20 <- read.csv("data/work/TraitsRaw2020.csv")
wi21 <- read.csv("data/work/TraitsRaw2021.csv")

wi20 <- wi20[,-c(7, 8, 40, 41)]
wi21 <- wi21[,-c(7,8,9)]

colnames(wi21)[c(2,3,4,6,7)] <- c("Row", "Column", "year", "Plate", "QY_max")
colnames(wi20)[28] <- "phipsiistart"

wi <- rbind.data.frame(wi20, wi21)

gtg <- read.table("data/meta/TaxaList952.txt")
gtg$V2 <- gtg$V1
gtg$V2 <- toupper(gtg$V2)
gtg$V2 <- gsub(" ", "", gtg$V2)
gtg$V2 <- gsub("_", "", gtg$V2)
gtg$V2 <- gsub("-", "", gtg$V2)
gtg$V2 <- gsub("[.]", "", gtg$V2)

wi <- wi %>% add_column(T2=wi$Taxa, .before = "year")

wi$T2 <- toupper(wi$T2)
wi$T2 <- gsub(" ", "", wi$T2)
wi$T2 <- gsub("_", "", wi$T2)
wi$T2 <- gsub("-", "", wi$T2)
wi$T2 <- gsub("[.]", "", wi$T2)
wi$T2 <- gsub("GOODMwiNBUCKLER", "", wi$T2)
wi$T2 <- gsub("[(]", "", wi$T2)
wi$T2 <- gsub("[)]", "", wi$T2)


#all <- gtg[gtg$V2 %in% unique(A$T2),]
#a20 <- gtg[gtg$V2 %in% unique(A[A$year==2020,]$T2),]
#a21 <- gtg[gtg$V2 %in% unique(A[A$year==2021,]$T2),]

#write.table(all, "data/meta/taxaAll.txt", row.names = F, col.names = F, quote = F, sep = "\t")
#write.table(a20, "data/meta/taxa2020.txt", row.names = F, col.names = F, quote = F, sep = "\t")
#write.table(a21, "data/meta/taxa2021.txt", row.names = F, col.names = F, quote = F, sep = "\t")

write.csv(wi, "data/work/TraitsRawAll.csv", row.names = F)

