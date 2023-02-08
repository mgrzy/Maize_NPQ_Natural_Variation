###########################
# Data preparation 2021 ###
###########################
library(readxl)
library(stringr)
library(tidyverse)
###########################

### Checks 
chrm <- read.csv("data/meta/checksToRemove2021.csv")
colnames(chrm) <- "Plot"

### Field map data
pl <- read.csv("data/meta/PlotMap2021.csv") # Filed map
pl <- pl[,3:5]
colnames(pl)[1] <- "Plot"
pl <- pl[!pl$Plot %in% chrm$Plot,]

### NPQ raw data
### Format data
A <- read_xlsx("data/raw/Wisconsin_field_2021_GWAS_datawith_changes_2022.xlsx", sheet = 2, na = ".")
A <- A[,-c(3)]
colnames(A)[4] <- "Plot"
colnames(A)[2] <- "Taxa"

#### Reomve missing data
#A <- A[1:5214,] 
#A <- A[!is.na(A$),]

### Remove checks 
A <- A[!A$Plot %in% chrm$Plot,]


A2 <- A[A$FILTER0615==0,]

keep <- colnames(A)[c(1,2,4,5,6,7,8, 111, 112, 131, 150, 169, 188, 226, 239:268)] # columns to keep
npq <- colnames(A)[c(1,2,4,5,6, 150:168)] # columns to keep
f <- colnames(A)[c(2,4,5,6,7,9,131:149)]

Af <- A[,f]
Anpq <- A[,npq]  # keep only selected columns 
A <- A[,keep] # keep only selected columns 

A <- A[,-c(17, 20, 24, 28, 38, 42)]

A <- merge(pl, A, by="Plot") # add filed columns and rows 
Anpq <- merge(pl, Anpq, by="Plot") # add filed columns and rows 
Af <- merge(pl, Af, by="Plot")

plate <- str_split_fixed(A$well, "", 2) # split 
plate <- as.data.frame(plate)

lookup <- setNames(seq_along(LETTERS), LETTERS)
plate$V1 <- as.numeric(lookup[plate$V1])
colnames(plate) <- c("Row_plate", "Column_plate")

write.csv(A,"data/work/TraitsRaw2021.csv", row.names = F)
write.csv(Anpq,"data/work/NpqRaw2021.csv", row.names = F)
write.csv(Af,"data/work/Qy_Raw2021.csv", row.names = F)
