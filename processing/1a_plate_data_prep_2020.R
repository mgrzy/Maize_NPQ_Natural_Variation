############################
# Data preparation 2020 ###
###########################
library(readxl)
library(stringr)
library(tidyverse)
############################

sc <- c(0.266333, 0.41333, 0.746667, 1.08, 2.08, 3.08, 4.08, 5.08 , 6.08, 7.08 , 8.08, 9.08, 
        10.08, 10.26333, 10.41333, 10.91333, 11.91333, 15.08, 20.08)

A <- read_xlsx("data/raw/Wisconsin_field_2020_GWAS datawith changes 2022.xlsx", sheet = 2)

pl <- read.csv("data/meta/PlotMap2020.csv") # Filed map
colnames(pl)[1] <- "Plot"

A <- A[1:4480,] #reomve missing data

A <- A[A$`FILTER fvfm065p1`==0,] # remove bad data

colnames(A)[2] <- c("Plot")

#data[data$Taxa=="CI_91B_Goodman-Buckler",]$Taxa <- "CI_91B"

keep <- colnames(A)[c(2,4,5,6,7,9, 112, 113, 132, 151, 170, 189, 227, 240:245, 247, 248, 250, 251, 252, 254:256, 258:262, 264:266, 268, 269)] # columns to keep
npq <- colnames(A)[c(2,4,5,6,7,9, 151:169)] # columns to keep
f <- colnames(A)[c(2,4,5,6,7,9,132:150)]

Af <- A[,f]
Anpq <- A[,npq]  # keep only selected columns 
A <- A[,keep] # keep only selected columns 
#colnames(A)[28] <- c("phippsiiend")

A <- merge(pl, A, by="Plot") # add filed columns and rows 
Anpq <- merge(pl, Anpq, by="Plot") # add filed columns and rows 
Af <- merge(pl, Af, by="Plot")

plate <- str_split_fixed(A$Well, "", 2) # split 
plate <- as.data.frame(plate)

lookup <- setNames(seq_along(LETTERS), LETTERS)
plate$V1 <- as.numeric(lookup[plate$V1])
colnames(plate) <- c("Row_plate", "Column_plate")

taxa <- read.csv("data/meta/TaxaPlot.csv")
colnames(taxa)[2] <- "Plot"

A <- merge(taxa, A, by="Plot")
Anpq <- merge(taxa, Anpq, by="Plot")
Af <- merge(taxa, Af, by="Plot")

A <- cbind(A, plate)
Anpq <- cbind(Anpq, plate)
Af <- cbind(Af, plate)

A <- A[,-c(3,6,7)]
Anpq <- Anpq[,-c(3,6,7)]
Af <- Af[,-c(3,6,7)]

A <- A %>%
  add_column(year=2020, .before = "Row")

Anpq <- Anpq %>%
  add_column(year=2020, .before = "Row")

Af <- Af %>%
  add_column(year=2020, .before = "Row")

write.csv(A,"data/work/TraitsRaw2020.csv", row.names = F)
write.csv(Anpq,"data/work/NpqRaw2020.csv", row.names = F)
write.csv(Af,"data/work/Qy_Raw2020.csv", row.names = F)
