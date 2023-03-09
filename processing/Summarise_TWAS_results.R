library(data.table)
library(tidyverse)

dat20 <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2020/TWAS_WiDiv503.csv.gz")

dat20 <- dat20[,-c(2:5)]

dat20 <- dat20 %>% pivot_longer(2:32)

dat20 %>% 
  filter(value<1.826684e-06)

dat21 <- fread("../Data/Maize_NPQ_Natural_Variation/bigResults/TWAS/2021/TWAS_WiDiv503_2021.csv.gz")

dat21 <- dat21[,-c(2:5)]

dat21 <- dat21 %>% pivot_longer(2:32)

a21 <- dat21 %>% 
  filter(value<1.826684e-06)

a20 <- dat20 %>% 
  filter(value<1.826684e-06)

a21$year <- 2021
a20$year <- 2020

a <- rbind.data.frame(a20, a21)

write.csv(a, "../Data/Maize_NPQ_Natural_Variation/results/data/TWAS.csv", row.names = F)
