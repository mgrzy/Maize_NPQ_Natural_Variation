library(tidyverse)
library(data.table)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

NPQ20 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2020gwasRNA.csv")
NPQ21 <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/blups/blup2021gwasRNA.csv")
NPQ20 <- NPQ20[,-c(2:3)]
NPQ21 <- NPQ21[,-c(2:3)]

Bugeater20 <- readxl::read_xls("../Data/Maize_NPQ_Natural_Variation/data/meta/Bugeater2020_merged_v5.xls")
Bugeater21 <- readxl::read_xlsx("../Data/Maize_NPQ_Natural_Variation/data/meta/Bugeater2021_merged_v6.2.xlsx")

Bugeater20 <- Bugeater20[,-c(1:4)]
colnames(Bugeater20)[1] <- "Taxa"

Bugeater21 <- Bugeater21[,-c(1,2,3,4,6)]
colnames(Bugeater21)[1] <- "Taxa"

Bugeater20 <- Bugeater20[,-c(7:8)]

Bugeater20 <- Bugeater20 %>%
  pivot_longer(2:31) %>%
  group_by(Taxa, name) %>%
  summarise(value=mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = Taxa, names_from = name, values_from = value)

Bugeater21 <- Bugeater21[,-c(4,37)]

Bugeater21 <- Bugeater21 %>%
  pivot_longer(2:31) %>%
  group_by(Taxa, name) %>%
  summarise(value=mean(value, na.rm=T)) %>%
  pivot_wider(id_cols = Taxa, names_from = name, values_from = value)

dat20 <- merge(NPQ20, Bugeater20, by="Taxa")

traits <- colnames(NPQ20)[2:32]

dat20 <- as.data.frame(dat20)
res20 <- c()

for (i in traits){
  z <- cor(dat20[,i], dat20[,33:62], use = "complete.obs")
  res20 <- rbind(res20, z)
  rm(z)
}
res20 <- as.data.frame(res20)
res20$trait <- traits

a <- res20 %>%
  pivot_longer(1:30)

dat21 <- merge(NPQ21, Bugeater21, by="Taxa")

dat21 <- as.data.frame(dat21)
res21 <- c()

for (i in traits){
  z <- cor(dat21[,i], dat21[,33:62], use = "complete.obs")
  res21 <- rbind(res21, z)
  rm(z)
}
res21 <- as.data.frame(res21)
res21$trait <- traits

a21 <- res21 %>%
  pivot_longer(1:30)
a$year <- 2020
a21$year <- 2021

a <- rbind.data.frame(a, a21)
a <- a %>%
  filter(year==2020)

ggplot(a, aes(x=value)) +
  geom_histogram(binwidth = 0.025) +
  xlab("Pearson correlation")

ggsave("../Data/Maize_NPQ_Natural_Variation/figures/corr.png", units = "mm", width = 75, height = 40, scale = 2)

