library(tidyverse)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

h2 <- read.csv("../Data/Maize_NPQ_Natural_Variation/results/data/Broad_Heritability.csv")
h2 <- h2 %>% 
  filter(year %in% c("2020", "2021"))


ggplot(h2, aes(reorder(trait, H2), H2, fill=year)) + 
  geom_col(stat="identity", position = position_dodge()) + 
  coord_flip() + 
  xlab("trait") + 
  theme(legend.position = "top")

h2$H2 <- round(h2$H2, 2)

write.csv(h2, "H2.csv", row.names = F)

h2 %>%
  group_by(year) %>%
  summarise(mean(H2))
