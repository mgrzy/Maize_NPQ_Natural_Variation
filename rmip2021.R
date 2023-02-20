library(data.table)
library(tidyverse)

# Set them to classic and axis text to black rather then grey
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

RMIP <- fread("../Data/Maize_NPQ_Natural_Variation/data/work/RMIPclean2021.csv")

RMIP$BPcum <- as.numeric(RMIP$BPcum)

axis.set <- RMIP %>%
  group_by(CHROM) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

RMIP$trait <- factor(RMIP$trait, levels = c("Induction", "Steady-state", "Relaxation", "Relaxed-stage", 
                                            "PSII", "PSIIrate"))

greeks <- list("Induction", "Steady-state", "Relaxation", "Relaxed-stage", bquote(phi[PSII]), bquote(phi[PSII_2]))

RMIP <- RMIP[RMIP$RMIP>1,]

ggplot(RMIP, aes(x = BPcum, y = RMIP/100)) +
  #geom_segment(data=gtm, aes(x=BPcum, xend=BPcum, y=0, yend=po)) + # black line for marking genes position
  #geom_text(data=gtm, aes(BPcum, y=po, label=Gen))  + # genes names
  #annotate("text", x = 727487348, y = 0.9, label = "italic(PsbS)", parse=TRUE, size=6) + # PsbS name
  #geom_vline(xintercept = 727487348, linetype=2, alpha=0.3) + # PsbS dashed line position
  geom_point(size=5, alpha=0.8, aes(colour=trait)) + # points for RMIP values
  geom_hline(yintercept = c(0.05, 0.2), color = c("grey40", "darkred"), linetype = "dashed") + # dashed significance thresholds
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + # chromosome numbers in the middle 
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.75)) + # y scale limits; do not change expand values
  labs(x = "Chromosome",  y = "RMIP") + 
  theme(legend.position = c(0.95, 0.85)) + 
  scale_color_manual(values = c("yellow", "orange", "darkgreen", "green",
                                        "grey", "grey10"),
                                        name="Trait related to", na.translate = F, label=greeks) # colour values for each trait
