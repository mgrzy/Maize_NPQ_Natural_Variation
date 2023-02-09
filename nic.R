w1 <- read.csv("MLM_TWAS_1.csv")
w2 <- read.csv("MLM_TWAS_2.csv")
w3 <- read.csv("MLM_TWAS_3.csv")
w4 <- read.csv("MLM_TWAS_4.csv")
w5 <- read.csv("MLM_TWAS_5.csv")
w6 <- read.csv("MLM_TWAS_6.csv")

w <- rbind.data.frame(w1, w3, w4, w5, w6)
colnames(w)[2] <- "MLM"

res <- read.csv("GLM_TWAS_1.csv")

w$MLM <- gsub("*", "" ,w$MLM)
w$MLM <- gsub("[**]", "" ,w$MLM)
w$MLM <- gsub("*", "" ,w$MLM)
w$MLM <- str_replace_all(string=w$MLM, pattern=" .", repl="")
w$MLM2 <- as.numeric(w$MLM)
wg3$mlm2 <- w$m

res <- res[res$chromosome %in% c(1:10),]
res$chromosome <- as.numeric(res$chromosome)

res <- rbind.data.frame(res[!res$chromosome==10,], res[res$chromosome==10,])

nCHR <- length(unique(res$chromosome))
res$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(res$chromosome))){
  nbp[i] <- max(res[res$chromosome == i,]$position_left)
  res[res$chromosome == i,"BPcum"] <- res[res$chromosome == i,"position_left"] + s
  s <- s + nbp[i]
}

axis.set <- res %>% 
  group_by(chromosome) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

axis.set$chromosome <- factor(axis.set$chromosome, levels = c(1:10))

res$mlm2 <- w$m

res2 <- merge(res, w, by="gene")
theme_set(theme_classic(base_size = 20, base_family = "NimbusSan"))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

res2$chromosome <- factor(res2$chromosome, levels = 1:10)

g1 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res)) , linetype=2) +
  geom_point(data=res2, aes(BPcum, -log10(p_val), colour=chromosome)) + 
  #scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome") + 
  ylim(0, 18)

g2 <- ggplot() + 
  geom_hline(yintercept = -log10(0.05/nrow(res)) , linetype=2) +
  geom_point(data=res2, aes(BPcum, -log10(MLM2), colour=chromosome)) + 
  #scale_color_d3() + 
  scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab("-log10(p-value)") +  
  xlab("Chromosome") + 
  ylim(c(0,18))

g3 <- ggplot() + 
  geom_point(data=res2, aes(
  y = -log10(sort(p_val, decreasing = F)), 
  x = -log10(ppoints(length(p_val))))) + 
  geom_point(data=res2, aes(
    y = -log10(sort(MLM2, decreasing = F)), 
    x = -log10(ppoints(length(p_val)))), colour="red") +
  geom_abline(slope = 1) +
  xlim(0, max(-log10(res2$MLM2) + 0.5)) +
  ylim(0, max(-log10(res2$MLM2) + 0.5)) + 
  ylab("Observed -log10(p)") + 
  xlab("Expected -log10(p)")
