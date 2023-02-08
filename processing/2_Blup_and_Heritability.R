library(lme4)
library(bestNormalize)
library(tidyverse)

A <- read.csv("data/work/TraitsRawAll.csv")

A <- A %>%
  add_column(.after = "Plot", 
             Block=str_split_fixed(A$Plot, "", 4)[,1])

trait <- colnames(A)[9:39]
taxa <- unique(A$T2)
taxa2020 <- unique(A[A$year==2020,]$T2)
taxa2021 <- unique(A[A$year==2021,]$T2)

blupAll <- data.frame(Taxa=taxa)
blupAll.TrueScale <- blupAll

blup2020 <- data.frame(Taxa=taxa2020)
blup.TrueScale2020 <- blup2020

blup2021 <- data.frame(Taxa=taxa2021)
blup.TrueScale2021 <- blup2021

block.blup <- data.frame(Taxa=as.character(), 
                         b1=as.numeric(), 
                         b2=as.numeric(), 
                         trait=as.character(), 
                         year=as.character())

block.cor <- data.frame(cor=as.numeric(), 
                        trait=as.character(), 
                        year=as.character())

H2 <- data.frame(trait=as.character(), 
                 year=as.character(), 
                 H2=as.numeric())

P <- data.frame(trait=as.character(), 
                Pval=as.numeric(),
                AICF=as.numeric(), 
                AICS=as.numeric(),
                LogLikF=as.numeric(), 
                LogLikS=as.numeric())

for (i in 1:31){

  ###########################
  ### Model for full data ###
  ###########################
  
  ### Subset 
  B <- A[,c("T2", "year", "Plate",trait[i])]
  
  #Find and remove outlier
  B$out <- rstatix::is_outlier(B[,4])
  B <- B[B$out==FALSE,]

  ### Transform and normalize data
  bn <- bestNormalize(B[,4], standardize = T)
  B$y <- bn$x.t

  ### Fit model with lme4
  mF <- lmer(y ~ (1|T2) + (1|year) + (1|Plate:year), B)

  ### Calculate heritability with "standard method"
  v <- as.data.frame(VarCorr(mF))
  h <- v$vcov[1] / ((v$vcov[1]) + (v$vcov[2]/3) + (v$vcov[3]/2) + (v$vcov[4]/6))
  
  ### Merge H2 into one data frame
  h <- data.frame(trait=trait[i],
                  year="All", 
                  H2=h)
  H2 <- rbind(H2, h)
  
  ### Obtain blups plus coefficinet and merge to one data frame
  b <- data.frame(Taxa=rownames(coef(mF)$T2),
                  a=coef(mF)$T2[,1])
  colnames(b)[2] <- trait[i]
  blupAll <- plyr::join(blupAll, b, by="Taxa")
  
  ### Back transform blups to oryginal unit and save for visualtization 
  b[,2] <- predict(bn, b[,2], inverse=T)
  blupAll.TrueScale <- plyr::join(blupAll.TrueScale, b)
  
  ### Test significance of year
  mF <- lmer(y ~ year + (1|T2) + (1|Plate:year), B)
  mS <- lmer(y ~  (1|T2) + (1|Plate:year), B)
  An <- anova(mF, mS)
  
  An <- data.frame(trait=trait[i], 
                   P=An$`Pr(>Chisq)`[2], 
                   AICF=An$AIC[2], 
                   AICS=An$AIC[1], 
                   LogLikF=An$logLik[1], 
                   LogLikS=An$logLik[2])
  
  P <- rbind(P, An)
  

  ###########################
  ### Model for 2020 data ###
  ###########################
  
  ## Subset 
  B <- A[A$year=="2020", c("T2", "year", "Plate","Block",trait[i])]
  
  ### Find and remove outlier
  B$out <- rstatix::is_outlier(B[,5])
  B <- B[B$out==FALSE,]
  
  ### Transform and normalize data
  bn <- bestNormalize(B[,5], standardize = T)
  B$y <- bn$x.t
  
  ### Fit model with lme4
  mF <- lmer(y ~ (1|T2) + (1|Plate), B)
  mF1 <- lmer(y ~ (1|T2) + (1|Plate), B[B$Block=="1",])
  mF2 <- lmer(y ~ (1|T2) + (1|Plate), B[B$Block=="2",])
  
  ### Calculate heritability with "standard method"
  v <- as.data.frame(VarCorr(mF))
  h <- v$vcov[1] / ((v$vcov[1]) + (v$vcov[2]/3) + (v$vcov[3]/3))
  
  ### Merge H2 into one data frame
  h <- data.frame(trait=trait[i],
                  year="2020", 
                  H2=h)
  H2 <- rbind(H2, h)

  ### Obtain blups plus coefficinet and merge to one data frame
  b <- data.frame(Taxa=rownames(coef(mF)$T2),
                  a=coef(mF)$T2[,1])
  colnames(b)[2] <- trait[i]
  blup2020 <- plyr::join(blup2020, b, by="Taxa")
  
  b1 <- data.frame(Taxa=rownames(coef(mF1)$T2),
                  b1=coef(mF1)$T2[,1], 
                  trait=trait[i])
  
  b2 <- data.frame(Taxa=rownames(coef(mF2)$T2),
                   b2=coef(mF2)$T2[,1], 
                   trait=trait[i])

  bb <- merge(b1, b2, by=c("Taxa", "trait"))
  bb$year <- "2020"
  
  bc <- data.frame(cor=cor(bb$b1, bb$b2), 
                          trait=trait[i], 
                          year="2020")
  block.cor <- rbind.data.frame(block.cor, bc)
  
  
  block.blup <- rbind(block.blup, bb)
  
  b[,2] <- predict(bn, b[,2], inverse=T)
  blup.TrueScale2020 <- plyr::join(blup.TrueScale2020, b)
  ###########################
  ### Model for 2021 data ###
  ###########################
  
  ## Subset 
  B <- A[A$year=="2021", c("T2", "year", "Plate", "Block", trait[i])]
  
  ### Find and remove outlier
  B$out <- rstatix::is_outlier(B[,5])
  B <- B[B$out==FALSE,]
  
  ### Transform and normalize data
  bn <- bestNormalize(B[,5], standardize = T)
  B$y <- bn$x.t
  
  ### Fit model with lme4
  mF <- lmer(y ~ (1|T2) + (1|Plate), B)
  mF1 <- lmer(y ~ (1|T2) + (1|Plate), B[B$Block=="1",])
  mF2 <- lmer(y ~ (1|T2) + (1|Plate), B[B$Block=="2",])
  
  ### Calculate heritability with "standard method"
  v <- as.data.frame(VarCorr(mF))
  h <- v$vcov[1] / ((v$vcov[1]) + (v$vcov[2]/3) + (v$vcov[3]/3))
  
  ### Merge H2 into one data frame
  h <- data.frame(trait=trait[i],
                  year="2021", 
                  H2=h)
  H2 <- rbind(H2, h)
  
  ### Obtain blups plus coefficinet and merge to one data frame
  b <- data.frame(Taxa=rownames(coef(mF)$T2),
                  a=coef(mF)$T2[,1])
  colnames(b)[2] <- trait[i]
  blup2021 <- plyr::join(blup2021, b, by="Taxa")
  
  b1 <- data.frame(Taxa=rownames(coef(mF1)$T2),
                   b1=coef(mF1)$T2[,1], 
                   trait=trait[i])
  
  b2 <- data.frame(Taxa=rownames(coef(mF2)$T2),
                   b2=coef(mF2)$T2[,1], 
                   trait=trait[i])
  
  bb <- merge(b1, b2, by=c("Taxa", "trait"))
  bb$year <- "2021"
  
  bc <- data.frame(cor=cor(bb$b1, bb$b2), 
                   trait=trait[i], 
                   year="2021")
  block.cor <- rbind.data.frame(block.cor, bc)
  
  
  block.blup <- rbind(block.blup, bb)
  
  
  b[,2] <- predict(bn, b[,2], inverse=T)
  blup.TrueScale2021 <- plyr::join(blup.TrueScale2021, b)
  
rm(b, bn, mF, B, v, h, An, b1, b2, bb, bc, mF1, mF2, mS)
print(paste(trait[i], "Done!!!"))
}

#cy <- block.blup %>%
#  group_by(trait) %>%
#  summarise(cor(b1, b2))
#cy <- data.frame(cor=cy$`cor(b1, b2)`, trait=cy$trait, year="All")
#a  <- rbind(block.cor, cy)

widiv752 <- read.table("data/meta/WiDiv752.geno.ind")
colnames(widiv752) <- "Taxa"
widiv752$widiv752 <- widiv752$Taxa

widiv800 <- read.table("data/meta/WiDiv800.geno.ind")
colnames(widiv800) <- "Taxa"
widiv800$widiv800 <- widiv800$Taxa

widiv797 <- read.table("data/meta/WiDiv797")
colnames(widiv797) <- "Taxa"
widiv797$widiv797 <- widiv797$Taxa

blup.TrueScale2020$widiv797 <- blup.TrueScale2020$Taxa
blup.TrueScale2021$widiv797 <- blup.TrueScale2021$Taxa


widiv797$widiv797 <- toupper(widiv797$widiv797)
widiv797$widiv797 <- gsub(" ", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("_", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("-", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("[.]", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("282SET", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("DK", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("[(]", "", widiv797$widiv797)
widiv797$widiv797 <- gsub("[)]", "", widiv797$widiv797)

widiv800$widiv800 <- toupper(widiv800$widiv800)
widiv800$widiv800 <- gsub(" ", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("_", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("-", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("[.]", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("282SET", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("DK", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("[(]", "", widiv800$widiv800)
widiv800$widiv800 <- gsub("[)]", "", widiv800$widiv800)

widiv752$widiv752 <- toupper(widiv752$widiv752)
widiv752$widiv752 <- gsub(" ", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("_", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("-", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("[.]", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("282SET", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("DK", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("[(]", "", widiv752$widiv752)
widiv752$widiv752 <- gsub("[)]", "", widiv752$widiv752)

colnames(widiv800)[2] <- "widiv797"
colnames(widiv752)[2] <- "widiv797"


blup.TrueScale2020$widiv797 <- gsub("DK", "", blup.TrueScale2020$widiv797)
blup.TrueScale2020$widiv797 <- gsub("GOODMANBUCKLER", "", blup.TrueScale2020$widiv797)

blup.TrueScale2021$widiv797 <- gsub("DK", "", blup.TrueScale2021$widiv797)


blup2020gwasV5 <- plyr::join(widiv797, blup.TrueScale2020, by="widiv797")
blup2021gwasV5 <- plyr::join(widiv797, blup.TrueScale2021, by="widiv797")

blup2020gwasV5 <- blup2020gwasV5[,-c(1:2)]
blup2021gwasV5 <- blup2021gwasV5[,-c(1:2)]

blup2020gwasV4 <- plyr::join(widiv800, blup.TrueScale2020, by="widiv797")
blup2021gwasV4 <- plyr::join(widiv800, blup.TrueScale2021, by="widiv797")
blup2021gwasV4 <- blup.TrueScale2021[!duplicated(blup2021gwasV4$Taxa),]

blup2020gwasRNA <- plyr::join(widiv752, blup.TrueScale2020, by="widiv797")
blup2021gwasRNA <- plyr::join(widiv752, blup.TrueScale2021, by="widiv797")

write.csv(blup2020gwasV5, "data/work/blups/blup2020gwasV5.csv",row.names = F)
write.csv(blup2021gwasV5, "data/work/blups/blup2021gwasV5.csv",row.names = F)
write.csv(blup2020gwasV4, "data/work/blups/blup2020gwasV4.csv",row.names = F)
write.csv(blup2021gwasV4, "data/work/blups/blup2021gwasV4.csv",row.names = F)
write.csv(blup2020gwasRNA, "data/work/blups/blup2020gwasRNA.csv",row.names = F)
write.csv(blup2021gwasRNA, "data/work/blups/blup2021gwasRNA.csv", row.names = F)

write.csv(H2, "results/data/Broad_Heritability.csv", row.names = F)
write.csv(P, "results/data/Year_Effect.csv", row.names = F)
