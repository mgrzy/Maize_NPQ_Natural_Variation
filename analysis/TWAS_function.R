TWAS <- function(phenotype, expressionMatrix, covariate)
{
 
  commonName <- colnames(phenotype)[1]

  X <- merge(phenotype, expressionMatrix, by=commonName)
  X <- as.data.frame(X)
  
  toRemove <- which(as.vector(is.na(X[,2])))

  X <- X[-c(toRemove),]  
  covariate <- covariate[-c(toRemove),]
  
  ph <- as.numeric(X[,2])
  
  genExp <- X[,-c(1,2)]
  
  genToKeep <- apply(genExp, 2, FUN = function(x) {
    b <- table(x > 0)
    z <- as.numeric(b["TRUE"]) > (nrow(genExp)/2)
    return(z)}
  )
  
  genToKeep <- which(genToKeep==TRUE)
  genExp <- genExp[,genToKeep]
  genExp <- as.data.frame(genExp)
  
  covariate <- as.matrix(covariate[,-c(1)])
  
  res <- apply(genExp, 2, function(x) 
    summary(lm(ph ~ covariate + x))$coefficients[7,4])
  
  res <- data.frame(gene=names(res), p_val=res)
  
  res <- merge(res, gene)
  
  res <- res[res$chromosome %in% c(1:10),]
  res$chromosome <- as.numeric(res$chromosome)
  res$ist <- p.adjust(res$p_val) < 0.05
  
  gQQ <- ggplot(res, aes(
    y = -log10(sort(p_val, decreasing = F)), 
    x = -log10(ppoints(length(p_val))))) + 
    geom_point() +
    geom_abline(slope = 1) +
    xlim(0, max(-log10(res$p_val) + 0.5)) +
    ylim(0, max(-log10(res$p_val) + 0.5)) + 
    ylab("Observed -log10(p)") + 
    xlab("Expected -log10(p)")
  
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
  
  gMan <- ggplot() + 
    geom_hline(yintercept = -log10(0.05/nrow(res)) , linetype=2) +
    geom_point(data=res, aes(BPcum, -log10(p_val), colour=as.character(chromosome))) + 
    #scale_color_d3() + 
    scale_x_continuous(label = axis.set$chromosome, breaks = axis.set$center) + 
    theme(legend.position = "none") + 
    ylab("-log10(p-value)") +  
    xlab("Chromosome")
    ggtitle(colnames(dat)[1])
  
}


library(sommer)

K <- as.matrix(bigmemory::attach.big.matrix("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.kin.desc"))
ind <- read.table("../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/rMVP/WiDiv752.geno.ind")

rownames(K) <- ind$V1
colnames(K) <- ind$V1

m <- mmer(NPQmaxme~1 + covariate,
     random=~vsr(Taxa,Gu=K),
     rcov=~units, nIters=3,
     data=z, verbose = FALSE)

m2 <- mmer(NPQmaxme ~  covariate + Zm00001d042697,
          random=~vsr(Taxa,Gu=K),
          rcov=~units, nIters=3,
          data=z, verbose = FALSE)


res2 <- apply(genExp[,1:5], 2, function(x){
  m2 <- mmer(NPQmaxme ~  covariate + x,
             random=~vsr(Taxa,Gu=K),
             rcov=~units, nIters=3,
             data=z, verbose = FALSE)})
