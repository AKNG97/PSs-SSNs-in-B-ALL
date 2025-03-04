#This script swhows how to calculate the
#LIONESS networks in an expression matrix

library(dplyr)
library(tidyr)
library(Hmisc)
library(future)
library(future.apply)

Cancer_ExprM.TARGET <- readRDS("NormData/BALL_TARGET_norm.RDS")
Cancer_GCN.TARGET <- readRDS("GCN_BALL_TARGET.RDS")
Cancer_ExprM.TARGET <- Cancer_ExprM.TARGET[,1:132]
Cancer_GCN.TARGET$Sp_corr <- as.numeric(Cancer_GCN.TARGET$Sp_corr)

Cancer_ExprM.MP2 <- readRDS("NormData/BALL_MP2PRT_norm.RDS")
Cancer_GCN.MP2 <- readRDS("GCN_BALL_MP2PRT.RDS")
Cancer_ExprM.MP2 <- Cancer_ExprM.MP2[,1:1284]
Cancer_GCN.MP2$Sp_corr <- as.numeric(Cancer_GCN.MP2$Sp_corr)

Lioness_ParSamples <- function(sample){
  
  sample_n <- which(sample.list == sample)
  
  n_samples <- ncol(Cancer_ExprM)
  ExprM2 <- Cancer_ExprM[,-sample_n]
  
  for(j in 1:nrow(Cancer_GCN)){
    
    Source <- Cancer_GCN[j,"Source"] 
    Target <- Cancer_GCN[j,"Target"]
    
    j_cor <- stats::cor(ExprM2[rownames(ExprM2) == Source,],
                        ExprM2[rownames(ExprM2) == Target,], method="spearman")
    
    SSN_cor <- n_samples*(Cancer_GCN[j,"Sp_corr"] - j_cor) + j_cor
    
    if(j == 1){
      
      SSN <- SSN_cor

    } else {

      SSN <- c(SSN, SSN_cor)
      
    }
  }
  
  SSN <- tibble(Edge_ID = Cancer_GCN$ID,
                Sample = SSN)
  
  colnames(SSN) = c("Edge_ID", colnames(Cancer_ExprM)[sample_n])

  return(SSN)
  
}

sample.list <- colnames(Cancer_ExprM.TARGET)
future::plan(multisession, workers = 45)
results <- future_lapply(sample.list, FUN = Lioness_ParSamples, future.seed = TRUE)
SSN <- Reduce(function(x, y) inner_join(x, y, by = "Edge_ID"), results)
dim(SSN)
write.csv(SSN, "Lioness_BALL_TARGET.csv", quote = FALSE, row.names = F)

sample.list <- colnames(Cancer_ExprM.MP2)
future::plan(multisession, workers = 45)
results <- future_lapply(sample.list, FUN = Lioness_ParSamples, future.seed = TRUE)
SSN <- Reduce(function(x, y) inner_join(x, y, by = "Edge_ID"), results)
dim(SSN)
write.csv(SSN, "Lioness_BALL_MP2.csv", quote = FALSE, row.names = F)

              
