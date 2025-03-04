#This script swhows how to calculate the
#LIONESS networks in an expression matrix

Cancer_ExprM <- readRDS("NormData/BALL_TARGET_norm.RDS")
Cancer_GCN <- readRDS("BALL_TARGET.RDS")
Cancer_ExprM <- Cancer_ExprM[,1:132]
Cancer_GCN$Sp_corr <- as.numeric(Cancer_GCN$Sp_corr)

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

sample.list <- colnames(Cancer_ExprM)
future::plan(multisession, workers = 45)
results <- future_lapply(sample.list, FUN = Lioness_ParSamples, future.seed = TRUE)
SSN <- Reduce(function(x, y) inner_join(x, y, by = "Edge_ID"), results)
dim(SSN)
write.csv(SSN, "Lioness_BALL_TARGET.csv", quote = FALSE, row.names = F)
