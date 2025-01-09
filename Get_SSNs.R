
Cancer_ExprM <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/NormData/MP2PRT_BM_norm_FULL.RDS")
Cancer_GCN <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/BALL_MP2_GCN_SignifEdges_2ndRun.RDS")
Cancer_ExprM <- Cancer_ExprM[,1:1284]
#Cancer_GCN <- data.table::fread("BALL_TARGET_and_NBM_GCN_SignifEdges.tsv")
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
      
      SSN <- data.frame(Edge_ID = paste(pmin(Source, Target), pmax(Source, Target), sep ="_"),
                        Sample = round(SSN_cor, digits = 4))
      colnames(SSN) = c("Edge_ID", colnames(Cancer_ExprM)[sample_n])
    } else {
      SSN[nrow(SSN) + 1,] <- c(paste(pmin(Source, Target), pmax(Source, Target), sep ="_"),
                               round(SSN_cor, digits = 4))
    }
  }
  
  write.csv(SSN, paste0("SSNs_check_MP2/", sample, "SSN.csv"), quote = FALSE)
  return(SSN)
  
}

Lioness_ParSamples2 <- function(sample){
  
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
      
      # SSN <- data.frame(Edge_ID = paste(pmin(Source, Target), pmax(Source, Target), sep ="_"),
      #                   Sample = round(SSN_cor, digits = 4))
      # colnames(SSN) = c("Edge_ID", colnames(Cancer_ExprM)[sample_n])
    } else {
      # SSN[nrow(SSN) + 1,] <- c(paste(pmin(Source, Target), pmax(Source, Target), sep ="_"),
      #                          round(SSN_cor, digits = 4))
      
      SSN <- c(SSN, SSN_cor)
      
    }
  }
  
  SSN <- tibble(Edge_ID = Cancer_GCN$ID,
                Sample = SSN)
  
  colnames(SSN) = c("Edge_ID", colnames(Cancer_ExprM)[sample_n])
  
  write.csv(SSN[1:5, ], paste0("SSNs_check_MP2/", sample, "SSN.csv"), quote = FALSE)
  return(SSN)
  
}

sample.list <- colnames(Cancer_ExprM)
future::plan(multisession, workers = 45)
results <- future_lapply(sample.list, FUN = Lioness_ParSamples2, future.seed = TRUE)
SSN <- Reduce(function(x, y) inner_join(x, y, by = "Edge_ID"), results)
dim(SSN)
write.csv(SSN, "LionessPS_BALL_MP2_Full_CheckError.csv", quote = FALSE, row.names = F)
