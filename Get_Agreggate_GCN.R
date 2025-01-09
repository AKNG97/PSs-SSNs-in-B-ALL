library(dplyr)
library(tidyr)
library(Hmisc)
library(future)
library(future.apply)

dir.create("GCNs/")

#Load expression matrices

BALL_TARGET_norm <- readRDS("NormData/BALL_TARGET_norm.RDS")
BALL_MP2PRT_norm <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/NormData/MP2PRT_BM_norm_FULL.RDS")
BALL_MP2PRT_norm <- BALL_MP2PRT_norm[,1:1284]

#Paralelo
par_Sp_Corr <- function(expr_m, x){
  
  gene1 <- rownames(expr_m)[x]
  
  for(y in 1:nrow(expr_m)){
    
    if(y != x){
      
      gene2 <- rownames(expr_m)[y]
      x_expr = expr_m[x, ]
      y_expr = expr_m[y, ]
      
      net <- rcorr(x_expr, y_expr, type = "spearman")
      
      if(net$P[1,2] < 1e-08){
        
        if(!exists("coexpr")){
          coexpr <- data.frame(Source = gene1,
                               Target = gene2,
                               Sp_corr = net$r[1,2],
                               p_value = net$P[1,2],
                               abs = abs(net$r[1,2]),
                               ID = paste(pmin(gene1,gene2),pmax(gene1,gene2),sep="_"))
        } else {
          coexpr[nrow(coexpr) + 1,] <- c(gene1, gene2, net$r[1,2], net$P[1,2], abs(net$r[1,2]),
                                         paste(pmin(gene1,gene2),pmax(gene1,gene2),sep="_"))
        }
      }
    }
  }
  if(exists("coexpr")){
    saveRDS(coexpr, paste0("GCNs_par2/GCN_", gene1, ".RDS"))
    #write.csv(coexpr, paste0("GCNs_par/GCN_", gene1, ".csv"), quote = F, row.names = F)
    return(coexpr)
  }
}
dir.create("GCNs_par2")
future::plan(multisession, workers = 60)
results <- future_lapply(expr_m = BALL_MP2PRT_norm, 1:nrow(BALL_MP2PRT_norm), 
                         FUN = par_Sp_Corr, future.seed = TRUE)
GCN <- do.call(rbind, results)
GCN <- GCN[!duplicated(GCN$ID),]
GCN <- GCN %>% as.data.frame() %>% arrange(desc(abs))
#write.csv(GCN[1:1e+05,], "BALL_GCN_Sp_100k.csv", quote = FALSE, row.names = F)
#write.csv(GCN[1:1e+07,], "BALL_GCN_Sp_10M.csv", quote = FALSE, row.names = F)
#write.csv(GCN[1:1e+04,], "BALL_GCN_Sp_10k.csv", quote = FALSE, row.names = F)
saveRDS(GCN, "BALL_MP2_GCN_SignifEdges_2ndRun.RDS")
#saveRDS(GCN %>% filter(p_value < 1e-08), "GCNs/BALL_TARGET_and_NBM_GCN_SignifEdges.RDS")
