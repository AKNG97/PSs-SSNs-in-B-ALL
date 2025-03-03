#This script shows how to calculate the single-sample edge weights 
#of new samples relative to a referene dataset.

#Since the only the RPL7P10-RPS3AP36 interaction showed 
#stable cross-cohorts results with this approach, this script only
#considers the calculation of this edge.

#### Load Norm libraries ####
library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(Hmisc)
library(tidyr)
library(ggbiplot)

####Load expression matrices ####

#Normalized TARGET matrix used for the biomarkers discovery
TARGET = readRDS("NormData/BALL_TARGET_norm.RDS")
#Raw matrix from TARGET
TARGET_RAW = readRDS("RawData/BALL_TARGET_raw.RDS")
#Raw matrix from MP2PRT
MP2_RAW <- readRDS("RawData/BALL_MP2PRT_raw.RDS")
CI.MP2 <- as.data.frame(colData(MP2_RAW))
#### NBM Raw matrix
NBM_raw <- readRDS("RawData/NBM_raw.RDS")

#### Filter data ####
#Load annotation file
annot <- readRDS("annot_normalized_TARGET.RDS")

MP2_RAW <- assay(MP2_RAW)

MP2_RAW <- MP2_RAW[rownames(MP2_RAW) %in% annot$gene_id,]
TARGET_RAW <- TARGET_RAW[rownames(TARGET_RAW) %in% annot$gene_id,]
NBM_raw <- NBM_raw[rownames(NBM_raw) %in% annot$gene_id,]

length(intersect(rownames(TARGET), annot$HGNC_symbol))

rownames(MP2_RAW) <- annot$HGNC_symbol[match(rownames(MP2_RAW), annot$gene_id)]
rownames(TARGET_RAW) <- annot$HGNC_symbol[match(rownames(TARGET_RAW), annot$gene_id)]
rownames(NBM_raw) <- annot$HGNC_symbol[match(rownames(NBM_raw), annot$gene_id)]

MP2_RAW <- MP2_RAW[rownames(MP2_RAW) %in% rownames(TARGET),]
TARGET_RAW <- TARGET_RAW[rownames(TARGET_RAW) %in% rownames(TARGET),]
NBM_raw <- NBM_raw[rownames(NBM_raw) %in% rownames(TARGET),]

#### Norm for each Sample in MP2, save only PS of interest, par, save Norm per MP2 ####
relevantPS <- read.csv("RPL7P10", "RPS3AP36")

TARGET_NBM <- cbind(TARGET_RAW, NBM_raw)

annot <- annot[annot$HGNC_symbol %in% rownames(TARGET_RAW),]

annot <- annot[match(rownames(TARGET_NBM), annot$HGNC_symbol),]

CI_BALL <- data.frame(Sample = colnames(TARGET_NBM), 
                      Group = c(rep("BALL", ,132), 
                                rep("Normal Bone Marrow", 70)))

#Get GlobalCI

CI.BALL_TCGA <- read.delim("OS_files/CI_BALL.tsv")
CI.BALL_MP2 <- read.delim("OS_files/CI_BALL_MP2.tsv")
CI.BALL_MP2$sample_submitter_id <- gsub("-", ".", CI.BALL_MP2$sample_submitter_id)
CI.BALL_TCGA$patient <- gsub("-", ".", CI.BALL_TCGA$patient)

CI.BALL_MP2 <- CI.BALL_MP2 %>% dplyr::select(Sample, VitalStatus, OS_time)

CI.BALL_TCGA <- CI.BALL_TCGA %>% dplyr::select(Sample, VitalStatus, OS_time)

CI.BALL_TCGA$VitalStatus[CI.BALL_TCGA$VitalStatus == "Alive"] <- 0
CI.BALL_TCGA$VitalStatus[CI.BALL_TCGA$VitalStatus == "Dead"] <- 1
CI.BALL_TCGA$VitalStatus <- as.numeric(CI.BALL_TCGA$VitalStatus)

GlobalCI <- rbind(CI.BALL_MP2, CI.BALL_TCGA) 

GlobalCI$VitalStatus[GlobalCI$DLF > 1825] <- 0
GlobalCI$DLF[GlobalCI$DLF > 1825] <- 1825

colnames(TARGET_NBM) <- gsub("-", ".", colnames(TARGET_NBM))
colnames(MP2_RAW) <- gsub("-", ".", colnames(MP2_RAW))

norm <- function(x, y, z) {
  
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$Group))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  rnas2 <- rnas2[rownames(rnas2) %in% relevantPS$PS,]
  return(rnas2)
}

write.csv(GlobalCI, "GlobalCI.csv", quote = F, row.names = F)

# >>>>> LOOP STARTS HERE  <<<<<<<<
#Dirs
ExprMs <- "ExprM_MP2/"
dir.create(ExprMs)

CI_BALL <- data.frame(Sample = colnames(TARGET_NBM), 
                      Group = c(rep("BALL", 132), rep("Normal Bone Marrow", 70)))

IterateNorm <- function(i) {
  Sample_i <- colnames(MP2_RAW)[i]
  
  TARGET_NBM_MP2i <- cbind(TARGET_NBM, MP2_RAW[,i])
  colnames(TARGET_NBM_MP2i)[203] <- Sample_i
  CI_BALL_i <- rbind(CI_BALL, c(Sample_i, "BALL"))
  rownames(CI_BALL_i) <- CI_BALL_i$Sample
  
  TARGET_NBM_MP2i_norm <- norm(TARGET_NBM_MP2i, annot, CI_BALL_i)
  TARGET_NBM_MP2i_norm <- TARGET_NBM_MP2i_norm[,c(1:132,203)]
  
  saveRDS(TARGET_NBM_MP2i_norm, paste0(ExprMs, "Expr_MP2_Sample_", i, ".RDS"))
  
}
library(future.apply)
future::plan(multisession, workers = 40)
results <- future_lapply(1:ncol(MP2_RAW), FUN = IterateNorm, future.seed = TRUE)

plan(sequential)
#### Get SSN ######

library(dplyr)
library(tidyr)
ExprMs <- "ExprM_MP2/"
GCNs <- "GCNs_MP2/"
SSNs <- "SSNs_MP2/"
dir.create(GCNs)
dir.create(SSNs)

relevantEdges <- c("RPL7P10-RPS3AP36")

#Get GCN
GCN <- data.frame(Edge = relevantEdges$Edge, 
                  Rho = NA) %>% dplyr::group_by(Edge) %>%
  dplyr::mutate(Source = stringr::str_split(Edge,"-")[[1]][1], 
                Target = stringr::str_split(Edge,"-")[[1]][2])

IterateGCN_Lioness <- function(i) {
  
  TARGET_NBM_MP2i_norm <- readRDS(paste0(ExprMs, "Expr_MP2_Sample_", i, ".RDS"))
  
  GCN_i <- GCN %>% 
    dplyr::mutate(Rho = cor(TARGET_NBM_MP2i_norm[rownames(TARGET_NBM_MP2i_norm) == Source,],
                            TARGET_NBM_MP2i_norm[rownames(TARGET_NBM_MP2i_norm) == Target,],
                            method = "spearman"))
  
  write.csv(GCN_i, paste0(GCNs, "GCN_MP2_Sample_", i, ".csv"), quote = F, row.names = F)
  n_samples <- 133
  
  #Get Lioness
  for(j in 1:nrow(GCN_i)){
    
    Source <- GCN_i$Source[j]
    Target <- GCN_i$Target[j]
    alpha <- GCN_i$Rho[j]
    SSN_j <- numeric()
    
    for(x in 1:ncol(TARGET_NBM_MP2i_norm)){
      
      j_cor <- stats::cor(TARGET_NBM_MP2i_norm[rownames(TARGET_NBM_MP2i_norm) == Source,-x],
                          TARGET_NBM_MP2i_norm[rownames(TARGET_NBM_MP2i_norm) == Target,-x], 
                          method="spearman")
      SSN_cor <- round(n_samples*(alpha - j_cor) + j_cor, 4)
      #names(SSN_cor) <- paste(Source, Target, sep = "_")
      
      if(length(SSN_j) == 0){
        SSN_j <- SSN_cor
      } else {
        SSN_j <- c(SSN_j, SSN_cor)
      }
    }
    
    if(j == 1){
      SSN_i <- SSN_j
    } else {
      SSN_i <- rbind(SSN_i, SSN_j)
    }
  }
  
  SSN_i <- as.data.frame(SSN_i)
  rownames(SSN_i) <- GCN_i$Edge
  colnames(SSN_i) <- colnames(TARGET_NBM_MP2i_norm)
  
  write.csv(SSN_i, paste0(SSNs,"SSN_MP2_Sample_", i, ".csv"), quote = F, row.names = T)
  
}

library(future.apply)
future::plan(multisession, workers = 20)
results <- future_lapply(1:1284, FUN = IterateGCN_Lioness, future.seed = TRUE)

future::plan(sequential)

#### Get full matrix of MP2PRT data ####
library(dplyr)
library(tidyr)
library(M3C)
library(gtsummary)
library(ggsurvfit)
library(survival)
library(caret)

SSNs <- "SSNs_MP2/"

ClassifyMP2_Samples <- function(i){
  
  SSN_i <- read.csv(paste0(SSNs,"SSN_MP2_Sample_", i, ".csv")) %>% as.data.frame()

  Sample_i <- colnames(SSN_i)[133]
  
  #Exclude de test sample (MP2PRT) from the mean and sd calculation
  original_data_relevant <- SSN_i[, 1:132]
  row_means <- rowMeans(original_data_relevant)
  row_sds <- apply(original_data_relevant, 1, sd)
  
  #Scale all samples, including the MP2PRT sample
  SSN_i_scaled <- (SSN_i - row_means) / row_sds
  
  #Identify the risk category based on the edge weight
  RiskCategoryi <- ifelse(SSN_i_scaled[1,133] > -0.741, "High risk", "Low risk")
  
  MP2_strat <- data.frame(Sample = Sample_i,
                          RPL7P10_RPS3AP36 = SSN_i_scaled[1,133],
                          Category = RiskCategoryi)
  
  return(MP2_strat)
  
}

library(future.apply)
future::plan(multisession, workers = 3)
results <- future_lapply(1:1284, FUN = ClassifyMP2_Samples, future.seed = TRUE)
Relative_Strata <- do.call(rbind, results)
future::plan(sequential)

write.csv(Relative_Strata, "MP2PRT_Stratification_Relative_to_TARGET.csv", quote = F, row.names = F)



