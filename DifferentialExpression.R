
library(dplyr)
library(rstatix)
library(future.apply)
library(NOISeq)
library(edgeR)
library(tidyr)
library(exactRankTests)
library(DESeq2)
library(SummarizedExperiment)

BALL_TARGET_RAW = readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/RawData/BALL_primaryBM_raw.RDS")

BALL_MP2_RAW = assay(readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/MP2PRT_raw_BM_FULL.RDS"))

AML_NBM <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/AML_NBM_download_061023.RDS")
NBM <- AML_NBM[,AML_NBM$sample_type == "Bone Marrow Normal"]
table(as.factor(NBM$primary_diagnosis))
# Acute myeloid leukemia, NOS 
# 145 
#Hence, there are "normal" samples from AML patients
which(duplicated(NBM$patient))
#integer(0)
table(as.factor(NBM$vital_status))
NBM <- NBM[, is.na(NBM$primary_diagnosis)]
dim(NBM)
#[1] 60660    70
NBM_raw <- assay(NBM)
CI_NBM <- as.data.frame(colData(NBM)) %>% mutate(primary_diagnosis = "Normal_Bone_Marrow")
colnames(NBM_raw) <- CI_NBM$patient

annot <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/annot_TCGA.RDS")

BALL_TARGET <- cbind(BALL_TARGET_RAW, NBM_raw)
BALL_TARGET <- BALL_TARGET[rownames(BALL_TARGET) %in% annot$gene_id,]
colnames(BALL_TARGET) <- gsub("-", ".", colnames(BALL_TARGET))
rownames(BALL_TARGET) <- annot$HGNC_symbol[match(rownames(BALL_TARGET), annot$gene_id)]

BALL_MP2 <- cbind(BALL_MP2_RAW, NBM_raw)
BALL_MP2 <- BALL_MP2[rownames(BALL_MP2) %in% annot$gene_id,]
colnames(BALL_MP2) <- gsub("-", ".", colnames(BALL_MP2))
rownames(BALL_MP2) <- annot$HGNC_symbol[match(rownames(BALL_MP2), annot$gene_id)]

keep <- filterByExpr(BALL_TARGET, group = c(rep("BALL", 132), rep("NBM", 70)))
BALL_TARGET_raw <- BALL_TARGET[keep, ]
BALL_TARGET_tmm <- tmm(BALL_TARGET_raw, long = 1000, lc = 0, k = 0)
dim(BALL_TARGET_tmm)

keep <- filterByExpr(BALL_MP2, group = c(rep("BALL", 1284), rep("NBM", 70)))
BALL_MP2_raw <- BALL_MP2[keep, ]
BALL_MP2_tmm <- tmm(BALL_MP2_raw, long = 1000, lc = 0, k = 0)
dim(BALL_MP2_tmm)

saveRDS(BALL_MP2, "BALL_MP2_tmm_expr.RDS")
saveRDS(BALL_TARGET_tmm, "BALL_target_tmm_expr.RDS")
