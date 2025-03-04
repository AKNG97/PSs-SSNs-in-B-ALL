library(ComplexHeatmap)
library(dplyr)
library(ggsurvfit)
library(survival)
library(circlize)
library(NOISeq)
library(edgeR)
library(M3C)
library(data.table)
library(ggplot2)

##### TARGET #####

Lioness_BALL_TARGET = read.csv("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/Lioness_TARGET_BALL_FULL_final.csv", row.names = 1)

Cancer_GCN <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/Bootsrtrapping/BALL_GCN_Sp_sequential.RDS")

annot <- read.delim("gene_annotation_TCGA.tsv")

annot_PS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]

Cancer_GCN <- Cancer_GCN %>% filter(Source %in% annotPS$HGNC_symbol &
                                      Target %in% annotPS$HGNC_symbol)

#Get network without Ps-Ps edges
Lioness_BALL_TARGET1 = Lioness_BALL_TARGET[!(rownames(Lioness_BALL_TARGET) %in%
                                              Cancer_GCN$ID),]

#Get network with Ps-Ps edges
Lioness_BALL_TARGET2 = Lioness_BALL_TARGET[rownames(Lioness_BALL_TARGET) %in%
                                               Cancer_GCN$ID,]

#Filter network1
filtered_results_TARGET1 <- featurefilter(Lioness_BALL_TARGET1, percentile=1508*100/39082,
                                         method='var', topN=5)
#Filter network2
filtered_results_TARGET2 <- featurefilter(Lioness_BALL_TARGET2, percentile=25,
                                         method='var', topN=5)

Lioness_BALL_TARGET_scaled1 <- as.data.frame(t(scale(t(filtered_results_TARGET1$filtered_data))))
Lioness_BALL_TARGET_scaled2 <- as.data.frame(t(scale(t(filtered_results_TARGET1$filtered_data))))

r.hc_TARGET <- M3C(Lioness_BALL_TARGET_scaled2, method=1, clusteralg = "hc", cores = 30, maxK = 10, seed = 123)

ggsave("pval_M3C_BALL_Target_PS.png", plot = r.hc_TARGET$plots[[3]], dpi = 300)
ggsave("RCSI_M3C_BALL_Target_PS.png", plot = r.hc_TARGET$plots[[4]], dpi = 300)

#TARGET
assigments.TARGET <- r.hc_TARGET$realdataresults[[2]]$ordered_annotation
assigments.TARGET <- assigments.TARGET %>% mutate(Sample = rownames(assigments.TARGET))
saveRDS(assigments.TARGET, "assigmentsTARGET_PS_25pMostVariable.RDS")

Lioness_BALL_TARGET_scaled <- as.data.frame(t(scale(t(Lioness_BALL_TARGET))))
assigments.TARGET <- readRDS("assigmentsTARGET_PS_25pMostVariable.RDS")
colnames(Lioness_BALL_TARGET_scaled) <- gsub("-", ".", colnames(Lioness_BALL_TARGET_scaled))
Lioness_BALL_TARGET_pivot <- Lioness_BALL_TARGET_scaled %>% mutate(Edge = rownames(Lioness_BALL_TARGET_scaled)) %>% 
  pivot_longer(!Edge, names_to = "Sample", values_to = "Z_score_Lioness") %>% 
  inner_join(assigments.TARGET, by = c("Sample" = "Sample"))

diff_TARGET <- Lioness_BALL_TARGET_pivot %>% 
  group_by(Edge) %>% 
  dplyr::summarize(Median1 = median(Z_score_Lioness[consensuscluster == 1]),
                   Median2 = median(Z_score_Lioness[consensuscluster == 2]),
                   Median3 = median(Z_score_Lioness[consensuscluster == 3]),
                   Diff_1vs2 = median(Z_score_Lioness[consensuscluster == 1]) - median(Z_score_Lioness[consensuscluster == 2]),
                   Diff_1vs3 = median(Z_score_Lioness[consensuscluster == 1]) - median(Z_score_Lioness[consensuscluster == 3]),
                   Diff_2vs3 = median(Z_score_Lioness[consensuscluster == 2]) - median(Z_score_Lioness[consensuscluster == 3]),
                  pval.Cls1_vs2 = wilcox.exact(Z_score_Lioness[consensuscluster == 1],
                                                Z_score_Lioness[consensuscluster == 2])$p.value,
                   pval.Cls1_vs3 = wilcox.exact(Z_score_Lioness[consensuscluster == 1],
                                                Z_score_Lioness[consensuscluster == 3])$p.value,
                   pval.Cls2_vs3 = wilcox.exact(Z_score_Lioness[consensuscluster == 2],
                                                Z_score_Lioness[consensuscluster == 3])$p.value
  ) %>%
 mutate(adj.pval_1vs2 =  p.adjust(pval.Cls1_vs2, method = "fdr"),
         adj.pval_1vs3 =  p.adjust(pval.Cls1_vs3, method = "fdr"),
         adj.pval_2vs3 =  p.adjust(pval.Cls2_vs3, method = "fdr"))

diff_TARGET %>% filter(adj.pval_2vs3 < 0.05) %>% group_by(Edge) %>% 
  summarise(diff = Median2 - Median3) %>%
  dplyr::count(pos = diff > 0.5,
               neg = diff < -0.5)

diff_TARGET <- diff_TARGET %>% inner_join(Cancer_GCN[,c(1,2,4)], by = c("Edge" = "ID"))

write.csv(diff_TARGET, "DiffCoExpr_TARGET_25p.csv", quote = F, row.names = F)


