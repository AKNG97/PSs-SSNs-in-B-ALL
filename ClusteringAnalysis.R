#This script shows how to perform the clustering analysis of samples
#using the LIONESS networks

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

Lioness_BALL_TARGET = read.csv("Lioness_BALL_TARGET.csv", row.names = 1)

Cancer_GCN <- readRDS("GCN_BALL_TARGET.RDS")

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
Lioness_BALL_TARGET_scaled2 <- as.data.frame(t(scale(t(filtered_results_TARGET2$filtered_data))))

#Run Clustering algorithm with each filterd and scaled network, here we run it with the Ps-Ps network
r.hc_TARGET <- M3C(Lioness_BALL_TARGET_scaled2, method=1, clusteralg = "hc", cores = 30, maxK = 10, seed = 123)

ggsave("pval_M3C_BALL_Target_PS.png", plot = r.hc_TARGET$plots[[3]], dpi = 300)
ggsave("RCSI_M3C_BALL_Target_PS.png", plot = r.hc_TARGET$plots[[4]], dpi = 300)

#TARGET
assigments.TARGET <- r.hc_TARGET$realdataresults[[3]]$ordered_annotation
assigments.TARGET <- assigments.TARGET %>% mutate(Sample = rownames(assigments.TARGET))
saveRDS(assigments.TARGET, "assigmentsTARGET_PS_25pMostVariable.RDS")

CI.BALL_TARGET <- read.delim("OS_files/CI_BALL.tsv")
#KM
CI.BALL_TARGET$Sample <- gsub("-", ".", CI.BALL_TARGET$Sample)
CI.BALL_TARGET <- CI.BALL_TARGET %>% inner_join(assigments.TARGET, by = c("Sample" = "Sample"))

surv.plot.BALL.TARGET <- survfit2(Surv(OS_time, VitalStatus) ~ consensuscluster, data = CI.BALL_TARGET) %>%
  ggsurvfit(size = 1.5) +
  labs(
    x = "Days to last follow up",
    y = "Overall survival probability"
  ) + add_pvalue() + theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=14)) +
  scale_color_manual(values = c("#F8766D", "#619CFF", "#00BA38")) +
  add_risktable(
    risktable_height = 0.2,
    size = 4, # increase font size of risk table statistics
    theme =   # increase font size of risk table title and y-axis label
      list(
        theme_risktable_default(axis.text.y.size = 12, 
                                plot.title.size = 12),
        theme(plot.title = element_text(face = "bold"))
      )
  )

ggsave("KM_TARGET_PS.png", surv.plot.BALL.TARGET, height = 10, width = 7, dpi = 300)

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

##### MP2PRT #####
Lioness_BALL_MP2 = data.table::fread("/Lioness_BALL_MP2.csv")
Lioness_BALL_MP2 <- Lioness_BALL_MP2 %>% as.data.frame
rownames(Lioness_BALL_MP2) <- Lioness_BALL_MP2$Edge_ID
Lioness_BALL_MP2 <- Lioness_BALL_MP2[,-1]
Cancer_GCN <- readRDS("GCN_BALL_MP2PRT.RDS")

annot <- read.delim("gene_annotation_TCGA.tsv")

annot_PS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]

#Get network without Ps-Ps edges
Lioness_BALL_MP21 = Lioness_BALL_TARGET[!(rownames(Lioness_BALL_MP2) %in%
                                              Cancer_GCN$ID),]

#Get network with Ps-Ps edges
Lioness_BALL_MP22 = Lioness_BALL_TARGET[rownames(Lioness_BALL_MP2) %in%
                                               Cancer_GCN$ID,]

#Filter network1
filtered_results_MP21 <- featurefilter(Lioness_BALL_MP21, percentile=4811*100/524894,
                                         method='var', topN=5)
#Filter network2
filtered_results_MP22 <- featurefilter(Lioness_BALL_MP22, percentile=25,
                                         method='var', topN=5)

Lioness_BALL_MP2_scaled1 <- as.data.frame(t(scale(t(filtered_results_MP21$filtered_data))))
Lioness_BALL_MP2_scaled2 <- as.data.frame(t(scale(t(filtered_results_MP22$filtered_data))))

#The Clustering algorithm should be performed with each filterd and scaled network, here we only run it with the Ps-Ps network
r.hc_MP2 <- M3C(Lioness_BALL_MP2_scaled2, method=1, clusteralg = "hc", cores = 30, maxK = 10, seed = 123)

assigments.MP2 <- r.hc_MP2$realdataresults[[2]]$ordered_annotation
assigments.MP2 <- assigments.MP2 %>% mutate(Sample = rownames(assigments.MP2))
saveRDS(assigments.MP2, "assigmentsMP2_PS_25pMostVariable.RDS")

CI.BALL_MP2 <- read.delim("OS_files/CI_BALL_MP2.tsv")
#KM
CI.BALL_MP2$Sample <- gsub("-", ".", CI.BALL_MP2$Sample)
CI.BALL_MP2 <- CI.BALL_MP2 %>% inner_join(assigments.MP2, by = c("Sample" = "Sample"))

surv.plot.BALL.MP2 <- survfit2(Surv(OS_time, VitalStatus) ~ consensuscluster, data = CI.BALL_MP2) %>%
  ggsurvfit(size = 1.5) +
  labs(
    x = "Days to last follow up",
    y = "Overall survival probability"
  ) + add_pvalue() + theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=14)) +
  scale_color_manual(values = c("#F8766D", "#619CFF", "#00BA38")) +
  add_risktable(
    risktable_height = 0.2,
    size = 4, # increase font size of risk table statistics
    theme =   # increase font size of risk table title and y-axis label
      list(
        theme_risktable_default(axis.text.y.size = 12, 
                                plot.title.size = 12),
        theme(plot.title = element_text(face = "bold"))
      )
  )

ggsave("KM_MP2_PS.png", surv.plot.BALL.TARGET, height = 10, width = 7, dpi = 300)

Lioness_BALL_MP2_scaled <- as.data.frame(t(scale(t(Lioness_BALL_MP2))))
assigments.MP2 <- readRDS("assigmentsMP2_PS_25pMostVariable.RDS")
colnames(Lioness_BALL_MP2_scaled) <- gsub("-", ".", colnames(Lioness_BALL_MP2_scaled))
Lioness_BALL_MP2_pivot <- Lioness_BALL_MP2_scaled %>% mutate(Edge = rownames(Lioness_BALL_MP2_scaled)) %>% 
  pivot_longer(!Edge, names_to = "Sample", values_to = "Z_score_Lioness") %>% 
  inner_join(assigments.TARGET, by = c("Sample" = "Sample"))

diff_MP2 <- Lioness_BALL_MP2_pivot %>% 
  group_by(Edge) %>% 
  dplyr::summarize(Median1.MP2 = median(Z_score_Lioness[consensuscluster == 1]),
                   Median2.MP2 = median(Z_score_Lioness[consensuscluster == 2]),
                   Diff_1vs2.MP2 = median(Z_score_Lioness[consensuscluster == 1]) - median(Z_score_Lioness[consensuscluster == 2]),
                   pval.Cls1_vs2.MP2 = wilcox.exact(Z_score_Lioness[consensuscluster == 1],
                                                    Z_score_Lioness[consensuscluster == 2])$p.value) %>% 
  mutate(adj.pval_1vs2.MP2 =  p.adjust(pval.Cls1_vs2.MP2, method = "fdr"))

write.csv(diff_MP2, "DiffCoExpr_MP2_25p.csv", quote = F, row.names = F)





