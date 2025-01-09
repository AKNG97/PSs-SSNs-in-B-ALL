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

annotPS <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/annot_PS.RDS")
annot <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/annot_TCGA.RDS")
Cancer_GCN <- Cancer_GCN %>% filter(Source %in% annotPS$HGNC_symbol &
                                      Target %in% annotPS$HGNC_symbol)

Lioness_BALL_TARGET1 = Lioness_BALL_TARGET[!(rownames(Lioness_BALL_TARGET) %in%
                                              Cancer_GCN$ID),]

Lioness_BALL_TARGET2 = Lioness_BALL_TARGET[rownames(Lioness_BALL_TARGET) %in%
                                               Cancer_GCN$ID,]

filtered_results_TARGET <- featurefilter(Lioness_BALL_TARGET1, percentile=1508*100/39082,
                                         method='var', topN=5)

filtered_results_TARGET <- featurefilter(Lioness_BALL_TARGET2, percentile=25,
                                         method='var', topN=5)

Lioness_BALL_TARGET_scaled <- as.data.frame(t(scale(t(filtered_results_TARGET$filtered_data))))

r.hc_TARGET <- M3C(Lioness_BALL_TARGET_scaled, method=1, clusteralg = "hc", cores = 30, maxK = 10, seed = 123)

ggsave("pval_M3C_BALL_Target_PS_paper.png", plot = r.hc_TARGET$plots[[3]], dpi = 300)
ggsave("RCSI_M3C_BALL_Target_PS_paper_2.png", plot = r.hc_TARGET$plots[[4]], dpi = 300)

CI.BALL_TARGET <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/CI_primaryBM_BALL.RDS")
#TARGET
assigments.TARGET <- r.hc_TARGET$realdataresults[[2]]$ordered_annotation
assigments.TARGET <- assigments.TARGET %>% mutate(Sample = rownames(assigments.TARGET))
#saveRDS(assigments.TARGET, "Final/assigmentsTARGET_noPS_1508.RDS")

# 
CI.BALL_TARGET$patient <- gsub("-", ".", CI.BALL_TARGET$patient)
# #ci.all$patient <- gsub("-", ".", ci.all$patient)
CI.BALL_TARGET <- CI.BALL_TARGET %>% inner_join(assigments.TARGET, by = c("patient" = "Sample"))

SurvCI <- data.frame(Sample = CI.BALL_TARGET$patient,
                     Vital_status = CI.BALL_TARGET$VS_update,
                     DLF = CI.BALL_TARGET$OS_time_update,
                     consensuscluster.Full = CI.BALL_TARGET$consensuscluster)

SurvCI <- data.frame(Sample = CI.BALL_TARGET$patient,
                     Vital_status = CI.BALL_TARGET$VS_update,
                     DLF = CI.BALL_TARGET$OS_time_update,
                     consensuscluster = CI.BALL_TARGET$consensuscluster)

SurvCI$Vital_status[SurvCI$Vital_status == "Alive"] <- 0
SurvCI$Vital_status[SurvCI$Vital_status == "Dead"] <- 1
SurvCI$Vital_status <- as.numeric(SurvCI$Vital_status)

SurvCI$Vital_status[SurvCI$DLF > 1825] <- 0
SurvCI$DLF[SurvCI$DLF > 1825] <- 1825

survfit_object <- survdiff(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI)
# 
# survminer::pairwise_survdiff(Surv(DLF, Vital_status) ~ consensuscluster,
#                             data = SurvCI)

print(survfit_object$strata)

colores <- c("consensuscluster=1" = "#619CFF", 
             "consensuscluster=2" = "#F8766D", 
             "consensuscluster=3" = "#00BA38")
surv.plot.BALL.lioness <- survfit2(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI) %>%
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

ggsave("KM_TARGET_Full_paper_2.png", surv.plot.BALL.lioness, height = 10, width = 7, dpi = 300)


default_colors <- ggplot_build(surv.plot.BALL.lioness)$data[[1]]$colour
print(default_colors)

# 

##### MP2 ##### 
#/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/new_analysis/analyze2
#module load R/4.2.1-foss-2022a

Lioness_BALL_MP2 = data.table::fread("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/LionessPS_BALL_MP2_Full_CheckError.csv")
Lioness_BALL_MP2 <- Lioness_BALL_MP2 %>% as.data.frame
rownames(Lioness_BALL_MP2) <- Lioness_BALL_MP2$Edge_ID
Lioness_BALL_MP2 <- Lioness_BALL_MP2[,-1]

Cancer_GCN <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/BALL_MP2_GCN_SignifEdges_2ndRun.RDS")

annotPS <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/annot_PS.RDS")
annot <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/annot_TCGA.RDS")
Cancer_GCN1 <- Cancer_GCN %>% filter(Source %in% annotPS$HGNC_symbol &
                                       Target %in% annotPS$HGNC_symbol)

Lioness_BALL_MP21 = Lioness_BALL_MP2[!(rownames(Lioness_BALL_MP2) %in%
                                         Cancer_GCN1$ID),]

Lioness_BALL_MP23 = Lioness_BALL_MP2[rownames(Lioness_BALL_MP2) %in%
                                       Cancer_GCN1$ID,]

filtered_results_MP2 <- featurefilter(Lioness_BALL_MP23, percentile=4811*100/19244, method='var', topN=5)
filtered_results_MP2 <- featurefilter(Lioness_BALL_MP21, percentile=4811*100/505650, method='var', topN=5)
filtered_results_MP2 <- featurefilter(Lioness_BALL_MP2, percentile=4811*100/524894, method='var', topN=5)

Lioness_BALL_MP23_scaled <- as.data.frame(t(scale(t(filtered_results_MP2$filtered_data))))
#Lioness_BALL_MP2_scaled <- as.data.frame(t(scale(t(Lioness_BALL_MP21))))
# 
# Lioness_BALL_MP23_scaled <- as.data.frame(t(scale(t(Lioness_BALL_MP23))))
# 
# filtered_results_MP21 <- featurefilter(Lioness_BALL_MP21, percentile=2000*100/505650, 
#                                        method='var', topN=5)
#Lioness_BALL_MP21_scaled <- as.data.frame(t(scale(t(filtered_results_MP21$filtered_data))))

r.hc_MP2 <- M3C(Lioness_BALL_MP23_scaled, method=1, clusteralg = "hc", cores = 70, maxK = 10, seed = 123)
#r.hc_MP2 <- M3C(filtered_results_MP2$filtered_data, method=1, clusteralg = "spectral", cores = 80, maxK = 5, seed = 123)

CI.BALL_MP2 <- readRDS("/storage/kuijjerarea/akng/PseudogenesAnalysis/HandleReps/ExtCI/MP2PRT/CorrectCI/Full_Correct/CI_MP2PRT_FULL.RDS")
assigments.MP2 <- r.hc_MP2$realdataresults[[2]]$ordered_annotation
assigments.MP2 <- assigments.MP2 %>% mutate(Sample = rownames(assigments.MP2))
#saveRDS(assigments.MP2, "Final/assigmentsMP2_Full_4811.RDS")
saveRDS(SurvCI, "assigmentsMP2_PS_25p.RDS")

#KM
CI.BALL_MP2$sample_submitter_id <- gsub("-", ".", CI.BALL_MP2$sample_submitter_id)
CI.BALL_MP2 <- CI.BALL_MP2 %>% inner_join(assigments.MP2, by = c("sample_submitter_id" = "Sample"))

SurvCI <- data.frame(Sample = CI.BALL_MP2$sample_submitter_id,
                     Vital_status = CI.BALL_MP2$OS_Status,
                     DLF = CI.BALL_MP2$OS_Time,
                     consensuscluster = CI.BALL_MP2$consensuscluster)

SurvCI$Vital_status[SurvCI$DLF > 1825] <- 0
SurvCI$DLF[SurvCI$DLF > 1825] <- 1825

survdiff(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI)

surv.plot.BALL.lioness <- survfit2(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI) %>% 
  ggsurvfit(size = 1.5) +
  labs(
    x = "Days to last follow up",
    y = "Overall survival probability"
  ) + add_pvalue() + theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=14)) +
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

ggsave("KM_MP2PRT_Full_4811.png", surv.plot.BALL.lioness, height = 10, width = 7, dpi = 300)

