
#Let's check if taking random sets of 1508 edges of the Full GCN allows for clustering samples with different OS
library(rlist)
library(M3C)
library(dplyr)
library(ggsurvfit)
library(survival)

#Load SSN
Lioness_BALL_TARGET = read.csv("Lioness_BALL_TARGET.csv", row.names = 1)

#Make 1000 random sets of 1508 edges
nedges <- 1508
num.seq.edges <- 1:nrow(Lioness_BALL_TARGET)

set.seed(123)
for (i in 1:1000) {
  random_edges <- sample(num.seq.edges, nedges, replace = FALSE)
  
  if(i == 1){
    samples.edges <- list(random_edges)
  } else {
    samples.edges <- list.append(samples.edges, random_edges)
  }
}

saveRDS(samples.edges,"samplesEdges_1000reps.RDS")
samples.edges <- readRDS("samplesEdges_1000reps.RDS")
dir.create("KMs_Boostrap")

for(i in 1:length(samples.edges)){
  
print(i)

Lioness_BALL_TARGET_scaled <- as.data.frame(t(scale(t(Lioness_BALL_TARGET[samples.edges[[i]],]))))
  
r.hc_TARGET <- M3C(Lioness_BALL_TARGET_scaled, method=1, clusteralg = "hc", cores = 30, maxK = 5, seed = 123)

scores <- r.hc_TARGET$scores

K <- scores %>% filter(NORM_P < 0.05,
                       RCSI == max(RCSI)) %>% select(K) %>% as.numeric()

CI.BALL_TARGET <- read.delim("OS_files/CI_BALL.tsv")
assigments.TARGET <- r.hc_TARGET$realdataresults[[K]]$ordered_annotation
assigments.TARGET <- assigments.TARGET %>% mutate(Sample = rownames(assigments.TARGET))

CI.BALL_TARGET$Sample <- gsub("-", ".", CI.BALL_TARGET$Sample)
CI.BALL_TARGET <- CI.BALL_TARGET %>% inner_join(assigments.TARGET, by = c("Sample" = "Sample"))

SurvCI <- data.frame(Sample = CI.BALL_TARGET$patient,
                     Vital_status = CI.BALL_TARGET$VitalStatus,
                     DLF = CI.BALL_TARGET$OS_time,
                     consensuscluster = CI.BALL_TARGET$consensuscluster)

SurvCI$Vital_status[SurvCI$Vital_status == "Alive"] <- 0
SurvCI$Vital_status[SurvCI$Vital_status == "Dead"] <- 1
SurvCI$Vital_status <- as.numeric(SurvCI$Vital_status)

SurvCI$Vital_status[SurvCI$DLF > 1825] <- 0
SurvCI$DLF[SurvCI$DLF > 1825] <- 1825

OS_pval <- survdiff(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI)$pvalue

if(OS_pval < 0.05){
  
  surv.plot.BALL.lioness <- survfit2(Surv(DLF, Vital_status) ~ consensuscluster, data = SurvCI) %>%
    ggsurvfit() +
    labs(
      x = "Days to last follow up",
      y = "Overall survival probability"
    ) +
    add_risktable() + add_pvalue()
  
  ggsave(paste0("KMs_Boostrap/KM_TARGET_Boostrap_Rep", i, ".png"), surv.plot.BALL.lioness, height = 10, width = 7, dpi = 300)
  
}

Summary_i <- data.frame(Rep = i,
                       OptimalK = K,
                       SurvPval = OS_pval)

if(i == 1){
  SummaryGlobal <- Summary_i } else {
    SummaryGlobal <- rbind(SummaryGlobal, Summary_i)
  }

}

saveRDS(SummaryGlobal, "SummaryGlobal_Bootstrap.RDS")
library(ggplot2)
ggplot(SummaryGlobal, aes(x = SurvPval)) + geom_histogram()

