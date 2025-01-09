library(SummarizedExperiment)
library(TCGAbiolinks)
require(dplyr)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(tidyr)
library(rDNAse)
library(text.alignment)
library(Hmisc)
library(infotheo)
library(ggplot2)

###### Let's check sequence similarities between pairs of pseudogenes #####

# 1. Load clusters of pseudogenes and retrieve sequences and
# get annotation from rowData of TCGA

Clusters = read.csv("PSCls_in_BALL.csv")

mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

sequences = getSequence(id=Clusters[, 1], type="ensembl_gene_id", seqType="gene_exon_intron", mart = mart)

annot_TCGA$ensembl <- gsub("\\.[0-9]*", "", annot_TCGA$gene_id)

sequences <- sequences %>% mutate(Gene_name = annot_TCGA$gene_name[match(sequences$ensembl_gene_id, annot_TCGA$ensembl)])
sequences$Gene_name[218] <- "RPSAP58"
rownames(sequences) <- sequences$Gene_name

# 2. Get MI values

network <- read.csv("Network_PS_1_2_3_19_09.csv")

network <- network %>% inner_join(Clusters, by = c("Source" = "Gene"))

network$Source <- sequences$Gene_name[match(network$Source, sequences$ensembl_gene_id)]
network$Target <- sequences$Gene_name[match(network$Target, sequences$ensembl_gene_id)]

network$Similarity <- NA

for(i in 1:nrow(network)){
  
  sim <- smith_waterman(sequences[sequences$Gene_name == network$Source[i],1] , 
                        sequences[sequences$Gene_name == network$Target[i],1] 
                        , type="characters")
  
  network$Similarity[i] <- sim[["similarity"]]
  
}

MIvsSim <- ggplot(network, aes(MI, Similarity)) + geom_point()

rcorr(network$MI, network$Similarity, type = "pearson")
# x    y
# x 1.00 0.37
# y 0.37 1.00
rcorr(network$MI, network$Similarity, type = "spearman")
# x    y
# x 1.00 0.09
# y 0.09 1.00
infotheo::mutinformation(discretize(network$MI), discretize(network$Similarity))
# [1] 0.06756366

# ANSW: No correlation between sequence similarity and coexpression values. If PS 
# coexpression was driven by a bias in mapping one would expect:
# a) Coexpression between parental gene and PSs (occurs rarely)
# b) High correlation between COEXP vals and SeqSim (does not happen)
# c) COEXP patterns similar in NBM to those found in cancer data

##### Are pseudogenes being expressed also in the normal bone marrow dataset? #####

# 1. Let's load data

BALL_BM_tpms <- assay(BALL_BM_raw_primary, 4)
Normal_BoneMarrow_tpms <- assay(Normal_BoneMarrow, 4)

rownames(BALL_BM_tpms) <- rowData(BALL_BM_raw_primary)$gene_name
rownames(Normal_BoneMarrow_tpms) <- rowData(Normal_BoneMarrow)$gene_name

PS_in_network <- sequences$Gene_name[-which(sequences$Gene_name %in% annot[annot$Type == "protein_coding", "HGNC_symbol"])]

PS_in_network_TPMs_BALL <- BALL_BM_tpms[rownames(BALL_BM_tpms) %in% PS_in_network,]
PS_in_network_TPMs_NBM <- Normal_BoneMarrow_tpms[rownames(Normal_BoneMarrow_tpms) %in% PS_in_network,]

PS_in_network_TPMs_BALL_pivot <- PS_in_network_TPMs_BALL %>% as.data.frame %>% mutate(Gene = rownames(PS_in_network_TPMs_BALL)) %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "TPMs") %>% mutate(Phenotype = "BALL")

PS_in_network_TPMs_NBM_pivot <- PS_in_network_TPMs_NBM %>% as.data.frame %>% mutate(Gene = rownames(PS_in_network_TPMs_NBM)) %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "TPMs") %>% mutate(Phenotype = "Normal Bone Marrow")

# Make boxplots

ggplot(PS_in_network_TPMs_NBM_pivot, aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

ggplot(PS_in_network_TPMs_BALL_pivot, aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

#
ggplot(PS_in_network_TPMs_NBM_pivot[PS_in_network_TPMs_NBM_pivot$Gene %in% c("RPL31P7", "RPL31P2"),], aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

ggplot(PS_in_network_TPMs_BALL_pivot[PS_in_network_TPMs_BALL_pivot$Gene %in% c("RPL31P7", "RPL31P2"),], aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))
#

# 2. Is the expression of a pseudogene correlated with its similarity to the PC sequence?
# Let's make several examples of it, calculate similarities and plot a scatter of similarity and TPMs, e.g. RPL31

RPL31_TPMs_BALL <- BALL_BM_tpms[grep("RPL31", rownames(BALL_BM_tpms)),]
RPL31_TPMs_NBM <- Normal_BoneMarrow_tpms[grep("RPL31", rownames(Normal_BoneMarrow_tpms)),]

RPL31_TPMs_BALL_pivot <- RPL31_TPMs_BALL %>% as.data.frame %>% mutate(Gene = rownames(RPL31_TPMs_BALL)) %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "TPMs") %>% mutate(Phenotype = "BALL")

RPL31_TPMs_NBM_pivot <- RPL31_TPMs_NBM %>% as.data.frame %>% mutate(Gene = rownames(RPL31_TPMs_NBM)) %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "TPMs") %>% mutate(Phenotype = "Normal Bone Marrow")

ggplot(RPL31_TPMs_BALL_pivot[RPL31_TPMs_BALL_pivot$Gene != "RPL31",], aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

ggplot(RPL31_TPMs_NBM_pivot[RPL31_TPMs_NBM_pivot$Gene != "RPL31",], aes(Gene, TPMs, color = Gene)) + geom_boxplot(show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

#Get sequences

as.data.frame(rowData(BALL_BM_raw_primary))[grep("RPL31", rownames(BALL_BM_tpms)),]

ensembl_RPL31 <- annot_TCGA$ensembl[match(rownames(RPL31_TPMs_BALL), annot_TCGA$gene_name)]

sequences_RPL31 = getSequence(id=ensembl_RPL31, type="ensembl_gene_id", seqType="gene_exon_intron", mart = mart)

sequences_RPL31 <- sequences_RPL31 %>% mutate(Gene_name = annot_TCGA$gene_name[match(sequences_RPL31$ensembl_gene_id, annot_TCGA$ensembl)])

for(i in 1:nrow(sequences_RPL31)){
  
  sim <- smith_waterman(sequences_RPL31[sequences_RPL31$Gene_name == "RPL31",1], 
                        sequences_RPL31[i,1] 
                        , type="characters")
  
  SeqSim.DF <- data.frame(Gene1 = sequences_RPL31[sequences_RPL31$Gene_name == "RPL31",3],
                  Gene2 = sequences_RPL31[i,3],
                  SeqSim = sim[["similarity"]])
  
  if(i == 1) {
    x = SeqSim.DF
  } else {
    x <- rbind(SeqSim.DF, x)
  }
  
}

RPL31_means_BALL <- RPL31_TPMs_BALL %>% as.data.frame() %>% mutate(Gene = rownames(RPL31_TPMs_BALL)) %>% relocate(Gene) %>% 
  rowwise() %>% group_by(Gene) %>% summarise(Mean_BALL = mean(c_across(where(is.numeric))))

RPL31_means_NBM <- RPL31_TPMs_NBM %>% as.data.frame() %>% mutate(Gene = rownames(RPL31_TPMs_NBM)) %>% relocate(Gene) %>% 
  rowwise() %>% group_by(Gene) %>% summarise(Mean_NBM = mean(c_across(where(is.numeric))))

RPL31_corrs_inputs <- x %>% inner_join(RPL31_means_BALL, by = c("Gene2" = "Gene")) %>% inner_join(RPL31_means_NBM, by = c("Gene2" = "Gene"))

rcorr(RPL31_corrs_inputs$SeqSim[-27],  RPL31_corrs_inputs$Mean_BALL[-27], type = "pearson")
# x    y
# x 1.00 0.17
# y 0.17 1.00
# 
# n= 34 
# 
# 
# P
# x      y     
# x        0.3451
# y 0.3451    

rcorr(RPL31_corrs_inputs$SeqSim[-27],  RPL31_corrs_inputs$Mean_NBM[-27], type = "pearson")

# x    y
# x 1.00 0.05
# y 0.05 1.00
# 
# n= 34 
# 
# 
# P
# x      y     
# x        0.7717
# y 0.7717      

rcorr(RPL31_corrs_inputs$SeqSim[-27],  RPL31_corrs_inputs$Mean_BALL[-27], type = "spearman")
# x   y
# x 1.0 0.5
# y 0.5 1.0
# 
# n= 34 
# 
# 
# P
# x      y     
# x        0.0024
# y 0.0024    

rcorr(RPL31_corrs_inputs$SeqSim[-27],  RPL31_corrs_inputs$Mean_NBM[-27], type = "spearman")
# x    y
# x 1.00 0.31
# y 0.31 1.00
# 
# n= 34 
# 
# 
# P
# x      y     
# x        0.0755
# y 0.0755  

mutinformation(discretize(RPL31_corrs_inputs$SeqSim[-27]), discretize(RPL31_corrs_inputs$Mean_BALL[-27]))
#0.1216422
mutinformation(discretize(RPL31_corrs_inputs$SeqSim[-27]), discretize(RPL31_corrs_inputs$Mean_NBM[-27]))
#0.1353224

ggplot(RPL31_corrs_inputs[-27,], aes(Mean_BALL, SeqSim )) + geom_point() + geom_smooth(method="lm")

ggplot(RPL31_corrs_inputs[-27,], aes(Mean_NBM, SeqSim )) + geom_point() + geom_smooth(method="lm")

#Okay, nice test, now let's run it with multiple pseudogenes
#2.1 Identify pseudogenes and parental genes in data

PS_in_Network <- rownames(PS_in_network_TPMs_BALL)
#remove RNU6-8, SNHG9, HSP90AB3P
PS_in_Network <- PS_in_Network[-c(14, 22, 198)]

PG <- unique(gsub("P[^P]*$", "", PS_in_Network))
# PG <- readRDS("PG.RDS")
saveRDS(PG, "PG.RDS")

BALL_BM_tpms_rm_PS <- BALL_BM_tpms[!rownames(BALL_BM_tpms) == "RPSAP58",]
Normal_BoneMarrow_tpms_rm_PS <- Normal_BoneMarrow_tpms[!rownames(Normal_BoneMarrow_tpms) == "RPSAP58",]

# rownames(BALL_BM_tpms_rm_PS[grep(paste("^", PG[i], "$|", PG[i], "P",  "\\d", sep = ""), rownames(BALL_BM_tpms_rm_PS)),])
# 
# rownames(BALL_BM_tpms_rm_PS[grep("^TUBB$|TUBBP\\d+", rownames(BALL_BM_tpms_rm_PS)),])

#Use counts
# on server normalized_counts <- readRDS("TMM_deseq_nom_counts.RDS")
BALL_BM_tpms_rm_PS <- normalized_counts[,1:128]
Normal_BoneMarrow_tpms_rm_PS <- normalized_counts[,-c(1:128)]

BALL_BM_tpms_rm_PS <- BALL_BM_tpms_rm_PS[-which(rownames(BALL_BM_tpms_rm_PS) %in% c("RPSAP58","HSP90AB3P")),]
Normal_BoneMarrow_tpms_rm_PS <- Normal_BoneMarrow_tpms_rm_PS[-which(rownames(Normal_BoneMarrow_tpms_rm_PS) %in% c("RPSAP58","HSP90AB3P")),]

###### Automatization #####
dir.create("/home/anakamura/SSN_pseudogenes/Redo_in_Oslo/scatters_SeqSim_vs_MeanExprCounts")
dir.scatters <- "/home/anakamura/SSN_pseudogenes/Redo_in_Oslo/scatters_SeqSim_vs_MeanExprCounts"

dir.create("/home/anakamura/SeqSims/scatters_SeqSim_vs_MeanExprCounts")
dir.scatters <- "/home/anakamura/SeqSims/scatters_SeqSim_vs_MeanExprCounts"
annot_TCGA<- readRDS("annot_TCGA_use.RDS") 
PG <- readRDS("PG.RDS")

for(i in 55:length(PG)){
  
  print(i)
  print(paste("Analayzing", PG[i]))
  
  PG_TPMs_BALL <- BALL_BM_tpms_rm_PS[grep(paste("^", PG[i], "$|", PG[i], "P",  "\\d+", sep = ""), rownames(BALL_BM_tpms_rm_PS)),]
  PG_TPMs_NBM <- Normal_BoneMarrow_tpms_rm_PS[grep(paste("^", PG[i], "$|", PG[i], "P",  "\\d+", sep = ""), rownames(Normal_BoneMarrow_tpms_rm_PS)),]
  
  PG_TPMs_BALL_pivot <- PG_TPMs_BALL %>% as.data.frame %>% mutate(Gene = rownames(PG_TPMs_BALL)) %>% 
    pivot_longer(!Gene, names_to = "Samples", values_to = "Counts") %>% mutate(Phenotype = "BALL")
  
  PG_TPMs_NBM_pivot <- PG_TPMs_NBM %>% as.data.frame %>% mutate(Gene = rownames(PG_TPMs_NBM)) %>% 
    pivot_longer(!Gene, names_to = "Samples", values_to = "Counts") %>% mutate(Phenotype = "Normal Bone Marrow")
  
  ensembl_PG <- annot_TCGA$ensembl[match(rownames(PG_TPMs_BALL), annot_TCGA$gene_name)]
  
  print("Getting sequences")
  sequences_PG = getSequence(id=ensembl_PG, type="ensembl_gene_id", seqType="gene_exon_intron", mart = mart)
  
  sequences_PG <- sequences_PG %>% mutate(Gene_name = annot_TCGA$gene_name[match(sequences_PG$ensembl_gene_id, annot_TCGA$ensembl)])
  
  print("Getting SeqSims")
  for(j in 1:nrow(sequences_PG)){
    
    sim <- smith_waterman(sequences_PG[sequences_PG$Gene_name == PG[i],1], 
                          sequences_PG[j,1] 
                          , type="characters")
    
    SeqSim.DF <- data.frame(Gene1 = sequences_PG[sequences_PG$Gene_name == PG[i],3],
                            Gene2 = sequences_PG[j,3],
                            SeqSim = sim[["similarity"]])
    
    if(j == 1) {
      x = SeqSim.DF
    } else {
      x <- rbind(SeqSim.DF, x)
    }
  }
  
  print("Make scatter")
  PG_means_BALL <- PG_TPMs_BALL %>% as.data.frame() %>% mutate(Gene = rownames(PG_TPMs_BALL)) %>% relocate(Gene) %>% 
    rowwise() %>% group_by(Gene) %>% summarise(Mean_BALL = mean(c_across(where(is.numeric))))
  
  PG_means_NBM <- PG_TPMs_NBM %>% as.data.frame() %>% mutate(Gene = rownames(PG_TPMs_NBM)) %>% relocate(Gene) %>% 
    rowwise() %>% group_by(Gene) %>% summarise(Mean_NBM = mean(c_across(where(is.numeric))))
  
  PG_corrs_inputs <- x %>% inner_join(PG_means_BALL, by = c("Gene2" = "Gene")) %>% inner_join(PG_means_NBM, by = c("Gene2" = "Gene"))
  
  PG_corrs_inputs <- PG_corrs_inputs[-which(PG_corrs_inputs$Gene1 == PG_corrs_inputs$Gene2),]
  
  scatter_BALL <- ggplot(PG_corrs_inputs, aes(Mean_BALL, SeqSim )) + geom_point() + geom_smooth(method="lm") +
    labs(title = paste(PG[i]),
         subtitle = "BALL")

  ggsave(filename =  paste0(dir.scatters, "/ScatterPlotCounts_",PG[i], "_BALL.png"), plot = scatter_BALL, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
  scatter_NBM <- ggplot(PG_corrs_inputs, aes(Mean_NBM, SeqSim )) + geom_point() + geom_smooth(method="lm") +
    labs(title = paste(PG[i]),
         subtitle = "Normal bone marrow")
  
  ggsave(filename =  paste0(dir.scatters, "/ScatterPlotCounts_",PG[i], "_NBM.png"), plot = scatter_NBM, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
  if(i == 1) {
    Full_PG_corrs_inputs = PG_corrs_inputs
  } else {
    Full_PG_corrs_inputs <- rbind(Full_PG_corrs_inputs, PG_corrs_inputs)
  }
  
  if(length(PG_corrs_inputs$Mean_BALL) > 4) {
    print("Making correlations")
    Corr_BALL <- rcorr(PG_corrs_inputs$Mean_BALL, PG_corrs_inputs$SeqSim, type = "spearman")
    Corr_NBM <- rcorr(PG_corrs_inputs$Mean_NBM, PG_corrs_inputs$SeqSim, type = "spearman")
    
    MI_BALL <- mutinformation(discretize(PG_corrs_inputs$Mean_BALL), discretize(PG_corrs_inputs$SeqSim))
    MI_NBM <- mutinformation(discretize(PG_corrs_inputs$Mean_NBM), discretize(PG_corrs_inputs$SeqSim))
    
    corr_df <- data.frame(Gene = PG[i],
                          CorrV_BALL = Corr_BALL$r[1,2],
                          CorrP_BALL = Corr_BALL$P[1,2],
                          CorrV_NBM = Corr_NBM$r[1,2],
                          CorrP_NBM = Corr_NBM$P[1,2],
                          MiV_BALL = MI_BALL,
                          MiV_NBM = MI_NBM)
    
    if(i == 1) {
      Full_corr_df = corr_df
    } else {
      Full_corr_df <- rbind(Full_corr_df, corr_df)
    } 
  } else {
    print(paste(PG[i],"dosen't have enough paired pseudogenes for correlation estimation"))
  }
}

saveRDS(Full_corr_df, "Full_corr_df.RDS")
saveRDS(Full_PG_corrs_inputs, "Full_PG_corrs_inputs.RDS")

Full_corr_df <- readRDS("Full_corr_df.RDS")
Full_PG_corrs_inputs <- readRDS("Full_PG_corrs_inputs.RDS")

#Analyze results
Relevant_corr_BALL <- Full_corr_df[Full_corr_df$CorrP_BALL < 0.05 & 
                                     Full_corr_df$CorrV_BALL >= 0.5,]

Relevant_corr_NBM <- Full_corr_df[Full_corr_df$CorrP_NBM < 0.05 & 
                                     Full_corr_df$CorrV_NBM >= 0.5,]

Global_relevant <- intersect(Relevant_corr_BALL$Gene, Relevant_corr_NBM$Gene)

Exclusive_BALL <- Relevant_corr_BALL$Gene[!(Relevant_corr_BALL$Gene %in% Relevant_corr_NBM$Gene)]

Exclusive_NBM <- Relevant_corr_NBM$Gene[!(Relevant_corr_NBM$Gene %in% Relevant_corr_BALL$Gene)]

Full_corr_df$Gene[!(Full_corr_df$Gene %in% Relevant_corr_BALL$Gene & 
  Full_corr_df$Gene %in% Relevant_corr_NBM$Gene)]

Full_PG_corrs_inputs$num <- rownames(Full_PG_corrs_inputs)
Full_PG_corrs_inputs <- Full_PG_corrs_inputs[!(Full_PG_corrs_inputs$num %in% c(139,425,1212,158,180,340,1137,1421,935,1177,4123,5117,1181)),-6]

dir.create("Relevant_plots")
for(i in 1:length(Global_relevant)){
  
  gene <- Global_relevant[i]
  
  corrs <- Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == gene, ]
  
  nbm_plot <- ggplot(corrs, 
         aes(Mean_NBM, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
       geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "NBM") + geom_text(angle = 45)
  
  ggsave(filename =  paste0( gene, "_NBM.png"), plot = nbm_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
  ball_plot <- ggplot(corrs, 
                     aes(Mean_BALL, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "BALL") + geom_text(angle = 45)
  
  ggsave(filename =  paste0(gene, "_BALL.png"), plot = ball_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
}

for(i in 1:length(Exclusive_BALL)){
  
  gene <- Exclusive_BALL[i]
  
  corrs <- Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == gene, ]
  
  nbm_plot <- ggplot(corrs, 
                     aes(Mean_NBM, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "NBM") + geom_text(angle = 45)
  
  ggsave(filename =  paste0(gene, "_NBM.png"), plot = nbm_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
  ball_plot <- ggplot(corrs, 
                      aes(Mean_BALL, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "BALL") + geom_text(angle = 45)
  
  ggsave(filename =  paste0(gene, "_BALL.png"), plot = ball_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
}

for(i in 1:length(Exclusive_NBM)){
  
  gene <- Exclusive_NBM[i]
  
  corrs <- Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == gene, ]
  
  nbm_plot <- ggplot(corrs, 
                     aes(Mean_NBM, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "NBM") + geom_text(angle = 45)
  
  ggsave(filename =  paste0( gene, "_NBM.png"), plot = nbm_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
  ball_plot <- ggplot(corrs, 
                      aes(Mean_BALL, SeqSim, label = corrs[,2] )) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = gene,
         subtitle = "BALL") + geom_text(angle = 45)
  
  ggsave(filename =  paste0(gene, "_BALL.png"), plot = ball_plot, units="in",
         width=10, height=10, dpi =300, bg = "white")
  
}

rcorr(Full_PG_corrs_inputs$Mean_BALL[Full_PG_corrs_inputs$Gene1 == "RPL10A"]
      , Full_PG_corrs_inputs$SeqSim[Full_PG_corrs_inputs$Gene1 == "RPL10A"], type = "pearson")


ggplot(Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == "RPL10A",], aes(SeqSim)) + geom_histogram() +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_PG_corrs_inputs$Mean_BALL[Full_PG_corrs_inputs$Gene1 == "RPL10A"]), sd = sd(Full_PG_corrs_inputs$SeqSim[Full_PG_corrs_inputs$Gene1 == "RPL10A"])))

ggplot(Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == "RPL10A",], aes(SeqSim)) + geom_histogram() +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_PG_corrs_inputs$SeqSim[Full_PG_corrs_inputs$Gene1 == "RPL10A"]), sd = sd(Full_PG_corrs_inputs$SeqSim[Full_PG_corrs_inputs$Gene1 == "RPL10A"])))

ggplot(Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == "RPL10A",], aes(Mean_NBM)) + geom_histogram() +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_PG_corrs_inputs$Mean_NBM[Full_PG_corrs_inputs$Gene1 == "RPL10A"]), sd = sd(Full_PG_corrs_inputs$Mean_NBM[Full_PG_corrs_inputs$Gene1 == "RPL10A"])))


ggplot(Full_PG_corrs_inputs[Full_PG_corrs_inputs$Gene1 == "RPL10A",], aes(Mean_BALL)) + geom_histogram() +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_PG_corrs_inputs$Mean_BALL[Full_PG_corrs_inputs$Gene1 == "RPL10A"]), sd = sd(Full_PG_corrs_inputs$Mean_BALL[Full_PG_corrs_inputs$Gene1 == "RPL10A"])))

ggplot(RPL31_corrs_inputs[-27,], aes(Mean_NBM)) + geom_histogram() +
stat_function(fun = dnorm, lwd = 2, col = 'red',
args = list(mean = mean(RPL31_corrs_inputs$Mean_NBM[-27]), sd = sd(RPL31_corrs_inputs$Mean_NBM[-27])))

ggplot(RPL31_corrs_inputs[-27,], aes(SeqSim)) + geom_histogram() +
stat_function(fun = dnorm, lwd = 2, col = 'red',
args = list(mean = mean(RPL31_corrs_inputs$SeqSim[-27]), sd = sd(RPL31_corrs_inputs$SeqSim[-27])))

Hist_full.BALL <- ggplot(Full_corr_df, aes(CorrV_BALL)) + geom_histogram(binwidth = 0.05) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_corr_df$CorrV_BALL), sd = sd(Full_corr_df$CorrV_BALL)))

Hist_full.NBM <- ggplot(Full_corr_df, aes(CorrV_NBM)) + geom_histogram(binwidth = 0.05) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Full_corr_df$CorrV_NBM), sd = sd(Full_corr_df$CorrV_NBM)))

#Find out number of pseudogenes per parental gene

for(i in 1:length(Full_corr_df$Gene)) {
  
  gene <- Full_corr_df$Gene[i]
  PG_TPMs_BALL <- BALL_BM_tpms_rm_PS[grep(paste("^", gene, "$|", gene, "P",  "\\d+", sep = ""), rownames(BALL_BM_tpms_rm_PS)),]
  
  if(length(rownames(PG_TPMs_BALL)) >= 10){
    if(!exists("relevant_genes")){
      relevant_genes <- gene
    } else {
      relevant_genes <- c(relevant_genes,gene)
    }
  }
}

#Genes with more than 10 genes
Filtered_Full_corr_df <- Full_corr_df[Full_corr_df$Gene %in% relevant_genes,]

Hist_RevGenes_pval.BALL <- ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$CorrP_BALL <= 0.05,], aes(CorrV_BALL)) + 
  geom_histogram(binwidth = 0.15) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Filtered_Full_corr_df$CorrV_BALL), sd = sd(Filtered_Full_corr_df$CorrV_BALL)))

Hist_RevGenes_pval.NBM <- ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$CorrP_NBM <= 0.05,], aes(CorrV_NBM)) + 
  geom_histogram(binwidth = 0.15) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Filtered_Full_corr_df$CorrV_NBM), sd = sd(Filtered_Full_corr_df$CorrV_NBM)))

#Now adj pval
Filtered_Full_corr_df <- Filtered_Full_corr_df %>% 
  mutate(adj.pval.BALL = p.adjust(Filtered_Full_corr_df$CorrP_BALL, method = "fdr", n = 2*length(Filtered_Full_corr_df$CorrP_BALL))) %>%
  mutate(adj.pval.NBM = p.adjust(Filtered_Full_corr_df$CorrP_NBM, method = "fdr", n = 2*length(Filtered_Full_corr_df$CorrP_NBM)))

Hist_RevGenes_pvalAdj.BALL <- ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$adj.pval.BALL <= 0.05,], aes(CorrV_BALL)) + geom_histogram(binwidth = 0.01) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Filtered_Full_corr_df$CorrV_BALL), sd = sd(Filtered_Full_corr_df$CorrV_BALL)))

Hist_RevGenes_pvalAdj.NBM = ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$adj.pval.NBM <= 0.05,], aes(CorrV_NBM)) + geom_histogram(binwidth = 0.01) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Filtered_Full_corr_df$CorrV_NBM), sd = sd(Filtered_Full_corr_df$CorrV_NBM)))

#Cut based on CorrV

  ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$CorrV_BALL >= 0.5,], aes(CorrV_BALL)) + geom_histogram(binwidth = 0.05) +
  stat_function(fun = dnorm, lwd = 2, col = 'red',
                args = list(mean = mean(Filtered_Full_corr_df$CorrV_BALL), sd = sd(Filtered_Full_corr_df$CorrV_BALL)))

  ggplot(Filtered_Full_corr_df[Filtered_Full_corr_df$CorrV_NBM >= 0.5,], aes(CorrV_NBM)) + geom_histogram(binwidth = 0.05) +
    stat_function(fun = dnorm, lwd = 2, col = 'red',
                  args = list(mean = mean(Filtered_Full_corr_df$CorrV_NBM), sd = sd(Filtered_Full_corr_df$CorrV_NBM)))
  
  
  
