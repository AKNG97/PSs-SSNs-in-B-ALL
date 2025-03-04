#This script shows the analysis of sequence similarity between pseudogenes and their parental genes
#using the data from TARGET-ALL-P2

library(SummarizedExperiment)
library(TCGAbiolinks)
require(dplyr)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(tidyr)
library(text.alignment)
library(Hmisc)
library(infotheo)
library(ggplot2)
library(igraph)
library(parallel)

#### Let's check if highly connected genes have high SeqSim between pairs ####
#1. Load GCN of PS
GCN <- readRDS("GCN_BALL_TARGET.RDS")

annot <- read.delim("gene_annotation_TCGA.tsv")

annotPS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]

#Filter GCN to select Pseudogene-Pseudogene edges

GCN <- GCN %>% filter(Source %in% annotPS$HGNC_symbol &
                                      Target %in% annotPS$HGNC_symbol)

#2. Get Communities

g <- graph_from_data_frame(GCN[,1:3], directed = FALSE)
set.seed(123)
lc <- cluster_louvain(g)

#3. Filter Communities with more than 10 elements
co <- unlist(communities(lc)) 
co <- communities(lc)[sapply(communities(lc), function(x) length(x) > 10)]

#4. Get Sequences
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version = 110)

annotPS <- annotPS[annotPS$HGNC_symbol %in% unique(c(GCN$Source, GCN$Target)), ]

sequences = getSequence(id=annotPS$Ensembl, type="ensembl_gene_id", seqType="gene_exon_intron", mart = mart)

sequences[,2] <- annotPS$HGNC_symbol[match(sequences[,2], annotPS$Ensembl)]

colnames(sequences) <- c("Seq", "GeneSymbol")

#5. Correlate SeqSim to EdgeWeight

network <- GCN[,1:3]
network$Similarity <- NA

for(i in 1:nrow(network)){
  
  print(i)
  
  sim <- smith_waterman(sequences[sequences$GeneSymbol == network$Source[i],1] , 
                        sequences[sequences$GeneSymbol == network$Target[i],1] 
                        , type="characters")
  
  network$Similarity[i] <- sim[["similarity"]]
  
}

#Function to calculate SeqSim
calculate_similarity <- function(i) {
  sim <- smith_waterman(sequences[sequences$GeneSymbol == network$Source[i], 1], 
                        sequences[sequences$GeneSymbol == network$Target[i], 1],
                        type = "characters")
  return(sim[["similarity"]])
}

similarity_results <- mclapply(1:nrow(network),
                               function(i) calculate_similarity(i), 
                               mc.cores = 20)

network$Similarity <- unlist(similarity_results)

saveRDS(network, "Network_SeqSimPairs.RDS")

rcorr(abs(network$Sp_corr), network$Similarity)
# x    y
# x 1.00 0.24
# y 0.24 1.00
# 
# n= 6032 
# 
# P
# x  y 
# x     0
# y  0   

GlobalCorrelation <- ggplot(network, aes(x = abs(Sp_corr), y = Similarity)) + 
  labs(x = "Edge Weight", y = "Sequence similarity", 
       title = paste0("Correlation between edge weight and sequence similarity, ρ = 0.24, p = 0"),
       subtitle = "Complete network") +
  geom_point()

dir.create("Plots")
ggsave("Plots/GlobalCorrelation_EdgeWeight_vs_SeqSim.png", GlobalCorrelation, height = 10, width = 10, dpi = 300)

saveRDS(co, "Communities_PS_TARGET_Network.RDS")

#6. Correlate SeqSim to EdgeWeight by Community

for(i in 1:length(co)){
  
  Comm <- network[network$Source %in% co[[i]] &
                    network$Target %in% co[[i]], ]
  
  correlation <- round(rcorr(Comm$Sp_corr, Comm$Similarity)$r[1,2], 2)
  pvalue <- round(rcorr(Comm$Sp_corr, Comm$Similarity)$P[1,2], 2)
  
  CorrelationPlot <- ggplot(Comm, aes(x = abs(Sp_corr), y = Similarity)) + 
    labs(x = "Edge Weight", y = "Sequence similarity", 
         title = paste0("Correlation between edge weight and sequence similarity, ", 
                                  "ρ = ", correlation, ", p = ",  pvalue),
         subtitle = paste0("Community ", i)) +
    geom_point()

  ggsave(paste0("Plots/Correlation_EdgeWeight_vs_SeqSim_Comm", i, ".png"), CorrelationPlot, 
         height = 10, width = 10, dpi = 300)
  
}


#### Now let's analyze pseudogenes families and see if their expression is correlated to their SeqSim with the Parental Gene ####

ALL <- readRDS("RawData/ALL_TARGET.RDS")

dim(ALL)
#[1] 60660   532
which(is.na(ALL$sample_type))
#integer(0)
which(is.na(ALL$definition))
#integer(0)
table(as.factor(ALL$primary_diagnosis))
# Acute lymphocytic leukemia 
# 532

table(as.factor(ALL$definition))
# Primary Blood Derived Cancer - Bone Marrow 
# 387 
# Primary Blood Derived Cancer - Peripheral Blood 
# 76 
# Recurrent Blood Derived Cancer - Bone Marrow 
# 68 
# Recurrent Blood Derived Cancer - Peripheral Blood 
# 1 

ALL_primaryBM <- ALL[,ALL$definition == "Primary Blood Derived Cancer - Bone Marrow"]

CI.ALL <- as.data.frame(colData(ALL_primaryBM))

#Load Clinical Information from TCGA
CI_1_ALL <- readxl::read_xlsx("ClinicalFiles/ALL_CI_1.xlsx")
CI_2_ALL <- readxl::read_xlsx("ClinicalFiles/ALL_CI_2.xlsx")
CI_3_ALL <- readxl::read_xlsx("ClinicalFiles/ALL_CI_3.xlsx")
CI_4_ALL <- readxl::read_xlsx("ClinicalFiles/ALL_CI_4.xlsx")
CI_5_ALL <- readxl::read_xlsx("ClinicalFiles/ALL_CI_5.xlsx")

CI_full_ALL <- rbind(CI_1_ALL, CI_2_ALL, CI_3_ALL, CI_4_ALL, CI_5_ALL)
table(CI_full_ALL$`Cell of Origin`)
# B Cell ALL B precursor (Non-T, Non-B ALL) 
# 517                             20 
# B-Precursor                     T Cell ALL 
# 271                            269 
table(CI_full_ALL$`Vital Status`)
# Alive    Dead Unknown 
# 870     195       7 

CI_full_ALL <- CI_full_ALL %>% dplyr::rename(VS_update = `Vital Status`,
                                             OS_time_update = `Overall Survival Time in Days`,
                                             patient = `TARGET USI`)
CI_full_ALL <- CI_full_ALL[!duplicated(CI_full_ALL$patient),]
CI_full_BALL <- CI_full_ALL[CI_full_ALL$`Cell of Origin` %in% c("B Cell ALL", "B-Precursor"),]

CI_BALL <- CI.ALL %>% inner_join(CI_full_BALL[,c("patient", "VS_update", "OS_time_update")],
                                 by = c("patient" = "patient"))
table(CI_BALL$VS_update)
# Alive    Dead Unknown 
# 74      65       3 

which(CI_BALL$VS_update == "Unknown")
#[1]   8  24 102

CI_BALL <- CI_BALL[,c("patient", "barcode", "vital_status", "days_to_death", "days_to_last_follow_up", "VS_update", "OS_time_update")]
CI_BALL[which(CI_BALL$VS_update == "Unknown"),]
# patient                  barcode vital_status days_to_death
# 8   TARGET-10-PANWEZ TARGET-10-PANWEZ-09B-01R        Alive            NA
# 24  TARGET-10-PAPHYN TARGET-10-PAPHYN-09B-01R        Alive            NA
# 102 TARGET-10-PASIZE TARGET-10-PASIZE-09B-01R      Unknown            NA
# days_to_last_follow_up VS_update OS_time_update
# 8                       NA   Unknown           3145
# 24                      NA   Unknown           3012
# 102                     NA   Unknown           1699

CI_BALL$VS_update[which(CI_BALL$VS_update == "Unknown")] <- CI_BALL$vital_status[which(CI_BALL$VS_update == "Unknown")]
table(CI_BALL$VS_update)
# Alive    Dead Unknown 
# 76      65       1  

CI_BALL <- CI_BALL[CI_BALL$VS_update %in% c("Alive", "Dead"),]
table(CI_BALL$VS_update)
# Alive  Dead 
# 76    65 
CI_BALL$primary_diagnosis = "BALL"
dim(CI_BALL)
#[1] 141   8

BALL <- ALL_primaryBM[,ALL_primaryBM$patient %in% CI_BALL$patient]
dim(BALL)
#[1] 60660   141

#Let's check for replicates
which(duplicated(BALL$patient))
#[1]  5  16  21  23  76  83 120 125 126

#Prepare to add them with DESeq2

which((1:ncol(BALL) == match(colnames(BALL), CI_BALL$barcode)) == FALSE)
#integer(0). Hence, samples in expression matrix and metadata are in the same order

BALL_TPM <- assay(BALL, 4)

annot <- read.delim("gene_annotation_TCGA.tsv")

annotPS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]
                             
BALL_TPM <- BALL_TPM[rownames(BALL_TPM) %in% annot$gene_id, ]
rownames(BALL_TPM) <- annot$HGNC_symbol[match(rownames(BALL_TPM), annot$gene_id)]

PS_in_Network <- rownames(BALL_TPM[rownames(BALL_TPM) %in% unique(c(GCN$Source, GCN$Target)), ])

PG <- unique(gsub("P[^P]*$", "", PS_in_Network))

#We need to run this 124 families were excluded for a mistake in line 239

mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version = 110)

dir.scatters <- "Scatter_SeqSim"

dir.create(dir.scatters)

#We will skip 31 = H3, don't know which is the parental gene
#MT-CO3 (i=65) is not present on the expression matrix, let's skip
#MT-ND4L (i = 66) is not present on the expression matrix, let's skip
#RNA5S ( i = 67) is not present on the expression matrix, let's skip
#103 = TUBA, which referes to tublin alpha, let's use as symbol TUBA1A
#Problem with 105
#134 = MTCO1, symbol is not correct, could be MT-CO1, but it is not present either, let's skip
#141 = MTCYB, symbol is not correct, could be MT-CYB, but it is not present either, let's skip
#151 = MTND2, symbol is not correct and is not in annot, let's skip
#153 = MTND4, symbol is not correct and is not in annot, let's skip
#155 = MTND1, symbol is not correct and is not in annot, let's skip
#171 = MTND6, symbol is not correct and is not in annot, let's skip
#173 = MTATP6, symbol is not correct and is not in annot, let's skip
#187 = MTCO2, symbol is not correct and is not in annot, let's skip
#190 = MTND5, symbol is not correct and is not in annot, let's skip

PG[103] <- "TUBA1A"

for(i in 191:length(PG)){
  
  print(i)
  print(paste("Analayzing", PG[i]))
  
if(!(i %in% c(31, 65, 66, 67, 134)))  {
  
  PG_TPMs_BALL <- BALL_TPM[c(grep(paste("^", PG[i], "P\\d+$", sep = ""), rownames(BALL_TPM)),
                             which(rownames(BALL_TPM) == PG[i])),]
  # PG_TPMs_BALL <- BALL_TPM[c(grep(paste("^", PG[i], "P\\d+$", sep = ""), rownames(BALL_TPM)),
  #                            which(rownames(BALL_TPM) == "TUBA1A")),]
  # 
  
 if(!is.null(dim(PG_TPMs_BALL))){
    if(nrow(PG_TPMs_BALL) > 10){
    PG_TPMs_BALL_pivot <- PG_TPMs_BALL %>% as.data.frame %>% mutate(Gene = rownames(PG_TPMs_BALL)) %>% 
      pivot_longer(!Gene, names_to = "Samples", values_to = "TPMs") %>% mutate(Phenotype = "BALL")
    
    ensembl_PG <- annot$Ensembl[match(rownames(PG_TPMs_BALL), annot$HGNC_symbol)]
    
    print("Getting sequences")
    sequences_PG = getSequence(id=ensembl_PG, type="ensembl_gene_id", seqType="gene_exon_intron", mart = mart)
    
    sequences_PG <- sequences_PG %>% mutate(Gene_name = annot$HGNC_symbol[match(sequences_PG$ensembl_gene_id, annot$Ensembl)])
    
    print("Getting SeqSims")
    for(j in 1:nrow(sequences_PG)){
      
    if(sequences_PG[j,2] != PG[i]) {
    #if(sequences_PG[j,2] != "TUBA1A") {
      
      sim <- smith_waterman(sequences_PG[sequences_PG$Gene_name == PG[i],1], 
                            sequences_PG[j,1] 
                            , type="characters")
      
      SeqSim.DF <- data.frame(Gene1 = sequences_PG[sequences_PG$Gene_name == PG[i],3],
                              Gene2 = sequences_PG[j,3],
                              SeqSim = sim[["similarity"]])
      
      if(!exists("x")){
        x = SeqSim.DF
      } else {
        x <- rbind(SeqSim.DF, x)
      }
      }
    }
    
    print("Make scatter")
    PG_means_BALL <- PG_TPMs_BALL %>% as.data.frame() %>% mutate(Gene = rownames(PG_TPMs_BALL)) %>% relocate(Gene) %>% 
      rowwise() %>% group_by(Gene) %>% summarise(Mean_BALL = mean(c_across(where(is.numeric))))
    
    PG_corrs_inputs <- x %>% inner_join(PG_means_BALL, by = c("Gene2" = "Gene"))
    
    PG_corrs_inputs <- PG_corrs_inputs[-which(PG_corrs_inputs$Gene1 == PG_corrs_inputs$Gene2),]
    
    scatter_BALL <- ggplot(PG_corrs_inputs, aes(x = Mean_BALL, y = SeqSim)) + geom_point() + geom_smooth(method="lm") +
      labs(title = paste(PG[i]),
           subtitle = "BALL", 
           x = "Mean expression (TPMs)", y = "Sequence similarity")
    
    ggsave(filename =  paste0(dir.scatters, "/ScatterPlotCounts_",PG[i], "_BALL.png"), plot = scatter_BALL, units="in",
           width=10, height=10, dpi =300, bg = "white")
    
    if(!exists("Full_PG_corrs_inputs")){
      Full_PG_corrs_inputs = PG_corrs_inputs
    } else {
      Full_PG_corrs_inputs <- rbind(Full_PG_corrs_inputs, PG_corrs_inputs)
    }
    
    print("Making correlations")
    Corr_BALL <- rcorr(PG_corrs_inputs$Mean_BALL, PG_corrs_inputs$SeqSim, type = "spearman")
    
    Mi_BALL <- mutinformation(discretize(PG_corrs_inputs$Mean_BALL), discretize(PG_corrs_inputs$SeqSim))
    
    corr_df <- data.frame(Gene = PG[i],
                          Rho = Corr_BALL$r[1,2],
                          P = Corr_BALL$P[1,2],
                          MiV_BALL = Mi_BALL)
    
    if(!exists("Full_corr_df")){
      Full_corr_df = corr_df
    } else {
      Full_corr_df <- rbind(Full_corr_df, corr_df)
    }
  } }}
}

saveRDS(Full_corr_df, "Full_corr_df.RDS")
saveRDS(Full_PG_corrs_inputs, "Full_PG_corrs_inputs.RDS")

#Analyze results
Relevant_corr_BALL <- Full_corr_df[Full_corr_df$P < 0.05 & 
                                     Full_corr_df$Rho >= 0.5,]


  
  
