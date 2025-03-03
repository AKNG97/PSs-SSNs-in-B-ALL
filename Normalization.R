#### 1 Load dependencies and functions ####

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

# 1.1 Load functions

Get_raw_matrix <- function(x,y) {
  z <- cbind(x, y)
  print(dim(z))
  print(head(rownames(z)))
  # rownames(z) <- rowData(x)$gene_name
  # print(head(rownames(z)))
  z <- z[rownames(z) %in% annot_TCGA$gene_id,]
  rownames(z) <- annot_TCGA$HGNC_symbol[match(rownames(z), annot_TCGA$gene_id)]
  print(dim(z))
  print(head(rownames(z)))
  print("Filter matrix")
  dataFilt <- TCGAanalyze_Filtering(tabDF = z,
                                    method = "quantile",
                                    qnt.cut = 0.25)
  threshold <- round(dim(z)[2]/2)
  ridx <- rowSums(dataFilt == 0) <= threshold
  dataFilt <- dataFilt[ridx, ]
  ridx <- rowMeans(dataFilt) >= 10
  dataFilt <- dataFilt[ridx, ]
  z <- z[rownames(z) %in% rownames(dataFilt), ]
  print(dim(z))
  return(z)
}

Get_factors_objects <- function(x,y){
  factors <- x[,c("patient", "primary_diagnosis")] %>%
    bind_rows(y[,c("patient", "primary_diagnosis")])
  rownames(factors) <- factors$barcode
  return(factors)
}

norm <- function(x, y, z) {
  y <- y[y$HGNC_symbol %in% rownames(x),]
  y <- y[match(rownames(x), y$HGNC_symbol),]
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$Phenotype))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}




#### 2 Get Biomart annotation  #####

ensembl <- useEnsembl(biomart = "genes", version = 110,
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)
annot <- annot[annot$HGNC_symbol != "",]
annot <- annot[!duplicated(annot$ensembl_gene_id),]
dim(annot)

#### 3 Get RNA-seq data from TCGA #####
#### 3.1 ALL from TARGET ####

RES_DIR <- getwd()

DATA_DIR <- "/datos/ot/anakamura/PS_Analaysis"

dir.create(DATA_DIR)

setwd(DATA_DIR)

qry.ALL <- GDCquery(project = "TARGET-ALL-P2",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
GDCdownload(qry.ALL)
ALL <- GDCprepare(qry.ALL, summarizedExperiment = TRUE)

#### 3.2 ALL from MP2PRT #### 
qry.ALL_MP2 <- GDCquery(project = "MP2PRT-ALL",
                        data.category= "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")
GDCdownload(qry.ALL_MP2)
ALL_MP2 <- GDCprepare(qry.ALL_MP2, summarizedExperiment = TRUE)

#### 3.3 AML and Normal bone marrow #### 
qry.AML_NBM <- GDCquery(project = "TARGET-AML",
                        data.category= "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = c("Primary Blood Derived Cancer - Bone Marrow",
                                        "Bone Marrow Normal"))
# Warning: There are more than one file for the same case. Please verify query results.
# You can use the command View(getResults(query)) in rstudio
# Remove duplicated cases
qry.results.AML <- getResults(qry.AML_NBM)
qry.results.AML <- qry.results.AML[!duplicated(qry.results.AML$cases),]
qry.results.AML <- qry.results.AML[!duplicated(qry.results.AML$cases.submitter_id),]
# Warning: There are more than one file for the same case. Please verify query results.
# You can use the command View(getResults(query)) in rstudio

#Remove unusual names, then download files and get SummarizedExperiment object
cases <- qry.results.AML$cases
dash_count <- str_count(cases, "-")
qry.results.AML <- qry.results.AML[!(dash_count > 4),]
qry.AML_NBM[[1,1]] <- qry.results.AML

GDCdownload(qry.AML_NBM)
AML_NBM <- GDCprepare(qry.AML_NBM, summarizedExperiment = TRUE)

#### 3.4 Save raw data ####

setwd(RES_DIR)

dir.create("RawData")
saveRDS(ALL, "RawData/ALL_TARGET.RDS")
saveRDS(ALL_MP2, "RawData/ALL_MP2PRT.RDS")
saveRDS(AML_NBM, "RawData/AML_NBM_TARGET.RDS")

ALL <- readRDS("RawData/ALL_TARGET.RDS")
AML_NBM <- readRDS("RawData/AML_NBM_TARGET.RDS")

#### 4 Merge annotations from Biomart and TCGA, retain only conserved genes based on Ensembl Version and get PS annotation #####
annot_TCGA <- as.data.frame(rowData(ALL))
annot_TCGA$Ensembl <- gsub("\\.[0-9]*", "", annot_TCGA$gene_id)
dim(annot_TCGA)
annot_TCGA <- annot_TCGA %>% inner_join(annot[,c("ensembl_gene_id","Ensembl_ID_Version","HGNC_symbol", "HGNC_ID", "Type", "Chr", "GC", "Start", "End","Length")], 
                                        by = c("Ensembl" = "ensembl_gene_id"))
dim(annot)
#[1] 41110    10
dim(annot_TCGA)
#[1] 40840    20

which(duplicated(annot_TCGA$gene_name))
#[1] 35361 38048 38537 38751 39864 39865 39911 39940 39975 39980 40702 40785

which(duplicated(annot_TCGA$HGNC_symbol))
# [1] 35325 35361 35513 36343 37955 38537 38681 39864 39865 39911 39940 39975
# [13] 39980 40086 40678 40702 40785 40812

duplicated_GeneNames <- annot_TCGA$HGNC_symbol[which(duplicated(annot_TCGA$HGNC_symbol))]
dim(annot_TCGA)
#40704    19
annot_TCGA <- annot_TCGA[!(annot_TCGA$HGNC_symbol %in% duplicated_GeneNames),]

dim(annot_TCGA)
#[1] 40840    20

#Get Pseudogenes conserved on both Biomart and TCGA annotations

length(grep(".*pseudogene", annot_TCGA$gene_type))
#[1] 10363
length(grep(".*pseudogene", annot_TCGA$Type))
#[1] 10317

annot_PS <- annot_TCGA[intersect(grep(".*pseudogene", annot_TCGA$gene_type), 
                                 grep(".*pseudogene", annot_TCGA$Type)),]
dim(annot_PS)
#[1] 10280    20

write.table(annot_TCGA, "gene_annotation_TCGA.tsv", quote = F, row.names = F, sep = "\t")

#### 5 Filter BALL and NBM (Normal Bone Marrow) samples, retain cases with full OS data and check for replicates #####
dir.create("OS_files")
#### 5.1 BALL from TARGET-ALL2 ####

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
CI_1_ALL <- readxl::read_xlsx("Scripts/ClinicalFiles/ALL_CI_1.xlsx")
CI_2_ALL <- readxl::read_xlsx("Scripts/ClinicalFiles/ALL_CI_2.xlsx")
CI_3_ALL <- readxl::read_xlsx("Scripts/ClinicalFiles/ALL_CI_3.xlsx")
CI_4_ALL <- readxl::read_xlsx("Scripts/ClinicalFiles/ALL_CI_4.xlsx")
CI_5_ALL <- readxl::read_xlsx("Scripts/ClinicalFiles/ALL_CI_5.xlsx")

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
                                             EFS_status = `First Event`,
                                             EFS_time = `Event Free Survival Time in Days`,
                                             OS_time_update = `Overall Survival Time in Days`,
                                             patient = `TARGET USI`)
CI_full_ALL <- CI_full_ALL[!duplicated(CI_full_ALL$patient),]
CI_full_BALL <- CI_full_ALL[CI_full_ALL$`Cell of Origin` %in% c("B Cell ALL", "B-Precursor"),]

CI_BALL <- CI.ALL %>% inner_join(CI_full_BALL[,c("patient", "VS_update", "OS_time_update", "EFS_status", "EFS_time")],
                                 by = c("patient" = "patient"))
table(CI_BALL$VS_update)
# Alive    Dead Unknown 
# 74      65       3 

which(CI_BALL$VS_update == "Unknown")
#[1]   8  24 102

CI_BALL <- CI_BALL[,c("patient", "barcode", "vital_status", "EFS_status", "EFS_time",
                      "days_to_death", "days_to_last_follow_up", "VS_update", "OS_time_update")]
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

dds <- DESeqDataSetFromMatrix(countData = assay(BALL),
                              colData = CI_BALL,
                              design = ~ patient)

ddsColl <- collapseReplicates(dds, dds$patient, dds$patient)

#Check if counts were added
CI_BALL$barcode[126]
#[1] "TARGET-10-PARBVI-09A-01R"
#Collapsed data
head(counts(ddsColl[,"TARGET-10-PARBVI"]))
#                     TARGET-10-PARBVI
# ENSG00000000003.15               21
# ENSG00000000005.6                 0
# ENSG00000000419.13             2866
# ENSG00000000457.14              317
# ENSG00000000460.17              165
# ENSG00000000938.13              141
#Replicate 2
head(counts(dds[,"TARGET-10-PARBVI-09A-02R"]))
#                     TARGET-10-PARBVI-09A-02R
# ENSG00000000003.15                       17
# ENSG00000000005.6                         0
# ENSG00000000419.13                       34
# ENSG00000000457.14                       23
# ENSG00000000460.17                        1
# ENSG00000000938.13                       13
#Replicate 1
head(counts(dds[,"TARGET-10-PARBVI-09A-01R"]))
#                     TARGET-10-PARBVI-09A-01R
# ENSG00000000003.15                        4
# ENSG00000000005.6                         0
# ENSG00000000419.13                     2832
# ENSG00000000457.14                      294
# ENSG00000000460.17                      164
# ENSG00000000938.13                      128

#Data collapsed correctly

#Arrange samples metadata

BALL_TARGET_raw <- counts(ddsColl)
dim(BALL_TARGET_raw)
#[1] 60660   132
dim(CI_BALL)
#[1] 141   8
CI_BALL <- CI_BALL[!duplicated(CI_BALL$patient),]
dim(CI_BALL)
#[1] 132   8

length(which((1:ncol(BALL_TARGET_raw) == match(colnames(BALL_TARGET_raw), CI_BALL$patient)) == FALSE))
#[1] 132 Let's arrange in correct oder

CI_BALL <- CI_BALL[match(colnames(BALL_TARGET_raw), CI_BALL$patient),]

length(which((1:ncol(BALL_TARGET_raw) == match(colnames(BALL_TARGET_raw), CI_BALL$patient)) == FALSE))
#[1] 0

CI_BALL <- CI_BALL[,c("patient", "VS_update", "OS_time_update", "primary_diagnosis")]

colnames(CI_BALL) <- c("Sample", "VitalStatus", "OS_time", "Phenotype")

write.table(CI_BALL, "OS_files/CI_BALL.tsv", quote = F, row.names = F, sep = "\t")
saveRDS(BALL_TARGET_raw, "RawData/BALL_TARGET_raw.RDS")

#### 5.2 NBM from TARGET-AML ####
NBM <- AML_NBM[,AML_NBM$sample_type == "Bone Marrow Normal"]
table(as.factor(NBM$primary_diagnosis))
# Acute myeloid leukemia, NOS 
# 138 
#Hence, there are "normal" samples from AML patients
which(duplicated(NBM$patient))
#integer(0)
table(as.factor(NBM$vital_status))
# Alive  Dead 
# 90    48 
#From these "recovered patients" 55 died later.
#¿Should I manage these samples separately from samples without a primary diagnosis?
#Let's subset them, check on that later

NBM <- NBM[, is.na(NBM$primary_diagnosis)]
dim(NBM)
#[1] 60660    70
NBM_raw <- assay(NBM)
CI_NBM <- as.data.frame(colData(NBM)) %>% mutate(primary_diagnosis = "Normal_Bone_Marrow")

length(which((1:ncol(NBM_raw) == match(colnames(NBM_raw), CI_NBM$barcode)) == FALSE))
#[1] 0
colnames(NBM_raw) <- CI_NBM$patient

CI_NBM <- CI_NBM[,c("patient", "primary_diagnosis")]

colnames(CI_NBM) <- c("Sample", "Phenotype")

write.table(CI_NBM, "OS_files/CI_NBM.tsv", quote = F, row.names = F, sep = "\t")
saveRDS(NBM_raw, "RawData/NBM_raw.RDS")

#### 5.3 BALL from MP2PRT-ALL ####

CI.MP2PRT <- as.data.frame(colData(ALL_MP2))

#Data lacks full days_to_death data

length(which(duplicated(CI.MP2PRT$sample_submitter_id)))
#[1] 0

length(which(duplicated(CI.MP2PRT$submitter_id)))
#[1] 0

#No duplicated patients!!!

Ext_CI <- read.delim("ClinicalFiles/MP2PRT_Outcomes_1510.tsv")

CI.MP2PRT <- CI.MP2PRT %>% inner_join(Ext_CI, by = c("submitter_id" = "cases.submitter_id"))
CI.MP2PRT <- CI.MP2PRT[,c("sample_submitter_id", "sample_type", "submitter_id",
                          "primary_diagnosis", "ethnicity", "race",
                          "EFS_Time","EFS_Status",
                          "DFS_Time","DFS_Status",
                          "OS_Time", "OS_Status")]

CI.MP2PRT <- CI.MP2PRT[CI.MP2PRT$sample_type == "Blood Derived Cancer - Bone Marrow",]
# dim(CI.MP2PRT)
# [1] 1290   12

table(CI.MP2PRT$primary_diagnosis)
# Acute leukemia, NOS 
# 10 
# Acute lymphoblastic leukemia, NOS 
# 121 
# Blastic plasmacytoid dendritic cell neoplasm 
# 1 
# Chronic lymphocytic leukemia 
# 1 
# Common precursor B ALL 
# 1082 
# Leukemia, NOS 
# 48 
# Lymphoid leukemia, NOS 
# 1 
# Neoplasm, uncertain whether benign or malignant 
# 1 
# Osteosarcoma, NOS 
# 1 
# Precursor B-cell lymphoblastic leukemia 
# 9 
# Precursor cell lymphoblastic leukemia, NOS 
# 6 
# Precursor cell lymphoblastic leukemia, not phenotyped 
# 7 
# Subacute lymphoid leukemia 
# 1 
# Unknown 
# 1 

CI.MP2PRT <- CI.MP2PRT[CI.MP2PRT$primary_diagnosis %in% 
                         c("Common precursor B ALL", "Precursor B-cell lymphoblastic leukemia"),]

BALL_MP2_raw <- assay(ALL_MP2[,colnames(ALL_MP2) %in% CI.MP2PRT$sample_submitter_id])
CI.BALL_MP2 <- CI.MP2PRT[match(colnames(BALL_MP2_raw), CI.MP2PRT$sample_submitter_id),]

CI.BALL_MP2 <- CI.BALL_MP2[,c("submitter_id", "OS_Status", "OS_Time", "primary_diagnosis",
                              "EFS_Time", "EFS_Status", "DFS_Time", "DFS_Status")]

colnames(CI.BALL_MP2) <- c("Sample", "VitalStatus", "OS_time",
                           "Phenotype", "EFS_Time", "EFS_Status", 
                           "DFS_Time", "DFS_Status")

CI.BALL_MP2$Phenotype <- "BALL"

write.table(CI.BALL_MP2, "OS_files/CI_BALL_MP2.tsv", quote = F, row.names = F, sep = "\t")
colnames(BALL_MP2_raw) <- CI.BALL_MP2$Sample
saveRDS(BALL_MP2_raw, "RawData/BALL_MP2PRT_raw.RDS")

#### Filter expression matrices ####

BALL_MP2PRT_F <- Get_raw_matrix(BALL_MP2_raw, NBM_raw)

BALL_TARGET_F <- Get_raw_matrix(BALL_TARGET_raw, NBM_raw)

CI.BALL_MP2PRT_F <- rbind(CI.BALL_MP2[,c("Sample", "Phenotype")],
                          CI_NBM)
rownames(CI.BALL_MP2PRT_F) <- CI.BALL_MP2PRT_F$Sample

CI.BALL_TARGET_F <- rbind(CI_BALL[,c("Sample", "Phenotype")],
                          CI_NBM)
rownames(CI.BALL_TARGET_F) <- CI.BALL_TARGET_F$Sample

#### PCAs rqw expression ####

dir.create("PCAs")

#MP2PRT
before.pca <- prcomp(t(BALL_MP2PRT_F),center = TRUE,scale. = TRUE)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_MP2PRT_F$Phenotype)
ggsave("PCAs/PCA_raw_BALL_MP2PRT_BM.pdf")

BALL_MP2PRT_raw_PS <- BALL_MP2PRT_F[rownames(BALL_MP2PRT_F) %in% annot_PS$HGNC_symbol,]

before.pca <- prcomp(t(BALL_MP2PRT_raw_PS),center = TRUE,scale. = TRUE)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_MP2PRT_F$Phenotype)
ggsave("PCAs/PCA_raw_BALL_MP2PRT_BM_PS.pdf")

#TARGET
before.pca <- prcomp(t(BALL_TARGET_F),center = TRUE,scale. = TRUE)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_TARGET_F$Phenotype)
ggsave("PCAs/PCA_raw_BALL_TARGET_BM.pdf")

BALL_TARGET_raw_PS <- BALL_TARGET_F[rownames(BALL_TARGET_F) %in% annot_PS$HGNC_symbol,]

before.pca <- prcomp(t(BALL_TARGET_raw_PS),center = TRUE,scale. = TRUE)
ggbiplot(before.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_TARGET_F$Phenotype)
ggsave("PCAs/PCA_raw_BALL_TARGET_BM_PS.pdf")

#### Normalization ####
dir.create("NormData")

BALL_MP2PRT_norm <- norm(BALL_MP2PRT_F, annot_TCGA, CI.BALL_MP2PRT_F)
saveRDS(BALL_MP2PRT_norm, "NormData/BALL_MP2PRT_norm.RDS")

BALL_TARGET_norm <- norm(BALL_TARGET_F, annot_TCGA, CI.BALL_TARGET_F)
saveRDS(BALL_TARGET_norm, "NormData/BALL_TARGET_norm.RDS")

annot_TCGA<- annot_TCGA[annot_TCGA$HGNC_symbol %in% rownames(BALL_TARGET_norm), ]
annot_TCGA<- annot_TCGA[match(rownames(BALL_TARGET_norm),annot_TCGA$HGNC_symbol),]

saveRDS(annot_TCGA, "annot_TCGA.RDS")

#### PCAs normalized expression ####

#MP2PRT
after.pca <- prcomp(t(BALL_MP2PRT_norm),center = TRUE,scale. = TRUE)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_MP2PRT_F$Phenotype)
ggsave("PCAs/PCA_norm_BALL_MP2PRT_BM.pdf")

BALL_MP2PRT_norm_PS <- BALL_MP2PRT_norm[rownames(BALL_MP2PRT_norm) %in% annot_PS$HGNC_symbol,]
saveRDS(BALL_MP2PRT_norm_PS, "NormData/BALL_MP2PRT_norm_PS_2.RDS")

after.pca <- prcomp(t(BALL_MP2PRT_norm_PS),center = TRUE,scale. = TRUE)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_MP2PRT_F$Phenotype)
ggsave("PCAs/PCA_norm_BALL_MP2PRT_BM_PS.pdf")

#TARGET
after.pca <- prcomp(t(BALL_TARGET_norm),center = TRUE,scale. = TRUE)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_TARGET_F$Phenotype)
ggsave("PCAs/PCA_norm_BALL_TARGET_BM.pdf")

BALL_TARGET_norm_PS <- BALL_TARGET_norm[rownames(BALL_TARGET_norm) %in% annot_PS$HGNC_symbol,]
saveRDS(BALL_TARGET_norm_PS, "NormData/BALL_TARGET_norm_PS.RDS")

after.pca <- prcomp(t(BALL_TARGET_norm_PS),center = TRUE,scale. = TRUE)
ggbiplot(after.pca, var.axes=FALSE, ellipse=TRUE, groups=CI.BALL_TARGET_F$Phenotype)
ggsave("PCAs/PCA_norm_BALL_TARGET_BM_PS.pdf")
