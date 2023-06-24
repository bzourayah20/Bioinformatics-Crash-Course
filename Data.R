library(tidyverse)
library(readxl)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(ggpubr)


clinical <- read_delim("skcm_dfci_2015.tar/skcm_dfci_2015/data_clinical_patient.txt")

expr <- read_delim("skcm_dfci_2015.tar/skcm_dfci_2015/data_RNA_Seq_expression_median.txt")

expr$Entrez_Gene_Id <- as.character(expr$Entrez_Gene_Id)


expr <- expr %>%
  remove_rownames() %>%
  column_to_rownames(var = colnames(expr)[1])

# expr <- expr[!duplicated(expr$Entrez_Gene_Id),]

expr <- aggregate(expr[,2:41], by=list(name = expr$Entrez_Gene_Id), FUN = mean) 



clinical <- clinical[-c(1:4),]


clinical <- clinical %>%
  remove_rownames() %>%
  column_to_rownames(var = colnames(clinical)[1])

clinical <- clinical %>%
  filter(clinical$`Durable Clinical Benefit`!="X")

clinical$Response <- ifelse(clinical$`Durable Clinical Benefit` == "SD" | clinical$`Durable Clinical Benefit`=="PD", 
                            "NR", "R")
clinical <- clinical[rownames(clinical) %in% colnames(expr),]
expr  <- expr[,colnames(expr) %in% rownames(clinical)]


write.csv(expr, "expressionClean.csv")
