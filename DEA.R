library(tidyverse)
library(readxl)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(ggpubr)

all(colnames(expr) == rownames(clinical))

#Create a DESeq2 objects (dds)

dds <- DESeqDataSetFromMatrix(countData = round(expr), 
                              colData = clinical, 
                              design = ~Response)
#Quality Control
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

#Factor level 
dds$Response <- relevel(dds$Response, ref = "NR") 


#Run DESeq2 
dds <- DESeq(dds)
dds.res <- results(dds)

summary(dds.res)
dds.res <- as.data.frame(dds.res)

dds.res.filt <- dds.res %>%
  filter(pvalue <=0.05)


library(EnhancedVolcano)


EnhancedVolcano(dds.res,
                lab = rownames(dds.res),
                x = 'log2FoldChange',
                y = 'pvalue',
                #electLab = rownames(dds.res)[which(names(keyvals) %in% c('Upregulated', 'Downregulated'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 1.5,
                labSize = 2,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 1,
                legendPosition = 'top',
                #legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black', 
                title = "Immunotherapy DEGs",subtitle = "R vs. NR",
                legendLabSize = 11,
                legendLabels=c('NS','Log2FC','pvalue',
                               'pvalue&Log2FC'))

write.csv(dds.res, "DEAresults.csv")
