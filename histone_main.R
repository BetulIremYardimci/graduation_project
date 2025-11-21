
#import libraries
library(TCGAbiolinks)
library(DESeq2)
library(survival)
library(survminer)

# LOAD TCGA DATASET
query_small <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  barcode = results$cases[1:100]
)

GDCdownload(query_small)

results <- getResults(query_small)

#ALL DATASETS
#query_full <- GDCquery(
#  project = "TCGA-BRCA",
#  data.category = "Transcriptome Profiling",
#  data.type = "Gene Expression Quantification",
#  workflow.type = "STAR - Counts",
#  sample.type = c("Primary Tumor", "Solid Tissue Normal")
#)

#GDCdownload(query_full)


# DATA PREPERATION

genes_tcga <- c("H2AX", "MACROH2A1", "MACROH2A2", "H2AZ1")
gene_indices <- which(rowData(data)$gene_name %in% genes_tcga)
gene_indices

gene_data <- data[gene_indices, ]
dim(gene_data)
rowData(gene_data)$gene_name

head(assay(gene_data, "unstranded"))

#count matrix
counts <- assay(gene_data, "unstranded")

#get sample infos
colData <- colData(gene_data)
table(colData$sample_type)

#DESeq2 Analysis
dds <- DESeqDataSet(gene_data, design = ~ sample_type)
dds <- DESeq(dds)
results_deseq <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
results_deseq
