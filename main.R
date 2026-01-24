
#Libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

source("lib/tcga_data_extract.R")
source("lib/expression_analysis.R")
source("lib/methylation_analysis.R")
#source("lib/cnv_analysis.R")
#source("lib/mutation_analysis.R")
source("lib/survival_analysis.R")
source("lib/visualization.R")

config <- yaml::yaml.load_file("config.yaml")

data_category <- config$data_categories

cancer_types <- config$cancer_types
cat("Cancer types to analyze:", length(cancer_types), "\n")

# Pan Cancer
for (cancer in cancer_types){
      for (dc in data_category){
        extract_tcga_data(project_name = paste("TCGA", cancer, sep = "-"),
                          data_category = dc)
      }
}

# BRCA
brca_data <- extract_tcga_data(project_name = "TCGA-BRCA",
                               data_category = data_category$transcriptome)

histone_genes <- unlist(config$genes$tcga_gene_names)
cat("Genes to analyze:", paste(histone_genes, collapse = ", "), "\n")

