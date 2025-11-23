# TCGA Dataset Extraction Functions

library(TCGAbiolinks)

download_tcga_data <- function(cancer_type = "BRCA", sample_size = NULL) {
  
  project_name <- paste0("TCGA-", cancer_type)
  cat("Project:", project_name, "\n")
  
  # Full query
  query <- GDCquery(
    project = project_name,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
  # Subset if requested
  if (!is.null(sample_size)) {
    query_results <- getResults(query)
    barcodes_subset <- query_results$cases[1:sample_size]
    query <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal"),
      barcode = barcodes_subset
    )
  }
  
  GDCdownload(query)
  return(query)
}

prepare_tcga_data <- function(query) {
  library(TCGAbiolinks)
  data <- GDCprepare(query)
  return(data)
}