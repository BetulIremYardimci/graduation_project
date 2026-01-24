
download_tcga_data <- function(cancer_type = "BRCA", sample_size = NULL) {
  
  library(TCGAbiolinks)
  project_name <- paste0("TCGA-", cancer_type)
  
  query <- GDCquery(
    project = project_name,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
  query_results <- getResults(query)
  
  cat("  Total samples available:", nrow(query_results), "\n")
  cat("  Sample types:\n")
  print(table(query_results$sample_type))
  
  # Subset if requested (for testing)
  if (!is.null(sample_size)) {
    cat("  Subsetting to", sample_size, "samples for testing...\n")
    
    # Get barcode subset
    barcodes_subset <- query_results$cases[1:min(sample_size, nrow(query_results))]
    
    # Create subset query
    query <- GDCquery(
      project = project_name,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal"),
      barcode = barcodes_subset
    )
  }
  
  cat("  Downloading data from GDC...\n")
  cat("  (This may take several minutes depending on dataset size)\n")
  
  GDCdownload(query, method = "api", files.per.chunk = 10)
  
  cat("  Download complete!\n")
  
  return(query)
}


#' Prepare TCGA data for analysis
prepare_tcga_data <- function(query) {
  
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  
  cat("  Preparing data for analysis...\n")
  
  data <- GDCprepare(query, summarizedExperiment = TRUE)
  
  cat("  Data prepared:\n")
  cat("    Samples:", ncol(data), "\n")
  cat("    Genes:", nrow(data), "\n")
  cat("    Assays:", paste(assayNames(data), collapse = ", "), "\n")
  
  cat("  Sample type distribution:\n")
  print(table(colData(data)$sample_type))
  
  return(data)
}


#' Download multiple cancer types for pan-cancer analysis

download_multiple_cancers <- function(cancer_types, sample_size = NULL) {
  
  queries <- list()
  
  cat("Downloading", length(cancer_types), "cancer types\n")

  for (cancer in cancer_types) {
    cat("\n--- Processing:", cancer, "---\n")
    
    tryCatch({
      query <- download_tcga_data(
        cancer_type = cancer,
        sample_size = sample_size
      )
      queries[[cancer]] <- query
      
      # Small delay to avoid overwhelming GDC server
      Sys.sleep(2)
      
    }, error = function(e) {
      cat("✗ Failed to download", cancer, ":", e$message, "\n")
      queries[[cancer]] <- NULL
    })
  }
  
  return(queries)
}


#' Get sample size information for cancer type

get_cancer_sample_info <- function(cancer_type) {
  
  library(TCGAbiolinks)
  
  project_name <- paste0("TCGA-", cancer_type)
  
  # Query without downloading
  query <- GDCquery(
    project = project_name,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  results <- getResults(query)
  
  info <- data.frame(
    cancer_type = cancer_type,
    total_samples = nrow(results),
    primary_tumor = sum(results$sample_type == "Primary Tumor"),
    normal_tissue = sum(results$sample_type == "Solid Tissue Normal"),
    metastatic = sum(results$sample_type == "Metastatic"),
    other = nrow(results) - sum(results$sample_type %in% 
                                  c("Primary Tumor", "Solid Tissue Normal", "Metastatic"))
  )
  
  return(info)
}


#' Get sample information for all cancer types

get_all_cancer_info <- function(cancer_types) {
  
  info_list <- lapply(cancer_types, function(cancer) {
    cat("Querying", cancer, "...\n")
    tryCatch({
      get_cancer_sample_info(cancer)
    }, error = function(e) {
      data.frame(
        cancer_type = cancer,
        total_samples = NA,
        primary_tumor = NA,
        normal_tissue = NA,
        metastatic = NA,
        other = NA
      )
    })
  })
  
  info_table <- do.call(rbind, info_list)
  
  return(info_table)
}