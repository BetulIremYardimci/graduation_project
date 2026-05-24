# Queries the TCGA data depending on project name and Downloads

library(TCGAbiolinks)

library(TCGAbiolinks)

download_tcga_data <- function(project_name, data_category, force = FALSE) {

  safe_category <- gsub(" ", "_", data_category)
  query_file <- paste0("data/queries/", project_name, "_",
                       safe_category, "_query.RData")

  flag_file <- paste0("data/downloaded/", project_name, "_",
                      safe_category, ".flag")

  if (file.exists(flag_file) && !force) {
    message("Already downloaded: ", project_name, " - ", data_category)
    if (file.exists(query_file)) {
      load(query_file)
      return(query)
    }
  }

  message("Downloading: ", project_name, " - ", data_category)

  # Create query based on data category
  if (data_category == "Transcriptome Profiling") {
    query <- GDCquery(
      project = project_name,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
  } else if (data_category == "DNA Methylation") {
    query <- GDCquery(
      project = project_name,
      data.category = "DNA Methylation",
      platform = "Illumina Human Methylation 450",
      data.type = "Methylation Beta Value"
    )
  } else if (data_category == "Simple Nucleotide Variation") {
    query <- GDCquery(
      project = project_name,
      data.category = "Simple Nucleotide Variation",
      access = "open",
      data.type = "Masked Somatic Mutation"
    )
  } else if (data_category == "Copy Number Variation") {
    query <- GDCquery(
      project = project_name,
      data.category = "Copy Number Variation",
      data.type = "Gene Level Copy Number",
      access = "open"
    )
  } else {
    stop("Invalid data_category: ", data_category)
  }

  GDCdownload(query, method = "api", files.per.chunk = 10)

  dir.create("data/queries", showWarnings = FALSE, recursive = TRUE)
  save(query, file = query_file)

  dir.create("data/downloaded", showWarnings = FALSE, recursive = TRUE)
  writeLines(as.character(Sys.time()), flag_file)

  message("Download complete: ", project_name, " - ", data_category)

  return(query)
}

# GDCprepare
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


# Clinical Data - separately fetching clinical data
extract_clinical_data <- function(project_name, force = FALSE) {

  prepared_file <- paste0("data/processed/", project_name, "_clinical.RData")

  if (file.exists(prepared_file) && !force) {
    load(prepared_file)
    return(clinical_data)
  }

  message("Fetching clinical data: ", project_name)

  clinical_data <- GDCquery_clinic(project = project_name, type = "clinical")

  dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
  save(clinical_data, file = prepared_file)

  return(clinical_data)
}



