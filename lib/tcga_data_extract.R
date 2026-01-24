# Queries the TCGA data depending on project name and Downloads

library(TCGAbiolinks)

extract_tcga_data <- function(project_name, data_category, force = FALSE) {

  prepared_file <- paste0("data/processed/", project_name, "_",
                          gsub(" ", "_", data_category), "_prepared.RData")

  if (file.exists(prepared_file) && !force) {
    load(prepared_file)
    return(data_result)
  }

  # Transcriptome
  if (data_category == "Transcriptome Profiling") {
    query <- GDCquery(project = project_name,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

    GDCdownload(query, method = "api", files.per.chunk = 10)

    expr_matrix <- GDCprepare_read(query)
    #data_result <- GDCprepare(query) # Raw

    #save(data_result, file = prepared_file)
    return(data_result)
  }

  # Methylation
  if (data_category == "DNA Methylation") {
    query <- GDCquery(project = project_name,
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.type = "Methylation Beta Value" )

    GDCdownload(query, method = "api", files.per.chunk = 10)
    data_result <- GDCprepare(query)

    #save(data_result, file = prepared_file)
    return(data_result)
  }

  # Mutation
  if (data_category == "Simple Nucleotide Variation") {
    query <- GDCquery(project = project_name,
                      data.category = "Simple Nucleotide Variation",
                      access = "open",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

    GDCdownload(query, method = "api", files.per.chunk = 10)
    data_result <- GDCprepare(query)

    #save(data_result, file = prepared_file)
    return(data_result)
  }
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



