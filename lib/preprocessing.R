
#Batch Processing

library(SummarizedExperiment)
library(TCGAbiolinks)


process_in_batches_with_query <- function(query, batch_size = 100) {

  all_samples <- getResults(query)
  n_samples <- nrow(all_samples)

  message("Total samples: ", n_samples)
  message("Processing in batches of ", batch_size)

  all_batches <- list()

  for (i in seq(1, n_samples, by = batch_size)) {

    end_idx <- min(i + batch_size - 1, n_samples)
    batch_samples <- all_samples$cases[i:end_idx]

    message("Batch: samples ", i, " to ", end_idx)

    batch_query <- GDCquery(
      project = query$project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      barcode = batch_samples
    )

    batch_data <- GDCprepare(batch_query)

    counts <- assay(batch_data, "unstranded")
    rownames(counts) <- sub("\\..*$", "", rownames(counts))  # Ensembl .xx temizleme

    all_batches[[length(all_batches) + 1]] <- counts

    rm(batch_data, counts)
    gc()
  }

  final_data <- do.call(cbind, all_batches)
  return(final_data)
}

#################
sample_code_check <- function(counts_matrix, min_samples = 50,
                              keep_codes = c("01", "11")) {

  library(dplyr)

  message("\n=== Sample Code Analysis ===")

  sample_codes <- substr(colnames(counts_matrix), 14, 15)
  code_table <- table(sample_codes)

  message("\nSample codes found:")
  for (code in names(code_table)) {
    code_name <- case_when(
      code == "01" ~ "Primary Tumor",
      code == "02" ~ "Recurrent Tumor",
      code == "03" ~ "Primary Blood Derived Cancer",
      code == "06" ~ "Metastatic",
      code == "11" ~ "Solid Tissue Normal",
      TRUE ~ "Other"
    )
    message("  Code ", code, " (", code_name, "): ", code_table[code], " samples")
  }

  keep_samples <- sample_codes %in% keep_codes

  removed_codes <- unique(sample_codes[!keep_samples])
  removed_samples <- list()

  for (code in removed_codes) {
    code_samples <- colnames(counts_matrix)[sample_codes == code]
    removed_samples[[code]] <- list(
      code = code,
      count = length(code_samples),
      sample_ids = code_samples,
      reason = ifelse(
        code_table[code] < min_samples,
        paste0("Too few samples (", code_table[code], " < ", min_samples, ")"),
        "Not primary tumor or normal tissue"
      )
    )
  }
  # Filter data
  counts_filtered <- counts_matrix[, keep_samples]

  message("\n=== Filtering Summary ===")
  message("Original samples: ", ncol(counts_matrix))
  message("Kept samples: ", ncol(counts_filtered))
  message("Removed samples: ", sum(!keep_samples))

  if (length(removed_samples) > 0) {
    message("\nRemoved sample types:")
    for (code in names(removed_samples)) {
      info <- removed_samples[[code]]
      message("  Code ", info$code, ": ", info$count, " samples")
      message("    Reason: ", info$reason)
    }
  }

  final_codes <- substr(colnames(counts_filtered), 14, 15)
  message("\nFinal sample distribution:")
  print(table(final_codes))

  return(list(
    data_filtered = counts_filtered,
    removed_samples = removed_samples,
    filtering_report = data.frame(
      original_samples = ncol(counts_matrix),
      kept_samples = ncol(counts_filtered),
      removed_samples = sum(!keep_samples),
      codes_kept = paste(keep_codes, collapse = ", "),
      codes_removed = paste(removed_codes, collapse = ", ")
    )
  ))
}


######################
# Process and create SummarizedExperiment

library(SummarizedExperiment)

convert_to_se <- function(counts_matrix) {

  library(SummarizedExperiment)

  message("\n→ Converting to SummarizedExperiment (full transcriptome)...")

  ## ------------------------------------------------
  ## FIX: collapse duplicated Ensembl IDs
  ## ------------------------------------------------
  if (any(duplicated(rownames(counts_matrix)))) {

    message("  Collapsing duplicated Ensembl IDs by summing counts...")

    counts_matrix <- rowsum(
      counts_matrix,
      group = rownames(counts_matrix)
    )
  }

  ## ------------------------------------------------
  ## Gene metadata
  ## ------------------------------------------------
  gene_metadata <- data.frame(
    ensembl_id = rownames(counts_matrix),
    row.names  = rownames(counts_matrix)
  )

  ## ------------------------------------------------
  ## Sample metadata
  ## ------------------------------------------------
  sample_codes <- substr(colnames(counts_matrix), 14, 15)

  sample_metadata <- data.frame(
    sample_id   = colnames(counts_matrix),
    sample_code = sample_codes,
    sample_type = factor(
      ifelse(sample_codes == "01", "Tumor", "Normal"),
      levels = c("Normal", "Tumor")
    ),
    row.names = colnames(counts_matrix)
  )

  ## ------------------------------------------------
  ## Create SummarizedExperiment
  ## ------------------------------------------------
  se <- SummarizedExperiment(
    assays  = list(counts = counts_matrix),
    rowData = gene_metadata,
    colData = sample_metadata
  )

  ## ------------------------------------------------
  ## Summary
  ## ------------------------------------------------
  message("Genes: ", nrow(se))
  message("Samples: ", ncol(se))
  print(table(colData(se)$sample_type))

  return(se)
}


apply_combat_seq <- function(se, batch_var = "auto", min_samples = 3) {

  library(sva)

  cat("\n>>> Applying ComBat-seq for batch correction...\n")

  # Extract batch info
  if (batch_var == "auto") {
    batch <- substr(colnames(se), 6, 7)
    cat("   Using TSS (Tissue Source Site) as batch\n")
  } else {
    batch <- colData(se)[[batch_var]]
    cat("   Using", batch_var, "as batch\n")
  }

  cat("   Unique batches:", length(unique(batch)), "\n")

  # Filter small batches
  batch_counts <- table(batch)
  small_batches <- names(batch_counts[batch_counts < min_samples])

  if (length(small_batches) > 0) {
    cat("   Removing", length(small_batches),
        "batches with <", min_samples, "samples\n")
    cat("   Small batches:", paste(small_batches, collapse = ", "), "\n")

    keep_samples <- !batch %in% small_batches
    se <- se[, keep_samples]
    batch <- batch[keep_samples]

    cat("   Samples removed:", sum(!keep_samples), "\n")
    cat("   Samples kept:", sum(keep_samples), "\n")
  }

  cat("   Final batches:", length(unique(batch)), "\n")
  print(table(batch))

  # Get sample groups
  group <- colData(se)$sample_type

  # Apply ComBat-seq
  cat("\n   Running ComBat-seq...\n")
  counts_raw <- assay(se, "counts")

  counts_corrected <- ComBat_seq(
    counts = counts_raw,
    batch = batch,
    group = group
  )

  mode(counts_corrected) <- "integer"

  # Update SE
  assay(se, "counts") <- counts_corrected

  cat("✓ ComBat-seq correction complete!\n\n")

  return(se)
}
