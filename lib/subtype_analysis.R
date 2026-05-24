#Subtype Analysis
#Updated 21 May 2026

cancer <- "BRCA"
padj_threshold <- config$analysis$padj_threshold
log2fc_threshold <- config$analysis$log2fc_threshold

config <- yaml::yaml.load_file("config.yaml")

source("lib/differential_expression_updated.R")
source("lib/nonparametric_tests.R")

subtype_dir <- paste0("results/", cancer, "/subtype_analysis")
dir.create(subtype_dir, recursive = TRUE, showWarnings = FALSE)

#--------- SUBTYPE DATA EXTRACTION--------------
clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer), type = "clinical")
# to search interested col names:
grep("subtype|pam50|molecular|class", colnames(clinical),
     ignore.case = TRUE, value = TRUE)

#Subtype information get from PanCanAtlas
library(TCGAbiolinks)
brca_subtypes <- PanCancerAtlas_subtypes() %>%
  filter(cancer.type == "BRCA")

colnames(brca_subtypes)
table(brca_subtypes$Subtype_Selected, useNA = "ifany")

#Sample IDs must be checked to assign in SE
head(brca_subtypes$pan.samplesID, 20)

#load SE file
load("data/processed/BRCA_se_combat.RData")
se_samples <- colnames(se)

se_short <- substr(colnames(se), 1, 16)
pancan_short <- substr(brca_subtypes$pan.samplesID, 1, 16)
match_index <- match(se_short, pancan_short)

#subtype annotation
cat("Total SE samples:", ncol(se), "\n")
cat("Matched samples:", sum(!is.na(match_index)), "\n")
cat("Non matched samples:", sum(is.na(match_index)), "\n")

colData(se)$subtype <- brca_subtypes$Subtype_Selected[match_index]
table(colData(se)$sample_type, colData(se)$subtype, useNA = "ifany")
#useNA="ifany", shows NAs in a separate column

#Subtype analysis is proceeded with only tumor samples, bc information for normal
#samples not provided enough

se_tumor <- se[, colData(se)$sample_type == "Tumor" &
                 !is.na(colData(se)$subtype)]

table(colData(se_tumor)$subtype)

round(prop.table(table(colData(se_tumor)$subtype)) * 100, 1)

#---------DATA PREPARATION--------
if (cancer == "BRCA") {
  se_tumor_filtered <- se_tumor[, colData(se_tumor)$subtype != "BRCA.Normal"]
  cat("  Excluded BRCA.Normal subtype\n")
} else {
  se_tumor_filtered <- se_tumor
}

cat("  Final tumor sample count:", ncol(se_tumor_filtered), "\n")
cat("  Subtype distribution:\n")
print(table(colData(se_tumor_filtered)$subtype))
cat("\n")

table(colData(se_tumor_filtered)$subtype)

subtypes_to_analyze <- c("BRCA.LumA", "BRCA.LumB", "BRCA.Basal", "BRCA.Her2") #hardcoded
se_normal <- se[, colData(se)$sample_type == "Normal"]

subtype_counts <- table(colData(se_tumor_filtered)$subtype)
subtype_results <- list()

histone_genes <- unlist(config$genes$tcga_gene_names)
ensembl_ids <- unlist(config$genes$ensembl_ids)
gene_map <- setNames(histone_genes, ensembl_ids)

min_samples_threshold <- 20
#subtypes pass the threshold to analysis
subtypes_to_analyze <- names(subtype_counts[subtype_counts >= min_samples_threshold])

excluded <- names(subtype_counts[subtype_counts < min_samples_threshold])
if (length(excluded) > 0) {
  cat("\n  [INFO] Excluded subtypes (insufficient samples):\n")
  for (st in excluded) {
    cat("    -", st, ":", subtype_counts[st], "samples\n")
  }
}
cat("\n")

#final preprocessing
se_tumor_final <- se_tumor_filtered[,
                                    colData(se_tumor_filtered)$subtype %in% subtypes_to_analyze]

se_normal <- se[, colData(se)$sample_type == "Normal"]


#--------------------EXPRESSION ANALYSIS FOR EACH SUBTYPE-------------
all_results <- list()
all_pvalues <- c()  # For FDR correction

for (subtype_name in subtypes_to_analyze){
  cat("Analyzing:", subtype_name, "\n")

  se_subtype_tumor <- se_tumor_final[,
                                     colData(se_tumor_final)$subtype == subtype_name]

  se_subset <- cbind(se_subtype_tumor, se_normal)

  cat("  Tumor (", subtype_name, "):", ncol(se_subtype_tumor), "\n")
  cat("  Normal:", ncol(se_normal), "\n")

  deseq2_result <- run_deseq2_analysis(
    data = se_subset,
    padj_threshold = padj_threshold,
    log2fc_threshold = log2fc_threshold,
    use_sva = FALSE,
    n_sv = NULL,
    condition = "sample_type",
    ref_group = "Normal",
    case_group = "Tumor"
  )

  deseq2_histone <- deseq2_result$results %>%
    filter(ensembl_id %in% ensembl_ids)

  deseq2_histone$gene_name <- gene_map[deseq2_histone$ensembl_id]
  deseq2_histone$subtype <- subtype_name

  dds <- deseq2_result$dds
  vst_obj <- vst(dds, blind = FALSE)
  vst_counts_subset <- assay(vst_obj)

  wilcox_result <- run_nonparametric_tests(
    expr_matrix = vst_counts_subset,
    histone_df = deseq2_histone,
    cancer = paste0(cancer, "_", subtype_name),
    padj_method = "BH"
  )

  wilcox_result$subtype <- subtype_name

  all_results[[subtype_name]] <- list(
    deseq2 = deseq2_histone,
    wilcox = wilcox_result
  )

  all_pvalues <- c(all_pvalues, deseq2_histone$padj)

  write.csv(deseq2_histone,
            paste0(subtype_dir, "/", subtype_name, "_deseq2_histone.csv"),
            row.names = FALSE)

  write.csv(wilcox_result,
            paste0(subtype_dir, "/", subtype_name, "_wilcoxon_histone.csv"),
            row.names = FALSE)

}

adjusted_pvalues <- p.adjust(all_pvalues, method = "BH")

cat("Total tests:", length(all_pvalues), "\n")
cat("Significant (original padj < 0.05):", sum(all_pvalues < 0.05, na.rm = TRUE), "\n")
cat("Significant (global FDR padj < 0.05):", sum(adjusted_pvalues < 0.05, na.rm = TRUE), "\n\n")

combined_deseq2 <- do.call(rbind, lapply(all_results, function(x) x$deseq2))

combined_wilcox <- do.call(rbind, lapply(all_results, function(x) x$wilcox))

write.csv(combined_deseq2,
          paste0(subtype_dir, "/ALL_subtypes_deseq2_combined.csv"),
          row.names = FALSE)

write.csv(combined_wilcox,
          paste0(subtype_dir, "/ALL_subtypes_wilcoxon_combined.csv"),
          row.names = FALSE)


