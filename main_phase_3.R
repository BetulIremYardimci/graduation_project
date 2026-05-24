
# SURVIVAL ANALYSIS

library(survival)
library(survminer)
library(dplyr)
library(TCGAbiolinks)

source("lib/tcga_data_extract.R")

# 1. Clinical data çek
clinical_luad <- extract_clinical_data("TCGA-LUAD")
#clinical_brca <- extract_clinical_data("TCGA-BRCA")

# Yapıyı kontrol et
dim(clinical_luad)
colnames(clinical_luad)

# main_phase_3.R - Survival Analysis
# =====================================
library(survival)
library(survminer)
library(dplyr)

source("tcga_data_extract.R")

cancer_types <- c("LUSC", "PAAD", "THYM", "UCEC")

results_dir_lusc <- "results/LUSC"
results_dir_paad <- "results/PAAD"

histone_genes <- c("H2AX", "H2AZ1", "H2AFY", "H2AFY2")

# Ensembl ID'ler
histone_ensembl <- c(
  "H2AX"   = "ENSG00000188486",
  "H2AZ1"  = "ENSG00000106153",
  "H2AFY"  = "ENSG00000137309",
  "H2AFY2" = "ENSG00000099284"
)

#------------------------------------------------
# STEP 1: Clinical data hazırla
#------------------------------------------------
prepare_survival_data <- function(cancer) {

  project <- paste0("TCGA-", cancer)
  clinical <- extract_clinical_data(project)

  # OS hesapla: ölüm varsa days_to_death, yoksa days_to_last_follow_up
  clinical$OS_time <- ifelse(
    !is.na(clinical$days_to_death),
    clinical$days_to_death,
    clinical$days_to_last_follow_up
  )

  # Event: 1 = dead, 0 = alive/censored
  clinical$OS_event <- ifelse(
    clinical$vital_status == "Dead", 1, 0
  )

  # Stage basitleştir: I, II, III, IV
  clinical$stage_simple <- gsub("[ABC]$", "", clinical$ajcc_pathologic_stage)
  clinical$stage_simple[clinical$stage_simple %in%
                          c("Stage X", "Not Reported", "")] <- NA

  # Negatif OS_time'ları çıkar
  clinical <- clinical[!is.na(clinical$OS_time) & clinical$OS_time > 0, ]

  cat(cancer, "- Clinical samples:", nrow(clinical), "\n")
  return(clinical)
}

clinical_luad <- prepare_survival_data("LUAD")
clinical_brca <- prepare_survival_data("BRCA")

# Kontrol
table(clinical_luad$OS_event)
table(clinical_luad$stage_simple)

table(clinical_brca$OS_event)
table(clinical_brca$stage_simple)


#------------------------------------------------
# STEP 2: Expression data ile join
#------------------------------------------------
prepare_expression_survival <- function(cancer, clinical) {

  # Normalized counts yükle
  vst_mat <- readRDS(paste0("results/", cancer, "/normalized_counts.rds"))

  # Matrix'e çevir (eğer değilse)
  if (!is.matrix(vst_mat)) vst_mat <- as.matrix(vst_mat)

  # Histone genlerini çek
  histone_rows <- rownames(vst_mat)[rownames(vst_mat) %in% histone_ensembl]
  histone_expr <- as.data.frame(t(vst_mat[histone_rows, ]))

  # Ensembl ID → gene symbol
  colnames(histone_expr) <- names(histone_ensembl)[match(
    colnames(histone_expr), histone_ensembl
  )]

  # Barcode düzenle: TCGA-XX-XXXX formatına indir
  histone_expr$bcr_patient_barcode <- substr(rownames(histone_expr), 1, 12)

  # Sadece tumor sample'ları al (01)
  sample_types <- substr(rownames(histone_expr), 14, 15)
  histone_expr <- histone_expr[sample_types == "01", ]

  # Clinical ile join
  merged <- inner_join(clinical, histone_expr, by = "bcr_patient_barcode")

  cat(cancer, "- Merged samples:", nrow(merged), "\n")

  return(merged)
}

surv_luad <- prepare_expression_survival("LUAD", clinical_luad)
surv_brca <- prepare_expression_survival("BRCA", clinical_brca)


#------------------------------------------------
# STEP 3: Kaplan-Meier + Cox Analysis
#------------------------------------------------
run_survival_analysis <- function(surv_data, cancer, gene) {

  results_dir <- paste0("results/", cancer, "/plots")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  if (!gene %in% colnames(surv_data)) {
    cat("Skipping", gene, "- not found\n")
    return(NULL)
  }

  # High/Low grupları
  median_expr <- median(surv_data[[gene]], na.rm = TRUE)
  surv_data$group <- factor(
    ifelse(surv_data[[gene]] > median_expr, "High", "Low"),
    levels = c("Low", "High")
  )

  km_fit <- survfit(Surv(OS_time, OS_event) ~ group, data = surv_data)
  log_rank <- survdiff(Surv(OS_time, OS_event) ~ group, data = surv_data)
  pval <- round(1 - pchisq(log_rank$chisq, df = 1), 4)

  # KM Plot
  p_km <- ggsurvplot(
    km_fit,
    data = surv_data,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,   # Değişti
    palette = c("#457B9D", "#E63946"),
    title = paste0(gene, " — Overall Survival (", cancer, ")"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    legend.title = "",
    legend.labs = c("Low expression", "High expression"),
    surv.median.line = "hv",
    size = 1.2,
    ggtheme = theme_classic(base_size = 14),
    tables.theme = theme_cleantable()
  )

  # Arka plan siyahlığını düzelt
  p_km$plot <- p_km$plot +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  p_km$table <- p_km$table +
    theme(plot.background = element_rect(fill = "white", color = NA))

  png(paste0(results_dir, "/km_", tolower(gene), ".png"),
      width = 3000, height = 2400, res = 300)
  print(p_km)
  dev.off()

  surv_data[[paste0(gene, "_scaled")]] <- scale(surv_data[[gene]])[, 1]
  scaled_gene <- paste0(gene, "_scaled")

  # Univariate Cox
  cox_formula_uni <- as.formula(
    paste0("Surv(OS_time, OS_event) ~ ", scaled_gene)
  )
  cox_uni <- coxph(cox_formula_uni, data = surv_data)
  cox_summary <- summary(cox_uni)

  # Multivariate Cox
  surv_data_cox <- surv_data[!is.na(surv_data$stage_simple), ]
  surv_data_cox[[paste0(gene, "_scaled")]] <- scale(surv_data_cox[[gene]])[, 1]

  cox_formula_multi <- as.formula(
    paste0("Surv(OS_time, OS_event) ~ ", scaled_gene, " + stage_simple")
  )
  cox_multi <- coxph(cox_formula_multi, data = surv_data_cox)

  result <- data.frame(
    cancer    = cancer,
    gene      = gene,
    n         = nrow(surv_data),
    HR_uni    = round(cox_summary$conf.int[1, 1], 3),
    CI_low    = round(cox_summary$conf.int[1, 3], 3),
    CI_high   = round(cox_summary$conf.int[1, 4], 3),
    p_uni     = round(cox_summary$coefficients[1, 5], 4),
    p_logrank = pval
  )

  cat(cancer, gene, "- HR:", result$HR_uni, "p:", result$p_uni, "\n")

  return(list(result = result, cox_multi = cox_multi))
}

all_results <- list()

for (cancer in c("LUAD", "BRCA")) {
  surv_data <- if (cancer == "LUAD") surv_luad else surv_brca

  for (gene in histone_genes) {
    key <- paste0(cancer, "_", gene)
    all_results[[key]] <- run_survival_analysis(surv_data, cancer, gene)
  }
}

summary_df <- do.call(rbind, lapply(all_results, function(x) x$result))
print(summary_df)

write.csv(summary_df, "results/survival_summary.csv", row.names = FALSE)


#------------------------------------------------
# STEP 4: Forest Plot
#------------------------------------------------
# HR label ekle
summary_df$hr_label <- paste0(
  summary_df$HR_uni,
  " (", summary_df$CI_low, "-", summary_df$CI_high, ")"
)
summary_df$p_label <- ifelse(
  summary_df$p_uni < 0.001, "p<0.001",
  paste0("p=", summary_df$p_uni)
)

ggplot(summary_df, aes(x = HR_uni, y = reorder(gene, HR_uni))) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = significant),
                 height = 0.25, linewidth = 1) +
  geom_point(aes(color = significant), size = 4) +
  # HR değerini sağa yaz
  geom_text(aes(x = max(CI_high) + 0.05,
                label = hr_label),
            hjust = 0, size = 3.2, color = "gray30") +
  # p değerini sola yaz
  geom_text(aes(x = min(CI_low) - 0.05,
                label = p_label),
            hjust = 1, size = 3.2, color = "gray30") +
  scale_color_manual(
    values = c("FALSE" = "gray60", "TRUE" = "#E63946"),
    labels = c("p ≥ 0.05", "p < 0.05")
  ) +
  facet_wrap(~ cancer, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    limits = c(0.55, 1.7),
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4)
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Histone Variant Expression — Overall Survival",
    subtitle = "Hazard Ratio per 1 SD increase (univariate Cox)",
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    color = "Significance"
  ) +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    strip.text   = element_text(face = "bold", size = 12,
                                color = "white"),
    strip.background = element_rect(fill = "#1D3557"),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.spacing  = unit(1, "lines"),
    legend.position = "bottom"
  )

ggsave("results/survival_forest_plot_v2.png",
       width = 12, height = 7, dpi = 300)

