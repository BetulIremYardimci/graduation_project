
prepare_survival_data <- function(data, normalized_counts, genes) {
  
  library(SummarizedExperiment)
  
  cat("\n=== Preparing Survival Data ===\n")
  
  clinical <- as.data.frame(colData(data))
  
  required_cols <- c("days_to_death", "days_to_last_follow_up", "vital_status")
  missing_cols <- required_cols[!required_cols %in% colnames(clinical)]
  
  if (length(missing_cols) > 0) {
    cat("WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }
  
  survival_data <- data.frame(
    sample = rownames(clinical),
    
    os_time = ifelse(
      is.na(clinical$days_to_death),
      clinical$days_to_last_follow_up,
      clinical$days_to_death
    ),
    
    os_status = ifelse(
      clinical$vital_status == "Dead" | clinical$vital_status == "1", 
      1, 0
    ),
    
    os_months = ifelse(
      is.na(clinical$days_to_death),
      clinical$days_to_last_follow_up / 30.44,
      clinical$days_to_death / 30.44
    )
  )
  
  for (gene in genes) {
    gene_idx <- which(rownames(normalized_counts) == gene)
    if (length(gene_idx) > 0) {
      survival_data[[gene]] <- normalized_counts[gene_idx, survival_data$sample]
      
      gene_median <- median(survival_data[[gene]], na.rm = TRUE)
      survival_data[[paste0(gene, "_group")]] <- 
        ifelse(survival_data[[gene]] > gene_median, "High", "Low")
    }
  }
  
  clinical_vars <- c("age_at_diagnosis", "gender", "tumor_stage", 
                     "histological_grade", "er_status", "pr_status")
  
  for (var in clinical_vars) {
    if (var %in% colnames(clinical)) {
      survival_data[[var]] <- clinical[[var]]
    }
  }
  
  complete_idx <- complete.cases(survival_data[, c("os_time", "os_status")])
  survival_data <- survival_data[complete_idx, ]
  
  cat("Total samples with survival data:", nrow(survival_data), "\n")
  cat("Events (deaths):", sum(survival_data$os_status), "\n")
  cat("Censored:", sum(survival_data$os_status == 0), "\n")
  cat("Median follow-up:", 
      round(median(survival_data$os_months, na.rm = TRUE), 1), "months\n")
  
  return(survival_data)
}


run_univariate_cox <- function(survival_data, genes, cancer_type = "") {
  
  library(survival)
  library(dplyr)
  
  results <- data.frame()
  
  for (gene in genes) {
    if (!gene %in% colnames(survival_data)) {
      cat("Skipping", gene, "- not in dataset\n")
      next
    }
    
    cat("Analyzing", gene, "...\n")
    
    cox_formula <- as.formula(paste0("Surv(os_time, os_status) ~ ", gene))
    
    tryCatch({
      cox_model <- coxph(cox_formula, data = survival_data)
      cox_summary <- summary(cox_model)
      
      results <- rbind(
        results,
        data.frame(
          cancer_type = cancer_type,
          gene = gene,
          n_samples = cox_summary$n,
          n_events = cox_summary$nevent,
          HR = cox_summary$conf.int[1, 1],
          HR_lower = cox_summary$conf.int[1, 3],
          HR_upper = cox_summary$conf.int[1, 4],
          coef = cox_summary$coefficients[1, 1],
          se_coef = cox_summary$coefficients[1, 3],
          z_score = cox_summary$coefficients[1, 4],
          pvalue = cox_summary$coefficients[1, 5],
          concordance = cox_summary$concordance[1]
        )
      )
      
      cat(sprintf("  HR = %.2f (%.2f-%.2f), p = %.4f\n",
                  cox_summary$conf.int[1, 1],
                  cox_summary$conf.int[1, 3],
                  cox_summary$conf.int[1, 4],
                  cox_summary$coefficients[1, 5]))
      
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    })
  }
  
  if (nrow(results) > 0) {
    results$significant <- results$pvalue < 0.05
    
    cat("Genes analyzed:", nrow(results), "\n")
    cat("Significant (p < 0.05):", sum(results$significant), "\n")
  }
  
  return(results)
}


run_multivariate_cox <- function(survival_data, 
                                 genes, 
                                 covariates = c("age_at_diagnosis", "tumor_stage"),
                                 cancer_type = "") {
  
  library(survival)
  
  cat("\n=== Multivariate Cox Regression ===\n")
  cat("Covariates:", paste(covariates, collapse = ", "), "\n\n")
  
  available_covariates <- covariates[covariates %in% colnames(survival_data)]
  
  if (length(available_covariates) == 0) {
    cat("WARNING: No covariates available for adjustment\n")
    return(NULL)
  }
  
  cat("Using covariates:", paste(available_covariates, collapse = ", "), "\n")
  
  results <- data.frame()
  
  for (gene in genes) {
    if (!gene %in% colnames(survival_data)) {
      next
    }
    
    cat("\nAnalyzing", gene, "...\n")
    
    formula_str <- paste0("Surv(os_time, os_status) ~ ", gene, " + ",
                          paste(available_covariates, collapse = " + "))
    
    cox_formula <- as.formula(formula_str)
    
    tryCatch({
      cox_model <- coxph(cox_formula, data = survival_data)
      cox_summary <- summary(cox_model)
      
      gene_idx <- 1
      
      results <- rbind(
        results,
        data.frame(
          cancer_type = cancer_type,
          gene = gene,
          n_samples = cox_summary$n,
          n_events = cox_summary$nevent,
          HR_adjusted = cox_summary$conf.int[gene_idx, 1],
          HR_lower = cox_summary$conf.int[gene_idx, 3],
          HR_upper = cox_summary$conf.int[gene_idx, 4],
          pvalue_adjusted = cox_summary$coefficients[gene_idx, 5],
          covariates = paste(available_covariates, collapse = ";")
        )
      )
      
      cat(sprintf("  Adjusted HR = %.2f (%.2f-%.2f), p = %.4f\n",
                  cox_summary$conf.int[gene_idx, 1],
                  cox_summary$conf.int[gene_idx, 3],
                  cox_summary$conf.int[gene_idx, 4],
                  cox_summary$coefficients[gene_idx, 5]))
      
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    })
  }
  
  return(results)
}


#Kaplan-Meier survival curves

create_km_curve <- function(survival_data, 
                            gene, 
                            output_file = NULL,
                            cancer_type = "") {
  
  library(survival)
  library(survminer)
  
  gene_group <- paste0(gene, "_group")
  
  if (!gene_group %in% colnames(survival_data)) {
    cat("WARNING: Gene group not found:", gene_group, "\n")
    return(NULL)
  }
  
  km_formula <- as.formula(paste0("Surv(os_time, os_status) ~ ", gene_group))
  km_fit <- survfit(km_formula, data = survival_data)
  
  logrank <- survdiff(km_formula, data = survival_data)
  pval <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
  
  p <- ggsurvplot(
    km_fit,
    data = survival_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = paste0(cancer_type, " - ", gene, " Expression"),
    legend.title = gene,
    legend.labs = c("High", "Low"),
    palette = c("#E41A1C", "#377EB8"),
    ggtheme = theme_bw()
  )
  
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = print(p),
      width = 8,
      height = 8,
      dpi = 300
    )
    cat("Saved KM curve:", output_file, "\n")
  }
  
  return(list(
    fit = km_fit,
    plot = p,
    pvalue = pval
  ))
}


find_optimal_cutpoint <- function(survival_data, gene) {
  
  library(survminer)
  library(survival)
  
  if (!gene %in% colnames(survival_data)) {
    return(NULL)
  }
  

  optimal <- surv_cutpoint(
    survival_data,
    time = "os_time",
    event = "os_status",
    variables = gene,
    minprop = 0.1  # At least 10% in each group
  )

  
  survival_data[[paste0(gene, "_optimal")]] <- 
    ifelse(survival_data[[gene]] > optimal$cutpoint[[gene]], "High", "Low")
  
  cox_formula <- as.formula(
    paste0("Surv(os_time, os_status) ~ ", gene, "_optimal")
  )
  
  cox_model <- coxph(cox_formula, data = survival_data)
  cox_summary <- summary(cox_model)
  
  cat(sprintf("HR = %.2f (%.2f-%.2f), p = %.4f\n",
              cox_summary$conf.int[1, 1],
              cox_summary$conf.int[1, 3],
              cox_summary$conf.int[1, 4],
              cox_summary$coefficients[1, 5]))
  
  return(list(
    cutpoint = optimal$cutpoint[[gene]],
    statistic = optimal$statistic[[gene]],
    survival_data = survival_data,
    cox_model = cox_model
  ))
}


time_dependent_roc <- function(survival_data, 
                               gene, 
                               times = c(365, 1095, 1825)) {
  
  library(timeROC)
  
  if (!gene %in% colnames(survival_data)) {
    return(NULL)
  }
  
  
  roc_obj <- timeROC(
    T = survival_data$os_time,
    delta = survival_data$os_status,
    marker = survival_data[[gene]],
    cause = 1,
    times = times,
    iid = TRUE
  )
  

  for (i in seq_along(times)) {
    cat(sprintf("  %d days (%.1f years): %.3f\n", 
                times[i], 
                times[i] / 365.25,
                roc_obj$AUC[i]))
  }
  
  return(roc_obj)
}


subgroup_survival_analysis <- function(survival_data,
                                       gene,
                                       subgroup_var,
                                       cancer_type = "") {
  
  library(survival)
  
  if (!subgroup_var %in% colnames(survival_data)) {
    cat("WARNING: Subgroup variable not found:", subgroup_var, "\n")
    return(NULL)
  }
  
  subgroups <- unique(survival_data[[subgroup_var]])
  subgroups <- subgroups[!is.na(subgroups)]
  
  cat("Subgroups:", paste(subgroups, collapse = ", "), "\n\n")
  
  results <- data.frame()
  
  for (sg in subgroups) {
    
    sg_data <- survival_data[survival_data[[subgroup_var]] == sg, ]
    
    if (nrow(sg_data) < 10) {
      next
    }
    
    cox_formula <- as.formula(paste0("Surv(os_time, os_status) ~ ", gene))
    
    tryCatch({
      cox_model <- coxph(cox_formula, data = sg_data)
      cox_summary <- summary(cox_model)
      
      results <- rbind(
        results,
        data.frame(
          cancer_type = cancer_type,
          gene = gene,
          subgroup_var = subgroup_var,
          subgroup = sg,
          n_samples = cox_summary$n,
          n_events = cox_summary$nevent,
          HR = cox_summary$conf.int[1, 1],
          HR_lower = cox_summary$conf.int[1, 3],
          HR_upper = cox_summary$conf.int[1, 4],
          pvalue = cox_summary$coefficients[1, 5]
        )
      )
      
      cat(sprintf("  HR = %.2f (%.2f-%.2f), p = %.4f\n",
                  cox_summary$conf.int[1, 1],
                  cox_summary$conf.int[1, 3],
                  cox_summary$conf.int[1, 4],
                  cox_summary$coefficients[1, 5]))
      
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    })
  }
  
  return(results)
}


run_complete_survival_analysis <- function(data,
                                           normalized_counts,
                                           genes,
                                           cancer_type = "",
                                           output_dir = "results/survival") {
  

  cat("Complete Survival Analysis:", cancer_type, "\n")
  
  survival_data <- prepare_survival_data(data, normalized_counts, genes)
  
  if (is.null(survival_data) || nrow(survival_data) < 10) {
    cat("Insufficient survival data, skipping analysis\n")
    return(NULL)
  }
  
  univariate_results <- run_univariate_cox(survival_data, genes, cancer_type)
  
  multivariate_results <- run_multivariate_cox(
    survival_data, 
    genes,
    covariates = c("age_at_diagnosis", "tumor_stage", "gender"),
    cancer_type
  )
  
  dir.create(file.path(output_dir, cancer_type), 
             recursive = TRUE, showWarnings = FALSE)
  
  if (!is.null(univariate_results)) {
    sig_genes <- univariate_results$gene[univariate_results$significant]
    
    for (gene in sig_genes) {
      output_file <- file.path(output_dir, cancer_type, 
                               paste0(gene, "_KM_curve.pdf"))
      create_km_curve(survival_data, gene, output_file, cancer_type)
    }
  }
  
  return(list(
    survival_data = survival_data,
    univariate = univariate_results,
    multivariate = multivariate_results
  ))
}