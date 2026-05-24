# =============================================================================
# BRCA — Gen Başına GO/KEGG Advanced Görselleştirme (Temiz Versiyon)
# Genler: H2AX, H2AZ1, MACROH2A1, MACROH2A2
#
# Ön koşul: ana BRCA GSEA scripti çalıştırılmış olmalı
#           all_gsea_results workspace'de mevcut olmalı
# Bağımsız çalıştırma: PATH_STANDALONE$gsea_rds yolunu doldurun
#
# ÜRETİLEN GÖRSELLER (her gen için):
#   [1] GO NES Bar Plot      — NES gradyan bar, kırmızı/mavi
#   [2] GO Lollipop          — NES lollipop, renk = -log10(padj)
#   [3] GO Running Score     — enrichplot::gseaplot2, p-value tablosu
#   [4] KEGG Dot Plot        — activated | suppressed facet, GeneRatio
#   [5] KEGG NES Bar Plot    — NES gradyan bar
#   [6] KEGG Running Score   — enrichplot::gseaplot2, p-value tablosu
#   [7] GO ORA Bar Plot      — makale stili, BP/CC/MF renkli, Gene_Number
# =============================================================================

# Namespace çakışmalarını önlemek için tüm paketleri yüklemeden önce
# dplyr'ı en sona yükle — böylece select/filter/mutate dplyr'dan gelir
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(dplyr)   # EN SONA — select/filter override için
})

# =============================================================================
# KONFİGÜRASYON
# =============================================================================
PATH_STANDALONE <- list(
  gsea_rds = "results/BRCA/gsea_pathway_analysis/all_gsea_results_BRCA.rds"
)

OUTPUT_DIR     <- "results/BRCA/gsea_pathway_analysis/per_gene_advanced"
HISTONE_SYMBOLS <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# VERİ YÜKLEYİN
# =============================================================================
if (!exists("all_gsea_results")) {
  if (file.exists(PATH_STANDALONE$gsea_rds)) {
    cat("GSEA sonuçları yükleniyor:", PATH_STANDALONE$gsea_rds, "\n")
    all_gsea_results <- readRDS(PATH_STANDALONE$gsea_rds)
  } else {
    stop("all_gsea_results bulunamadı. Ana scripti çalıştırın veya ",
         "PATH_STANDALONE$gsea_rds yolunu doldurun.")
  }
}

# KEGG setleri (bağımsız çalışma için)
if (!exists("kegg_sets")) {
  cat("KEGG setleri yükleniyor...\n")
  all_col  <- msigdbr_collections()
  sub_col  <- if ("gs_subcollection" %in% names(all_col)) "gs_subcollection" else "gs_subcat"
  cat_col  <- if ("gs_collection"    %in% names(all_col)) "gs_collection"    else "gs_cat"
  kegg_subs <- all_col[[sub_col]][
    grepl("KEGG", all_col[[sub_col]], ignore.case = TRUE) &
      all_col[[cat_col]] == "C2"
  ]
  kegg_sub_use <- if ("CP:KEGG_LEGACY" %in% kegg_subs) "CP:KEGG_LEGACY" else kegg_subs[1]
  kegg_sets <- msigdbr(species = "Homo sapiens",
                       category = "C2", subcategory = kegg_sub_use) |>
    dplyr::select(gs_name, entrez_gene) |>
    dplyr::mutate(entrez_gene = as.character(entrez_gene))
  cat("  KEGG:", length(unique(kegg_sets$gs_name)), "set\n\n")
}

cat("=== BRCA Per-Gene GO/KEGG Advanced Görselleştirme ===\n\n")

# =============================================================================
# YARDIMCI FONKSİYONLAR
# =============================================================================

clean_term <- function(x, max_char = 55) {
  x <- gsub("^KEGG_|^REACTOME_|^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  ifelse(nchar(x) > max_char, paste0(substr(x, 1, max_char - 2), ".."), x)
}

wrap_term <- function(x, width = 38) stringr::str_wrap(x, width = width)

calc_gene_ratio <- function(df) {
  n_core <- lengths(strsplit(df$core_enrichment, "/"))
  n_core / df$setSize
}

# Renk sabitleri
NES_COLORS  <- c("#2166AC", "#6BAED6", "#D1E5F0", "white",
                 "#FDDBC7", "#EF6548", "#B2182B")
ONT_COLORS  <- c("Biological Process" = "#CCE2CB",
                 "Cellular Component" = "#BCCFDF",
                 "Molecular Function" = "#FBC4B6")
ONT_LABELS  <- c(BP = "Biological Process",
                 CC = "Cellular Component",
                 MF = "Molecular Function")

# Ortak theme
theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(
      plot.title       = element_text(face = "bold", size = base + 2, hjust = 0.5),
      plot.subtitle    = element_text(color = "grey40", size = base - 2, hjust = 0.5),
      axis.text.y      = element_text(size = base - 1.5, color = "grey20"),
      axis.text.x      = element_text(size = base - 2),
      axis.ticks.y     = element_blank(),
      legend.title     = element_text(face = "bold", size = base - 2),
      legend.text      = element_text(size = base - 2.5),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# =============================================================================
# GÖRSEL FONKSİYONLARI
# =============================================================================

# ── 1. NES Bar Plot (GO veya KEGG için ortak) ─────────────────────────────
plot_nes_bar <- function(df, title, top_n = 10) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  plot_df <- dplyr::bind_rows(
    dplyr::filter(df, NES > 0) |> dplyr::arrange(dplyr::desc(NES)) |> dplyr::slice_head(n = top_n),
    dplyr::filter(df, NES < 0) |> dplyr::arrange(NES)              |> dplyr::slice_head(n = top_n)
  ) |>
    dplyr::mutate(Term = factor(Term, levels = Term[order(NES)]))

  if (nrow(plot_df) == 0) return(NULL)

  lim <- max(abs(plot_df$NES)) * 1.08

  ggplot(plot_df, aes(x = NES, y = Term, fill = NES)) +
    geom_col(width = 0.72, alpha = 0.95) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", NES),
                  hjust = ifelse(NES > 0, -0.15, 1.15)),
              size = 2.9, color = "grey25") +
    scale_fill_gradientn(
      colors = NES_COLORS,
      values = scales::rescale(c(-lim, -1, -0.3, 0, 0.3, 1, lim)),
      limits = c(-lim, lim), name = "NES"
    ) +
    scale_x_continuous(limits = c(-lim * 1.25, lim * 1.25),
                       expand = expansion(mult = 0)) +
    labs(title = title, subtitle = "Normalized Enrichment Score (NES)",
         x = "NES", y = NULL) +
    theme_pub() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "grey50", linewidth = 0.4),
      legend.key.height  = unit(1.2, "cm")
    )
}

# ── 2. GO Lollipop ────────────────────────────────────────────────────────
plot_go_lollipop <- function(df, title, top_n = 30) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  plot_df <- dplyr::bind_rows(
    dplyr::filter(df, NES > 0) |> dplyr::slice_max(NES, n = ceiling(top_n / 2)),
    dplyr::filter(df, NES < 0) |> dplyr::slice_min(NES, n = floor(top_n / 2))
  ) |>
    dplyr::mutate(Term = factor(Term, levels = Term[order(NES)]))

  ggplot(plot_df, aes(x = NES, y = Term, color = log10padj)) +
    geom_segment(aes(xend = 0, yend = Term), color = "grey80", linewidth = 0.55) +
    geom_point(aes(size = setSize), alpha = 0.92) +
    geom_vline(xintercept = 0, color = "grey35", linewidth = 0.45, linetype = "dashed") +
    scale_color_gradientn(
      colors = c("#4393C3", "#92C5DE", "#FDDBC7", "#EF6548", "#B2182B"),
      name = expression(-log[10](adj.p))
    ) +
    scale_size_continuous(range = c(2.5, 9), name = "Set Size",
                          guide = guide_legend(override.aes = list(color = "grey50"))) +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    labs(title = title,
         subtitle = paste0("Top ", nrow(plot_df), " pathway | adj.p < 0.05"),
         x = "NES", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey92"))
}

# ── 3. KEGG Dot Plot (activated | suppressed facet) ───────────────────────
plot_kegg_dot <- function(df, title, top_n = 20) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  plot_df <- dplyr::bind_rows(
    dplyr::filter(df, NES > 0) |> dplyr::arrange(dplyr::desc(GeneRatio)) |> dplyr::slice_head(n = top_n),
    dplyr::filter(df, NES < 0) |> dplyr::arrange(dplyr::desc(GeneRatio)) |> dplyr::slice_head(n = top_n)
  ) |>
    dplyr::mutate(
      Count     = lengths(strsplit(core_enrichment, "/")),
      Direction = factor(ifelse(NES > 0, "activated", "suppressed"),
                         levels = c("activated", "suppressed")),
      TermWrap  = wrap_term(Term, width = 30),
      TermWrap  = factor(TermWrap, levels = unique(TermWrap[order(GeneRatio)]))
    )

  padj_lim <- c(min(plot_df$p.adjust, na.rm = TRUE),
                min(0.05, max(plot_df$p.adjust, na.rm = TRUE)))

  ggplot(plot_df, aes(x = GeneRatio, y = TermWrap, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.90) +
    facet_grid(. ~ Direction, scales = "free_x", space = "free_x") +
    scale_color_gradientn(
      colors = c("#B2182B", "#EF6548", "#FDDBC7", "#D1E5F0", "#6BAED6", "#2166AC"),
      limits = padj_lim, oob = scales::squish,
      name = "p.adjust", labels = scales::scientific_format(digits = 2),
      guide = guide_colorbar(barheight = unit(3.5, "cm"))
    ) +
    scale_size_continuous(range = c(3, 11), name = "Count",
                          breaks = c(10, 30, 60, 100),
                          guide = guide_legend(override.aes = list(color = "grey50"))) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(title = title,
         subtitle = paste0("GeneRatio | Top ", top_n, "/yön | adj.p < 0.05"),
         x = "GeneRatio", y = NULL) +
    theme_pub(base = 10.5) +
    theme(
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey94", color = "grey80"),
      panel.grid.major.y = element_line(color = "grey93")
    )
}

# ── 4. GO ORA Bar Plot — makale stili (BP/CC/MF renkli) ───────────────────
plot_go_ora_bar <- function(df, gene_name, direction = "Activated",
                            top_n_per_ont = 8, subtitle_extra = "") {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  plot_df <- df |>
    dplyr::filter(if ("Direction" %in% names(df)) Direction == direction else TRUE) |>
    dplyr::filter(!is.na(Gene_Number), Gene_Number > 0) |>
    dplyr::group_by(Ontology) |>
    dplyr::arrange(dplyr::desc(Gene_Number)) |>
    dplyr::slice_head(n = top_n_per_ont) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      TermWrap = wrap_term(Term, width = 40),
      Ontology = factor(Ontology, levels = c("Biological Process",
                                             "Cellular Component",
                                             "Molecular Function"))
    ) |>
    dplyr::arrange(Ontology, Gene_Number) |>
    dplyr::mutate(TermWrap = factor(TermWrap, levels = unique(TermWrap)))

  if (nrow(plot_df) == 0) return(NULL)

  x_max    <- max(plot_df$Gene_Number) * 1.22
  subtitle <- if (nchar(subtitle_extra) > 0) {
    paste0("TCGA-BRCA | ", direction, " | ", subtitle_extra)
  } else {
    paste0("TCGA-BRCA | ", direction, " GO terms | adj.p < 0.05")
  }

  ggplot(plot_df, aes(x = Gene_Number, y = TermWrap, fill = Ontology)) +
    geom_col(width = 0.72, alpha = 0.92, color = NA) +
    geom_text(aes(label = Gene_Number), hjust = -0.25,
              size = 3.0, color = "grey25") +
    scale_fill_manual(values = ONT_COLORS, name = "Type",
                      breaks = names(ONT_COLORS)) +
    scale_x_continuous(limits = c(0, x_max),
                       expand = expansion(mult = c(0, 0.02)),
                       name = "Gene_Number") +
    labs(title    = paste0("Most enriched GO terms \u2014 ", gene_name),
         subtitle = subtitle, y = "GO term") +
    theme_pub() +
    theme(
      plot.title         = element_text(face = "plain", size = 12,
                                        hjust = 0.5, color = "grey15"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.border       = element_rect(color = "grey75", linewidth = 0.5),
      legend.key.size    = unit(0.55, "cm"),
      legend.box.background = element_rect(color = "grey80", linewidth = 0.4,
                                           fill = "white"),
      plot.margin        = margin(12, 14, 10, 10)
    )
}

# =============================================================================
# ORA FONKSİYONU
# =============================================================================

run_go_ora <- function(gene_name, gsea_results, top_n = 10,
                       padj_cut = 0.05, cor_threshold = 0.20) {
  res <- gsea_results[[gene_name]]
  if (is.null(res) || is.null(res$gene_list)) return(NULL)

  gl         <- res$gene_list
  pos_genes  <- names(gl[gl >  cor_threshold])
  neg_genes  <- names(gl[gl < -cor_threshold])
  universe   <- names(gl)

  cat(sprintf("  ORA: %d pozitif + %d negatif gen (|r| > %.2f)\n",
              length(pos_genes), length(neg_genes), cor_threshold))

  run_one <- function(gene_set, ont, direction) {
    if (length(gene_set) < 10) return(NULL)
    tryCatch({
      obj <- enrichGO(
        gene          = gene_set,
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cut,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )
      if (is.null(obj) || nrow(obj@result) == 0) return(NULL)

      obj@result |>
        dplyr::filter(p.adjust < padj_cut) |>
        dplyr::mutate(
          Gene_Number = Count,
          Ontology    = ONT_LABELS[ont],
          Direction   = direction,
          Term        = Description
        ) |>
        dplyr::arrange(dplyr::desc(Gene_Number)) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::select(Term, Gene_Number, Ontology, Direction, p.adjust, Count, GeneRatio)
    }, error = function(e) {
      cat(sprintf("    enrichGO %s %s hatası: %s\n", ont, direction, conditionMessage(e)))
      NULL
    })
  }

  dplyr::bind_rows(
    run_one(pos_genes, "BP", "Activated"),
    run_one(pos_genes, "CC", "Activated"),
    run_one(pos_genes, "MF", "Activated"),
    run_one(neg_genes, "BP", "Suppressed"),
    run_one(neg_genes, "CC", "Suppressed"),
    run_one(neg_genes, "MF", "Suppressed")
  )
}

# =============================================================================
# KAYDETME YARDIMCISI
# =============================================================================

save_plot <- function(p, path, width, height, dpi = 240) {
  if (is.null(p)) return(invisible(NULL))
  ggplot2::ggsave(path, p, width = width, height = height,
                  dpi = dpi, units = "in", limitsize = FALSE)
  cat("  →", basename(path), "\n")
}

# =============================================================================
# ANA DÖNGÜ
# =============================================================================

summary_list <- list()

for (h_name in HISTONE_SYMBOLS) {

  cat("\n╔══════════════════════════════════════╗\n")
  cat(sprintf("║  %-36s║\n", paste0(h_name, " — GO/KEGG Analizi")))
  cat("╚══════════════════════════════════════╝\n")

  gene_dir <- file.path(OUTPUT_DIR, h_name)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  res      <- all_gsea_results[[h_name]]
  go_obj   <- res$GO
  kegg_obj <- res$KEGG

  # ── GO-BP hazırlık ────────────────────────────────────────────────────────
  go_df <- NULL
  if (!is.null(go_obj) && nrow(go_obj@result) > 0) {
    go_df <- go_obj@result |>
      dplyr::filter(p.adjust < 0.05) |>
      dplyr::mutate(
        GeneRatio = calc_gene_ratio(dplyr::pick(dplyr::everything())),
        Term      = clean_term(Description),
        log10padj = -log10(p.adjust)
      )
    cat(sprintf("  GO-BP: %d anlamlı pathway (+ %d / - %d)\n",
                nrow(go_df), sum(go_df$NES > 0), sum(go_df$NES < 0)))
  }

  # ── KEGG hazırlık ─────────────────────────────────────────────────────────
  kegg_df <- NULL
  if (!is.null(kegg_obj) && nrow(kegg_obj@result) > 0) {
    kegg_df <- kegg_obj@result |>
      dplyr::filter(p.adjust < 0.05) |>
      dplyr::mutate(
        GeneRatio = calc_gene_ratio(dplyr::pick(dplyr::everything())),
        Term      = clean_term(ID),
        log10padj = -log10(p.adjust)
      )
    cat(sprintf("  KEGG:  %d anlamlı pathway (+ %d / - %d)\n",
                nrow(kegg_df), sum(kegg_df$NES > 0), sum(kegg_df$NES < 0)))
  }

  # ── [1] GO NES Bar ────────────────────────────────────────────────────────
  cat("\n── Görseller ──\n")
  if (!is.null(go_df)) {
    p <- plot_nes_bar(go_df,
                      title = paste0("GSEA GO \u2014 Top Pathways (BRCA / ", h_name, ")"))
    n <- nrow(dplyr::bind_rows(
      dplyr::filter(go_df, NES > 0) |> dplyr::slice_head(n = 10),
      dplyr::filter(go_df, NES < 0) |> dplyr::slice_head(n = 10)))
    save_plot(p, file.path(gene_dir, paste0(h_name, "_GO_barplot_BRCA.jpeg")),
              width = 10, height = max(5.5, n * 0.36 + 2.0))
  }

  # ── [2] GO Lollipop ───────────────────────────────────────────────────────
  if (!is.null(go_df)) {
    p <- plot_go_lollipop(go_df,
                          title = paste0("BRCA / ", h_name, ": GO-BP GSEA"))
    n <- min(30, nrow(go_df))
    save_plot(p, file.path(gene_dir, paste0(h_name, "_GO_lollipop_BRCA.jpeg")),
              width = 10, height = max(5, n * 0.32 + 2.5))
  }

  # ── [3] GO Running Score ──────────────────────────────────────────────────
  if (!is.null(go_obj) && nrow(go_obj@result) >= 2) {
    top_ids <- go_obj@result |>
      dplyr::filter(p.adjust < 0.05) |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = 4) |>
      dplyr::pull(ID)

    if (length(top_ids) >= 2) {
      tryCatch({
        cols <- RColorBrewer::brewer.pal(max(3, length(top_ids)), "Set1")[1:length(top_ids)]
        p    <- enrichplot::gseaplot2(go_obj, geneSetID = top_ids, color = cols,
                                      pvalue_table = TRUE, ES_geom = "line",
                                      rel_heights = c(1.5, 0.4, 0.8))
        save_plot(p, file.path(gene_dir, paste0(h_name, "_GO_running_score_BRCA.jpeg")),
                  width = 11, height = 7, dpi = 220)
      }, error = function(e) cat("  GO running score atlandı:", conditionMessage(e), "\n"))
    }
  }

  # ── [4] KEGG Dot Plot ─────────────────────────────────────────────────────
  if (!is.null(kegg_df)) {
    p <- plot_kegg_dot(kegg_df,
                       title = paste0("GSEA KEGG Pathways \u2014 BRCA / ", h_name))
    n <- max(sum(kegg_df$NES > 0), sum(kegg_df$NES < 0))
    n <- min(n, 20)
    save_plot(p, file.path(gene_dir, paste0(h_name, "_KEGG_dotplot_BRCA.jpeg")),
              width = 13, height = max(6, n * 0.38 + 3.5))
  }

  # ── [5] KEGG NES Bar ─────────────────────────────────────────────────────
  if (!is.null(kegg_df)) {
    p <- plot_nes_bar(kegg_df,
                      title = paste0("GSEA KEGG \u2014 Top Pathways (BRCA / ", h_name, ")"))
    n <- nrow(dplyr::bind_rows(
      dplyr::filter(kegg_df, NES > 0) |> dplyr::slice_head(n = 10),
      dplyr::filter(kegg_df, NES < 0) |> dplyr::slice_head(n = 10)))
    save_plot(p, file.path(gene_dir, paste0(h_name, "_KEGG_barplot_BRCA.jpeg")),
              width = 10, height = max(5, n * 0.36 + 2.0))
  }

  # ── [6] KEGG Running Score ────────────────────────────────────────────────
  if (!is.null(kegg_obj) && nrow(kegg_obj@result) >= 1) {
    act_ids  <- kegg_obj@result |>
      dplyr::filter(p.adjust < 0.05, NES > 0) |>
      dplyr::arrange(p.adjust) |> dplyr::slice_head(n = 2) |> dplyr::pull(ID)
    supp_ids <- kegg_obj@result |>
      dplyr::filter(p.adjust < 0.05, NES < 0) |>
      dplyr::arrange(p.adjust) |> dplyr::slice_head(n = 2) |> dplyr::pull(ID)
    run_ids  <- c(act_ids, supp_ids)

    if (length(run_ids) >= 1) {
      tryCatch({
        cols <- RColorBrewer::brewer.pal(max(3, length(run_ids)), "Set1")[1:length(run_ids)]
        p    <- enrichplot::gseaplot2(kegg_obj, geneSetID = run_ids, color = cols,
                                      pvalue_table = TRUE, ES_geom = "line",
                                      rel_heights = c(1.5, 0.4, 0.8))
        save_plot(p, file.path(gene_dir, paste0(h_name, "_KEGG_running_score_BRCA.jpeg")),
                  width = 12, height = 7.5, dpi = 220)
      }, error = function(e) cat("  KEGG running score atlandı:", conditionMessage(e), "\n"))
    }
  }

  # ── [7] GO ORA Bar Plot — makale stili ────────────────────────────────────
  cat("  ORA çalışıyor...\n")
  ora_df <- run_go_ora(h_name, all_gsea_results, top_n = 10,
                       padj_cut = 0.05, cor_threshold = 0.20)

  if (!is.null(ora_df) && nrow(ora_df) > 0) {
    cat(sprintf("  ORA: %d term (Activated: %d, Suppressed: %d)\n",
                nrow(ora_df),
                sum(ora_df$Direction == "Activated"),
                sum(ora_df$Direction == "Suppressed")))

    write.csv(ora_df,
              file.path(gene_dir, paste0(h_name, "_GO_ORA_results.csv")),
              row.names = FALSE)

    # Activated
    p <- plot_go_ora_bar(ora_df, h_name, "Activated", 8,
                         "Positively correlated genes")
    n <- ora_df |> dplyr::filter(Direction == "Activated") |>
      dplyr::group_by(Ontology) |> dplyr::slice_head(n = 8) |> nrow()
    save_plot(p, file.path(gene_dir, paste0(h_name, "_GO_ORA_activated_BRCA.jpeg")),
              width = 9.5, height = max(5.5, n * 0.38 + 2.5))

    # Suppressed
    p <- plot_go_ora_bar(ora_df, h_name, "Suppressed", 8,
                         "Negatively correlated genes")
    n <- ora_df |> dplyr::filter(Direction == "Suppressed") |>
      dplyr::group_by(Ontology) |> dplyr::slice_head(n = 8) |> nrow()
    save_plot(p, file.path(gene_dir, paste0(h_name, "_GO_ORA_suppressed_BRCA.jpeg")),
              width = 9.5, height = max(5.5, n * 0.38 + 2.5))

    # Kombine facet
    comb_df <- ora_df |>
      dplyr::group_by(Direction, Ontology) |>
      dplyr::arrange(dplyr::desc(Gene_Number)) |>
      dplyr::slice_head(n = 5) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        TermWrap  = wrap_term(Term, 38),
        Ontology  = factor(Ontology, levels = c("Biological Process",
                                                "Cellular Component",
                                                "Molecular Function")),
        Direction = factor(Direction, levels = c("Activated", "Suppressed"))
      ) |>
      dplyr::arrange(Direction, Ontology, Gene_Number) |>
      dplyr::mutate(TermWrap = factor(TermWrap, levels = rev(unique(TermWrap))))

    if (nrow(comb_df) > 0) {
      x_max <- max(comb_df$Gene_Number) * 1.22

      p_comb <- ggplot(comb_df, aes(x = Gene_Number, y = TermWrap, fill = Ontology)) +
        geom_col(width = 0.70, alpha = 0.92, color = NA) +
        geom_text(aes(label = Gene_Number), hjust = -0.25, size = 2.8, color = "grey25") +
        facet_grid(Direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
        scale_fill_manual(values = ONT_COLORS, name = "Type", breaks = names(ONT_COLORS)) +
        scale_x_continuous(limits = c(0, x_max),
                           expand = expansion(mult = c(0, 0.02)),
                           name = "Gene_Number") +
        labs(title    = paste0("Most enriched GO terms \u2014 ", h_name),
             subtitle = "TCGA-BRCA | Activated & Suppressed correlated genes",
             y = "GO term") +
        theme_pub(base = 10.5) +
        theme(
          plot.title         = element_text(face = "plain", size = 12,
                                            hjust = 0.5, color = "grey15"),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
          panel.grid.major.y = element_blank(),
          panel.border       = element_rect(color = "grey75", linewidth = 0.5),
          strip.text.y       = element_text(face = "bold", size = 10, angle = 180),
          strip.placement    = "outside",
          strip.background   = element_rect(fill = "grey94", color = "grey80"),
          legend.key.size    = unit(0.55, "cm"),
          legend.box.background = element_rect(color = "grey80", linewidth = 0.4,
                                               fill = "white"),
          plot.margin        = margin(12, 14, 10, 10)
        )

      save_plot(p_comb,
                file.path(gene_dir, paste0(h_name, "_GO_ORA_combined_BRCA.jpeg")),
                width = 9.5, height = max(7, nrow(comb_df) * 0.34 + 3.5))
    }

  } else {
    cat("  ORA sonuç yok — gen listesi eksik veya eşik altında.\n")
  }

  # CSV
  if (!is.null(go_df))
    write.csv(go_df,   file.path(gene_dir, paste0(h_name, "_GO_results.csv")),   row.names = FALSE)
  if (!is.null(kegg_df))
    write.csv(kegg_df, file.path(gene_dir, paste0(h_name, "_KEGG_results.csv")), row.names = FALSE)

  # Özet
  summary_list[[h_name]] <- data.frame(
    Gene           = h_name,
    GO_total       = if (!is.null(go_df))   nrow(go_df)               else 0,
    GO_activated   = if (!is.null(go_df))   sum(go_df$NES > 0)        else 0,
    GO_suppressed  = if (!is.null(go_df))   sum(go_df$NES < 0)        else 0,
    KEGG_total     = if (!is.null(kegg_df)) nrow(kegg_df)             else 0,
    KEGG_activated = if (!is.null(kegg_df)) sum(kegg_df$NES > 0)      else 0,
    KEGG_suppressed= if (!is.null(kegg_df)) sum(kegg_df$NES < 0)      else 0
  )
}


summary_df <- dplyr::bind_rows(summary_list)
print(summary_df, row.names = FALSE)
write.csv(summary_df, file.path(OUTPUT_DIR, "BRCA_perGene_summary.csv"), row.names = FALSE)

