# =============================================================================
# BRCA — GO ORA 2×2 Facet Plot (Activated & Suppressed ayrı görseller)
#
# Veri kaynağı: per_gene_advanced klasöründeki *_GO_ORA_results.csv dosyaları
# GSEA veya RDS'e gerek YOK — direkt CSV okur.
#
# ÜRETİLEN GÖRSELLER:
#   [1] BRCA_GO_ORA_Activated_2x2.jpeg   — 4 varyant, Activated GO terms
#   [2] BRCA_GO_ORA_Suppressed_2x2.jpeg  — 4 varyant, Suppressed GO terms
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(patchwork)
})

# =============================================================================
# KONFİGÜRASYON
# =============================================================================
ORA_DIR         <- "results/BRCA/gsea_pathway_analysis/per_gene_advanced"
OUTPUT_DIR      <- "results/BRCA/gsea_pathway_analysis/per_gene_advanced"
HISTONE_SYMBOLS <- c("H2AX", "H2AZ1", "MACROH2A1", "MACROH2A2")
TOP_N_PER_ONT   <- 4

# Renk paletleri
COLORS_ACT <- c(
  "Biological Process" = "#52B788",
  "Cellular Component" = "#74C2E1",
  "Molecular Function" = "#F4A261"
)
COLORS_SUP <- c(
  "Biological Process" = "#9B5DE5",
  "Cellular Component" = "#F15BB5",
  "Molecular Function" = "#E63946"
)

# =============================================================================
# VERİ YÜKLE — CSV'lerden oku
# =============================================================================
cat("=== ORA CSV'leri yükleniyor ===\n")

ora_all <- list()
for (h in HISTONE_SYMBOLS) {
  csv_path <- file.path(ORA_DIR, h, paste0(h, "_GO_ORA_results.csv"))
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    # Sütun adlarını kontrol et ve standartlaştır
    # Beklenen: Term, Gene_Number, Ontology, Direction, p.adjust
    cat(sprintf("  ✓ %s: %d term (%d Activated, %d Suppressed)\n",
                h, nrow(df),
                sum(df$Direction == "Activated", na.rm = TRUE),
                sum(df$Direction == "Suppressed", na.rm = TRUE)))
    ora_all[[h]] <- df
  } else {
    cat(sprintf("  ✗ %s: CSV bulunamadı — %s\n", h, csv_path))
  }
}

# =============================================================================
# YARDIMCI FONKSİYONLAR
# =============================================================================
wrap_term <- function(x, width = 36) stringr::str_wrap(x, width = width)

make_panel <- function(ora_df, gene_name, direction, color_map,
                       top_n = TOP_N_PER_ONT) {

  df <- ora_df |>
    dplyr::filter(Direction == direction) |>
    dplyr::filter(!is.na(Gene_Number), Gene_Number > 0) |>
    dplyr::mutate(
      Ontology = factor(Ontology,
                        levels = c("Biological Process",
                                   "Cellular Component",
                                   "Molecular Function"))
    ) |>
    dplyr::group_by(Ontology) |>
    dplyr::arrange(dplyr::desc(Gene_Number)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("no significant terms"),
                 size = 3.5, color = "grey55", hjust = 0.5, vjust = 0.5) +
        theme_void() +
        labs(title = gene_name) +
        theme(
          plot.title      = element_text(face = "bold", size = 11,
                                         hjust = 0.5, color = "grey15"),
          plot.background = element_rect(fill = "white", color = NA)
        )
    )
  }

  # Sıralama: önce ontoloji grubu, sonra Gene_Number (artan)
  df <- df |>
    dplyr::arrange(Ontology, Gene_Number) |>
    dplyr::mutate(
      TermWrap = wrap_term(Term, width = 36),
      TermWrap = factor(TermWrap, levels = unique(TermWrap))
    )

  x_max <- max(df$Gene_Number, na.rm = TRUE) * 1.30

  ggplot(df, aes(x = Gene_Number, y = TermWrap, fill = Ontology)) +
    geom_col(width = 0.68, alpha = 0.93, color = NA) +
    geom_text(aes(label = Gene_Number),
              hjust = -0.22, size = 2.8, color = "grey22") +
    scale_fill_manual(values = color_map,
                      name   = "GO Ontology",
                      drop   = FALSE) +
    scale_x_continuous(
      limits = c(0, x_max),
      expand = expansion(mult = c(0, 0.02)),
      name   = "Gene Number"
    ) +
    labs(title = gene_name, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 11,
                                        hjust = 0.5, color = "grey15"),
      axis.text.y        = element_text(size = 7.8, color = "grey20"),
      axis.text.x        = element_text(size = 8),
      axis.title.x       = element_text(size = 8.5),
      axis.ticks.y       = element_blank(),
      panel.grid.major.x = element_line(color = "grey91", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey75", linewidth = 0.45),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8.5),
      legend.key.size    = unit(0.5, "cm"),
      plot.background    = element_rect(fill = "white", color = NA),
      plot.margin        = margin(8, 10, 6, 6)
    )
}

# =============================================================================
# PLOT ÜRETİCİ — direction ve renk paletini parametre al
# =============================================================================
build_2x2 <- function(ora_list, direction, color_map, main_title, subtitle) {

  panels <- lapply(HISTONE_SYMBOLS, function(h) {
    df <- ora_list[[h]]
    if (is.null(df)) {
      # Veri yoksa boş panel
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste0(h, "\nno data"),
                   size = 3.5, color = "grey55", hjust = 0.5, vjust = 0.5) +
          theme_void() +
          labs(title = h) +
          theme(plot.title = element_text(face = "bold", size = 11,
                                          hjust = 0.5, color = "grey15"))
      )
    }
    make_panel(df, h, direction, color_map)
  })
  names(panels) <- HISTONE_SYMBOLS

  # scale_fill_manual tüm panellerde aynı olacak şekilde güncelle
  # (patchwork guides="collect" için gerekli)
  panels <- lapply(panels, function(p) {
    p + scale_fill_manual(values = color_map,
                          name   = "GO Ontology",
                          drop   = FALSE) +
      theme(legend.position = "right",
            legend.title    = element_text(face = "bold", size = 9),
            legend.text     = element_text(size = 8.5),
            legend.key.size = unit(0.5, "cm"))
  })

  (panels[[1]] | panels[[2]]) /
    (panels[[3]] | panels[[4]]) +
    plot_annotation(
      title    = main_title,
      subtitle = subtitle,
      theme    = theme(
        plot.title      = element_text(face = "bold", size = 14,
                                       hjust = 0.5, color = "grey10"),
        plot.subtitle   = element_text(color = "grey40", size = 9,
                                       hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

# =============================================================================
# [1] ACTIVATED 2×2
# =============================================================================
cat("\n--- Activated 2×2 oluşturuluyor ---\n")

p_act <- build_2x2(
  ora_list   = ora_all,
  direction  = "Activated",
  color_map  = COLORS_ACT,
  main_title = "Most Enriched GO Terms — Activated",
  subtitle   = "TCGA-BRCA | Positively correlated genes (Spearman r > 0.20) | adj.p < 0.05"
)

out_act <- file.path(OUTPUT_DIR, "BRCA_GO_ORA_Activated_2x2.jpeg")
ggsave(out_act, p_act,
       width = 16, height = 13, dpi = 240,
       units = "in", limitsize = FALSE)
cat("  →", basename(out_act), "\n")

# =============================================================================
# [2] SUPPRESSED 2×2
# =============================================================================
cat("\n--- Suppressed 2×2 oluşturuluyor ---\n")

p_sup <- build_2x2(
  ora_list   = ora_all,
  direction  = "Suppressed",
  color_map  = COLORS_SUP,
  main_title = "Most Enriched GO Terms — Suppressed",
  subtitle   = "TCGA-BRCA | Negatively correlated genes (Spearman r < -0.20) | adj.p < 0.05"
)

out_sup <- file.path(OUTPUT_DIR, "BRCA_GO_ORA_Suppressed_2x2.jpeg")
ggsave(out_sup, p_sup,
       width = 16, height = 13, dpi = 240,
       units = "in", limitsize = FALSE)
cat("  →", basename(out_sup), "\n")

cat("\n✓ Tamamlandı.\n")
cat("Çıktılar:\n")
cat(" •", out_act, "\n")
cat(" •", out_sup, "\n")
