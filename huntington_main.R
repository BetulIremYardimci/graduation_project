BiocManager::install("GEOquery", force = TRUE)
library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)

#Expression

library(DESeq2)
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE64810&format=file&file=GSE64810_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  destfile = "GSE64810_raw_counts.tsv.gz",
  mode = "wb"
)

geo_id <- "GSE64810"

raw_counts <- read.table(gzfile("GSE64810_raw_counts.tsv.gz"),
                         header = TRUE,
                         row.names = 1,
                         sep = "\t")
write.table(raw_counts, "huntington/data/GSE64810_raw_counts.tsv")

dim(raw_counts)
head(raw_counts[1:3, 1:3])

sample_names <- colnames(raw_counts)
head(sample_names)

group <- ifelse(substr(colnames(counts), 1, 1) == "C", "Control", "HD")
table(group)

metadata <- data.frame(
  sample = colnames(raw_counts),
  condition = factor(group, levels = c("Control", "HD"))
)
rownames(metadata) <- colnames(raw_counts)
print(metadata)

raw_counts_mat <- as.matrix(raw_counts)
raw_counts_int <- round(raw_counts_mat)
storage.mode(raw_counts_int) <- "integer"

class(raw_counts_int)
head(raw_counts_int[1:3, 1:3])

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_int,
  colData   = metadata,
  design    = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
cat("Filtered genes remaining:", nrow(dds), "\n")

dds <- DESeq(dds)

res <- results(dds,
               contrast = c("condition", "HD", "Control"),
               alpha = 0.05)

summary(res)

histone_entrez <- c(
  H2AFX  = "3014",
  H2AFY  = "9555",
  H2AFY2 = "55506",
  H2AFZ  = "3015"
)

res_df <- as.data.frame(res)
res_df$entrez <- rownames(res_df)

histone_results <- res_df[res_df$entrez %in% histone_entrez, ]
histone_results$gene <- names(histone_entrez)[match(
  histone_results$entrez, histone_entrez)]

print(histone_results[, c("gene", "log2FoldChange", "pvalue", "padj")])

write.csv(res_df, "huntington/results/deseq2_full_results.csv")


# Normalize counts'tan histone değerlerini al
# counts objesi zaten DESeq2 normalize - ama row names Ensembl ID
# raw_counts'tan vst alalım

vst_counts <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_counts)

# Histone Entrez ID'leri ile filtrele
histone_entrez <- c(
  H2AFX  = "3014",
  H2AFY  = "9555",
  H2AFY2 = "55506",
  H2AFZ  = "3015"
)

# Histone satırlarını çek
histone_vst <- vst_mat[rownames(vst_mat) %in% histone_entrez, ]

# Satır isimlerini gen sembolleriyle değiştir
rownames(histone_vst) <- names(histone_entrez)[match(
  rownames(histone_vst), histone_entrez)]

# Kontrol et
print(histone_vst[, 1:5])

# Long format'a çevir
library(tidyr)
library(dplyr)

histone_long <- as.data.frame(histone_vst) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(group = ifelse(substr(sample, 1, 10) %in%
                          colnames(raw_counts)[metadata$condition == "Control"],
                        "Control", "HD"))

head(histone_long)

library(ggplot2)
library(gridExtra)

# İstatistik hesapla
histone_stats <- histone_long %>%
  group_by(gene, group) %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr   = sd(expression),
    n         = n(),
    .groups   = "drop"
  )

# Significance hesapla (DESeq2 padj kullan)
sig_df <- data.frame(
  gene  = c("H2AFX", "H2AFY", "H2AFY2", "H2AFZ"),
  padj  = c(0.9797,  0.9573,  0.7487,   0.0071),
  label = c("ns",    "ns",    "ns",      "**")
)

# Renk
colors <- c("Control" = "#A8C4D4", "HD" = "#F2A8A8")

# Sample sayıları
n_ctrl <- sum(metadata$condition == "Control")
n_hd   <- sum(metadata$condition == "HD")

# Her gen için plot
plot_list <- list()

plot_list <- list()

plot_list <- list()

histone_entrez_map <- c("3014", "9555", "55506", "3015")
names(histone_entrez_map) <- c("H2AFX", "H2AFY", "H2AFY2", "H2AFZ")

sig_data <- res_df %>%
  filter(entrez %in% histone_entrez_map) %>%
  mutate(gene = names(histone_entrez_map)[match(entrez, histone_entrez_map)]) %>%
  mutate(label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ "ns"
  ))



# res_df icinden histone p-degerlerini bulup etiketleme yapalim
sig_data <- res_df %>%
  filter(entrez %in% histone_entrez_map) %>%
  mutate(gene = names(histone_entrez_map)[match(entrez, histone_entrez_map)]) %>%
  mutate(label = case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ "ns"
  ))

# --- GUNCEL PLOT DONGUSU ---
plot_list <- list()
colors <- c("Control" = "#A8C4D4", "HD" = "#F2A8A8")

histone_entrez_map <- c("3014", "9555", "55506", "3015")
names(histone_entrez_map) <- c("H2AFX", "H2AFY", "H2AFY2", "H2AFZ")

for(g in names(histone_entrez_map)) {

  df_plot   <- histone_stats[histone_stats$gene == g, ]
  sig_label <- sig_data$label[sig_data$gene == g]
  p_val_raw <- sig_data$padj[sig_data$gene == g]

  # Y ekseni limitlerini p-degeri etiketi icin ayarla
  y_max     <- max(df_plot$mean_expr + df_plot$sd_expr)
  y_bracket <- y_max * 1.05
  y_star    <- y_max * 1.12

  p <- ggplot(df_plot, aes(x = group, y = mean_expr, fill = group)) +
    geom_bar(stat = "identity", width = 0.55, alpha = 0.9,
             color = "white", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_expr - sd_expr,
                      ymax = mean_expr + sd_expr),
                  width = 0.15, linewidth = 0.7, color = "gray30") +
    # Bracket (Anlamlilik çizgisi)
    annotate("segment", x = 1, xend = 2,
             y = y_bracket, yend = y_bracket, linewidth = 0.6) +
    annotate("segment", x = 1, xend = 1,
             y = y_bracket * 0.98, yend = y_bracket, linewidth = 0.6) +
    annotate("segment", x = 2, xend = 2,
             y = y_bracket * 0.98, yend = y_bracket, linewidth = 0.6) +
    # Anlamlilik Yıldızı veya 'ns'
    annotate("text", x = 1.5, y = y_star,
             label = sig_label, size = 6, fontface = "bold") +
    scale_fill_manual(values = colors) +
    # P-degerini alt bilgi olarak eklemek istersen (isteğe bağlı)
    labs(subtitle = paste0("adj.p = ", format.pval(p_val_raw, digits = 2)),
         y = "VST Normalized Expression") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(g)

  plot_list[[g]] <- p
}

# Birleştir
final_plot <- gridExtra::arrangeGrob(
  title_grob, subtitle_grob,
  gridExtra::arrangeGrob(plot_list[["H2AFX"]], plot_list[["H2AFY"]], ncol = 2),
  gridExtra::arrangeGrob(plot_list[["H2AFY2"]], plot_list[["H2AFZ"]], ncol = 2),
  legend_grob,
  heights = c(0.06, 0.04, 1, 1, 0.12)
)

# Kaydet
ggsave("huntington/plots/histone_expression_barplot.png",
       final_plot, width = 12, height = 10, dpi = 300)


#CORRELATION
# VST matrisinden H2AFZ'yi al
h2afz_entrez <- "3015"
idx <- which(rownames(vst_mat) == h2afz_entrez)

h2afz_expr <- vst_mat[idx, ]

# Control ve HD gruplarını ayır
ctrl_samples <- colnames(vst_mat)[metadata$condition == "Control"]
hd_samples   <- colnames(vst_mat)[metadata$condition == "HD"]

h2afz_ctrl <- h2afz_expr[ctrl_samples]
h2afz_hd   <- h2afz_expr[hd_samples]

vst_ctrl <- vst_mat[, ctrl_samples]
vst_hd   <- vst_mat[, hd_samples]

# Spearman korelasyon - Control
cor_ctrl <- apply(vst_ctrl, 1, function(x) {
  cor.test(x, h2afz_ctrl, method = "spearman", exact = FALSE)$estimate
})

pval_ctrl <- apply(vst_ctrl, 1, function(x) {
  cor.test(x, h2afz_ctrl, method = "spearman", exact = FALSE)$p.value
})

# Spearman korelasyon - HD
cor_hd <- apply(vst_hd, 1, function(x) {
  cor.test(x, h2afz_hd, method = "spearman", exact = FALSE)$estimate
})

pval_hd <- apply(vst_hd, 1, function(x) {
  cor.test(x, h2afz_hd, method = "spearman", exact = FALSE)$p.value
})

# Dataframe oluştur
cor_df <- data.frame(
  entrez   = rownames(vst_mat),
  cor_ctrl = cor_ctrl,
  cor_hd   = cor_hd,
  pval_ctrl = pval_ctrl,
  pval_hd   = pval_hd
)

cor_df$padj_ctrl <- p.adjust(cor_df$pval_ctrl, method = "BH")
cor_df$padj_hd   <- p.adjust(cor_df$pval_hd,   method = "BH")

# Filtrele - her iki grupta da anlamlı
sig_both <- cor_df[
  (abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05) |
    (abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05),
]

cat("Control'de anlamlı:", sum(abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05), "\n")
cat("HD'de anlamlı:", sum(abs(cor_df$cor_hd) > 0.4 & cor_df$padj_hd < 0.05), "\n")
cat("Her ikisinde:", nrow(sig_both), "\n")


library(org.Hs.eg.db)

# Entrez ID'den sembol dönüşümü
symbols <- mapIds(
  org.Hs.eg.db,
  keys = as.character(cor_df$entrez),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

cor_df$symbol <- symbols

# Grupları tanımla
cor_df$group <- "Neither"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05] <- "Control_only"
cor_df$group[abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05] <- "HD_only"
cor_df$group[abs(cor_df$cor_ctrl) > 0.4 & cor_df$padj_ctrl < 0.05 &
               abs(cor_df$cor_hd)   > 0.4 & cor_df$padj_hd   < 0.05] <- "Both"

table(cor_df$group)

# Kaydet
dir.create("results/HD", recursive = TRUE, showWarnings = FALSE)
write.csv(cor_df, "results/H2AFZ_correlation_results.csv", row.names = FALSE)

# HD'ye özgü genler
hd_specific <- cor_df[cor_df$group == "HD_only", ]
hd_specific <- hd_specific[order(abs(hd_specific$cor_hd), decreasing = TRUE), ]
cat("\nHD'ye özgü top 10 gen:\n")
print(head(hd_specific[, c("symbol", "cor_ctrl", "cor_hd", "padj_hd")], 10))

# HD_only ve Both gruplarının pathway analizi için Entrez ID listesi hazırla
hd_cor_genes <- cor_df$entrez[cor_df$group %in% c("HD_only", "Both")]
cat("Pathway analizi için gen sayısı:", length(hd_cor_genes), "\n")

# HD-specific ve Both gruplarındaki tüm anlamlı genleri tablo olarak kaydet
cor_table <- cor_df[cor_df$group %in% c("HD_only", "Both"), ]

# Sembolü olmayan satırları filtrele
cor_table <- cor_table[!is.na(cor_table$symbol), ]

# Anlamlı kolonlar + sırala
cor_table <- cor_table[, c("symbol", "entrez", "cor_ctrl", "cor_hd",
                           "padj_ctrl", "padj_hd", "group")]
cor_table <- cor_table[order(abs(cor_table$cor_hd), decreasing = TRUE), ]

# Kaydet
write.table(cor_table,
            "huntington/results/H2AFZ_HD_correlated_genes.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Toplam gen:", nrow(cor_table), "\n")
head(cor_table, 15)


#----------------
#PATHWAY

library(clusterProfiler)

# GO analizi
go_results <- enrichGO(
  gene          = as.character(hd_cor_genes),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# KEGG analizi
kegg_results <- enrichKEGG(
  gene          = as.character(hd_cor_genes),
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("GO anlamlı pathway:", nrow(go_results@result[go_results@result$p.adjust < 0.05,]), "\n")
cat("KEGG anlamlı pathway:", nrow(kegg_results@result[kegg_results@result$p.adjust < 0.05,]), "\n")

# Top 10 GO
head(go_results@result[, c("Description", "GeneRatio", "p.adjust")], 10)

# Top 10 KEGG
head(kegg_results@result[, c("Description", "GeneRatio", "p.adjust")], 10)



#PLOOTTTTSSS

# OAS3 scatter
oas3_idx <- which(rownames(vst_mat) == "4940")
oas3_expr <- vst_mat[oas3_idx, ]

scatter_df <- data.frame(
  H2AFZ = as.numeric(h2afz_expr),
  OAS3  = as.numeric(oas3_expr),
  group = metadata$condition
)

p_scatter <- ggplot(scatter_df, aes(x = H2AFZ, y = OAS3, color = group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Control" = "#2980B9", "HD" = "#C0392B")) +
  theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom") +
  labs(title = "H2AFZ vs OAS3 Expression",
       subtitle = paste0("Control: ρ=", round(cor_df[cor_df$entrez=="4940","cor_ctrl"],2),
                         " | HD: ρ=", round(cor_df[cor_df$entrez=="4940","cor_hd"],2)),
       x = "H2AFZ (VST)", y = "OAS3 (VST)", color = NULL)

# Correlation comparison plot
plot_df <- cor_df[cor_df$group != "Neither", ]
group_colors <- c("Both"="#8E44AD", "Control_only"="#2980B9", "HD_only"="#C0392B")

p_cor <- ggplot(plot_df, aes(x = cor_ctrl, y = cor_hd, color = group)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_hline(yintercept = c(-0.4, 0.4), linetype="dashed", color="gray50", linewidth=0.4) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype="dashed", color="gray50", linewidth=0.4) +
  scale_color_manual(values = group_colors,
                     labels = c(paste0("Both (n=",         sum(cor_df$group=="Both"), ")"),
                                paste0("Control only (n=", sum(cor_df$group=="Control_only"), ")"),
                                paste0("HD only (n=",      sum(cor_df$group=="HD_only"), ")"))) +
  theme_classic() +
  theme(legend.position="bottom",
        plot.title=element_text(size=13, face="bold")) +
  labs(title="H2AFZ Correlation: Control vs HD",
       x="Spearman ρ (Control)", y="Spearman ρ (HD)", color=NULL)

# KEGG dotplot
kegg_df <- kegg_results@result[kegg_results@result$p.adjust < 0.05, ]
kegg_df$GeneRatio_num <- sapply(kegg_df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
kegg_top <- head(kegg_df[order(kegg_df$p.adjust), ], 15)
kegg_top$Description <- factor(kegg_top$Description, levels=rev(kegg_top$Description))

p_kegg <- ggplot(kegg_top, aes(x=GeneRatio_num, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  scale_color_gradient(low="#C0392B", high="#2980B9", name="p.adjust") +
  scale_size_continuous(range=c(3,10), name="Gene Count") +
  theme_classic() +
  theme(axis.text.y=element_text(size=10),
        plot.title=element_text(size=13, face="bold")) +
  labs(title="KEGG Pathway — H2AFZ Correlated Genes (HD)",
       x="Gene Ratio", y=NULL)

# GO dotplot
go_df <- go_results@result[go_results@result$p.adjust < 0.05, ]
go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
go_top <- head(go_df[order(go_df$p.adjust), ], 15)
go_top$Description <- factor(go_top$Description, levels=rev(go_top$Description))

p_go <- ggplot(go_top, aes(x=GeneRatio_num, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  scale_color_gradient(low="#C0392B", high="#2980B9", name="p.adjust") +
  scale_size_continuous(range=c(3,10), name="Gene Count") +
  theme_classic() +
  theme(axis.text.y=element_text(size=10),
        plot.title=element_text(size=13, face="bold")) +
  labs(title="GO Biological Process — H2AFZ Correlated Genes (HD)",
       x="Gene Ratio", y=NULL)

# Kaydet
ggsave("huntington/plots/scatter_H2AFZ_OAS3.png",        p_scatter, width=7,  height=6, dpi=300)
ggsave("huntington/plots/correlation_comparison.png",     p_cor,     width=8,  height=7, dpi=300)
ggsave("huntington/plots/kegg_dotplot.png",               p_kegg,    width=10, height=7, dpi=300)
ggsave("huntington/plots/go_dotplot.png",                 p_go,      width=11, height=7, dpi=300)


#KEGG Huntington disease pathway genlerini al
hd_pathway_genes <- kegg_results@result["hsa05016", "geneID"]
hd_genes_list <- strsplit(hd_pathway_genes, "/")[[1]]

# GO top pathway - mitochondrial respirasome assembly genlerini al
mito_pathway_genes <- go_results@result["GO:0097250", "geneID"]
mito_genes_list <- strsplit(mito_pathway_genes, "/")[[1]]

# HD-specific korelasyon genlerinden pathway'de olanları bul
# cor_df'te symbol kullan
hd_specific_genes <- cor_df[cor_df$group == "HD_only" & !is.na(cor_df$symbol), ]
hd_specific_genes <- hd_specific_genes[order(abs(hd_specific_genes$cor_hd), decreasing=TRUE), ]

# Hem HD-specific hem Huntington pathway'de olan
overlap_hd <- hd_specific_genes[hd_specific_genes$symbol %in% hd_genes_list, ]
cat("HD-specific + Huntington pathway:", nrow(overlap_hd), "\n")
print(head(overlap_hd[, c("symbol", "cor_ctrl", "cor_hd")], 10))

# Hem HD-specific hem mito pathway'de olan
overlap_mito <- hd_specific_genes[hd_specific_genes$symbol %in% mito_genes_list, ]
cat("HD-specific + Mito pathway:", nrow(overlap_mito), "\n")
print(head(overlap_mito[, c("symbol", "cor_ctrl", "cor_hd")], 10))


library(ggplot2)
library(gridExtra)

# 9 gen için Entrez ID'leri
target_genes <- overlap_mito[, c("symbol", "entrez", "cor_ctrl", "cor_hd")]

plot_list <- list()

for(i in 1:nrow(target_genes)) {

  gene_sym    <- target_genes$symbol[i]
  gene_entrez <- as.character(target_genes$entrez[i])
  rho_ctrl    <- round(target_genes$cor_ctrl[i], 2)
  rho_hd      <- round(target_genes$cor_hd[i], 2)

  gene_idx  <- which(rownames(vst_mat) == gene_entrez)
  gene_expr <- as.numeric(vst_mat[gene_idx, ])

  df <- data.frame(
    H2AFZ     = as.numeric(h2afz_expr),
    gene_expr = gene_expr,
    group     = metadata$condition
  )

  p <- ggplot(df, aes(x = H2AFZ, y = gene_expr, color = group)) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    scale_color_manual(values = c("Control" = "#2980B9", "HD" = "#C0392B")) +
    theme_classic() +
    theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "gray40"),
      axis.title    = element_text(size = 8),
      axis.text     = element_text(size = 7),
      legend.position = "none"
    ) +
    labs(
      title    = gene_sym,
      subtitle = paste0("Control: ρ=", rho_ctrl, " | HD: ρ=", rho_hd),
      x        = "H2AFZ (VST)",
      y        = paste0(gene_sym, " (VST)")
    )

  plot_list[[gene_sym]] <- p
}

# Legend için ayrı plot
legend_plot <- ggplot(data.frame(group = c("Control", "HD"), x = 1, y = 1),
                      aes(x = x, y = y, color = group)) +
  geom_point() +
  scale_color_manual(values = c("Control" = "#2980B9", "HD" = "#C0392B"),
                     name = "") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 11))

g_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
legend_grob <- g_legend(legend_plot)

# Title
title_grob <- grid::textGrob(
  "H2AFZ vs Mitochondrial Genes — HD-Specific Correlations",
  gp = grid::gpar(fontsize = 13, fontface = "bold")
)
subtitle_grob <- grid::textGrob(
  "Genes correlated with H2AFZ exclusively in HD (GO: mitochondrial respirasome assembly)",
  gp = grid::gpar(fontsize = 10, col = "gray40")
)

# 3x3 grid
final_plot <- gridExtra::arrangeGrob(
  title_grob,
  subtitle_grob,
  gridExtra::arrangeGrob(
    plot_list[[1]], plot_list[[2]], plot_list[[3]],
    plot_list[[4]], plot_list[[5]], plot_list[[6]],
    plot_list[[7]], plot_list[[8]], plot_list[[9]],
    ncol = 3
  ),
  legend_grob,
  heights = c(0.08, 0.05, 1, 0.08)
)

ggsave("huntington/plots/H2AFZ_mito_scatter_grid.png",
       final_plot, width = 12, height = 13, dpi = 300)
cat("Kaydedildi.\n")
