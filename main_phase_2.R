
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(clusterProfiler)

source("lib/filter_correlated_genes.R")

config <- yaml::yaml.load_file("config.yaml")
cancer <- "LUAD"

plots_dir <- "results/LUAD/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#----------------------------------------------
#Correlation Analysis

#step1: DGCA (differential) - Tumor specific findings - needs to be corrected
# library(DGCA)
# dgca_results <- run_dgca_analysis(
#   counts = vst_counts,
#   groups = sample_type,
#   histone_genes = ensembl_ids
# )

#step2: combined correlation
correlation_results <- run_correlation_analysis(
  dds = deseq2_results$dds,
  histone_ensembl_ids = ensembl_ids,
  correlation_threshold = 0.5,
  padj_threshold = 0.001,
  variance_percentile = 0.5,
  use_vst = TRUE
)

# Correlation results
saveRDS(correlation_results,
        paste0(results_dir, "/correlation_results.rds"))
write.csv(correlated_genes,
          paste0(results_dir, "/correlated_genes_summary.csv"),
          row.names = FALSE)

# Extract unique correlated genes
correlated_genes <- extract_histone_correlated_genes(correlation_results)
gene_list <- correlated_genes$gene_id

#not to run
# correlated <- read.csv("results/LUAD/plots/correlated_genes_summary.csv")
# gene_list <- correlated$gene_id

# ENSEMBL → ENTREZ
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = gene_list,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

entrez_ids <- entrez_ids[!is.na(entrez_ids)]
print(paste("Mapped to ENTREZ:", length(entrez_ids)))

# GO enrichment
go_results <- enrichGO(
  gene = entrez_ids,           # 1,272 ENTREZ ID
  OrgDb = org.Hs.eg.db,        # Human database
  ont = "BP",                   # Biological Process
  pAdjustMethod = "BH",        # FDR correction
  pvalueCutoff = 0.05,         # Significance threshold
  qvalueCutoff = 0.2           # FDR threshold
)

#how many pathways
dim(go_results)
go_df <- as.data.frame(go_results)
head(go_df)

#plots
#dot plot
p <- dotplot(
  go_results,
  showCategory = 20,
  font.size = 12,
  title = "GO Enrichment - LUAD Histone-Correlated Genes"
) +
  theme(axis.text.y = element_text(size = 10))

ggsave("results/LUAD/plots/go_dotplot_custom.png", p, width = 12, height = 10)

#bar plot
pdf("results/LUAD/plots/pathway_barplot.pdf", width = 10, height = 8)
barplot(go_results, showCategory = 15)
dev.off()

#gene concept network
pdf("results/LUAD/plots/pathway_cnetplot.pdf", width = 14, height = 12)
cnetplot(go_results, showCategory = 5)
dev.off()

#--------------------------------------------
# KEGG enrichment
kegg_results <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",  # human
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

#dot plot
dotplot(kegg_results,
        showCategory = 20,
        font.size = 12,
        title = "KEGG Pathway Enrichment: LUAD Histone-Correlated Genes") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold")
  ) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Gene Ratio",
    y = "KEGG Pathway",
    color = "p.adjust",
    size = "Count"
  )

ggsave("results/LUAD/plots/kegg_dotplot.png", width = 12, height = 10, dpi = 300)

# bar plot
barplot(kegg_results,
        showCategory = 15,
        title = "Top 15 KEGG Pathways") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 11)
  )

ggsave("results/LUAD/plots/kegg_barplot.png", width = 12, height = 10, dpi = 300)

#emapp
kegg_results <- pairwise_termsim(kegg_results)

emapplot(kegg_results,
         showCategory = 20,
         layout = "kk") +
  theme(legend.position = "right")

ggsave("results/LUAD/kegg_emapp.png", width = 12, height = 10, dpi = 300)

#-------------------------------------------------
#GSEA

de_results <- read.csv("results/LUAD/plots/deseq2_full_results.csv")
nrow(de_results)
head(de_results)

de_results$padj_clean <- de_results$padj
de_results$padj_clean[de_results$padj == 0] <- 1e-300

de_clean <- de_results[!is.na(de_results$log2FoldChange) &
                         !is.na(de_results$padj_clean), ]

de_clean$rank_score <- sign(de_clean$log2FoldChange) *
  (-log10(de_clean$padj_clean))

summary(de_clean$rank_score)
head(de_clean[order(-de_clean$rank_score),
              c("ensembl_id", "log2FoldChange", "padj", "rank_score")], 10)

entrez_all <- mapIds(
  org.Hs.eg.db,
  keys = de_clean$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

de_clean$entrez <- entrez_all
de_clean <- de_clean[!is.na(de_clean$entrez), ]

gene_list_gsea <- setNames(de_clean$rank_score, de_clean$entrez)
gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)

#GSEA-GO
gsea_go <- gseGO(
  geneList = gene_list_gsea,      # Ranked list
  OrgDb = org.Hs.eg.db,            # Human database
  ont = "BP",                       # Biological Process
  pAdjustMethod = "BH",            # FDR correction
  pvalueCutoff = 0.05,             # Significance
  verbose = TRUE                    # Progress
)
nrow(gsea_go) #how many pathways

gsea_go_df <- as.data.frame(gsea_go)
gsea_go_df[order(-gsea_go_df$NES), c("Description", "NES", "pvalue", "p.adjust")][1:10, ]

top_up <- gsea_go_df[gsea_go_df$NES > 0, ][1:10, ]
top_down <- gsea_go_df[gsea_go_df$NES < 0, ][1:10, ]
plot_data <- rbind(top_up, top_down)

# Barplot
ggplot(plot_data, aes(x = reorder(Description, NES), y = NES, fill = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "GSEA GO: LUAD -Top 10 Up/Down Regulated Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "NES"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  )

ggsave("results/LUAD/plots/gsea_go_barplot.png", width = 12, height = 8, dpi = 300)

gseaplot2(gsea_go,
          geneSetID = 1:3,  # İlk 3 pathway
          pvalue_table = TRUE,
          ES_geom = "line",
          title = "Top 3 Upregulated Pathways")

ggsave("results/LUAD/gsea_enrichment_plot.png", width = 12, height = 10)

#ridge plot
ridgeplot(gsea_go,
          showCategory = 30,
          fill = "p.adjust",
          core_enrichment = TRUE) +
  theme(axis.text.y = element_text(size = 10))

ggsave("results/LUAD/plots/gsea_ridgeplot.png", width = 14, height = 12)

#heatmap plot
heatplot(gsea_go,
         showCategory = 15,  # Daha az pathway
         foldChange = gene_list_gsea) +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

ggsave("results/LUAD/plots/gsea_heatmap_fixed.png", width = 16, height = 10)

#network plot
gsea_go <- pairwise_termsim(gsea_go) #first calculate similarity

emapplot(gsea_go,
         showCategory = 30,
         layout = "nicely",
         cex_label_category = 0.6) +
  theme(legend.position = "right")

ggsave("results/LUAD/plots/gsea_network.png", width = 14, height = 12)

#volcano plot
ggplot(gsea_go_df, aes(x = NES, y = -log10(p.adjust))) +
  geom_point(aes(color = NES, size = setSize), alpha = 0.6) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "GSEA GO Pathway Volcano Plot",
    x = "Normalized Enrichment Score (NES)",
    y = "-log10(p.adjust)",
    color = "NES",
    size = "Gene Set Size"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("results/LUAD/plots/gsea_volcano.png", width = 12, height = 8, dpi = 300)

emapplot(gsea_go,
         showCategory = 30,
         layout = "nicely",  # veya "kk" veya "fr"
         cex.params = list(category_label = 0.5)) +
  theme(legend.position = "right")

ggsave("results/LUAD/plots/gsea_network_layout.png", width = 14, height = 12)

#GSEA-KEGG
gsea_kegg <- gseKEGG(
  geneList = gene_list_gsea,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = TRUE
)

gsea_kegg_df <- as.data.frame(gsea_kegg)

dotplot(gsea_kegg,
        showCategory = 20,
        split = ".sign",
        font.size = 10) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA KEGG Pathways") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/LUAD/plots/gsea_kegg_dotplot.pdf", width = 14, height = 10)

gseaplot2(gsea_kegg,
          geneSetID = 1:3,
          pvalue_table = TRUE,
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),  # Panel heights
          subplots = 1:3,
          pvalue_table_args = list(
            x = 0.7,   # Sağa kaydır (0-1 arası)
            y = 0.95   # Yukarı kaydır (0-1 arası)
          ))

ggsave("results/LUAD/plots/gsea_kegg_enrichment_fixed.pdf", width = 12, height = 8)

p <- gseaplot2(gsea_kegg,
               geneSetID = 1,
               title = "Cell cycle (NES=2.3, Upregulated)",
               color = "red",
               pvalue_table = TRUE,
               ES_geom = "line")


ggsave("results/LUAD/plots/gsea_kegg_cell_cycle.pdf", p, width = 10, height = 6)

cgmp_index <- which(grepl("cGMP-PKG", gsea_kegg_df$Description))

p_2 <- gseaplot2(gsea_kegg,
               geneSetID = cgmp_index,
               title = "cGMP-PKG signaling pathway (NES=-1.80, Downregulated)",
               pvalue_table = TRUE,
               color = "blue")  # Downregulated = blue

ggsave("results/LUAD/plots/gsea_kegg_cgmp_pkg.png", p_2, width = 10, height = 6)

#pathway control
gsea_kegg_df[1:10, c("ID", "Description", "NES", "p.adjust")]
gsea_kegg_df[2, "Description"]

#---------------------------------------------------
#PPI -Protein Protein Interaction - STRINGdb
install.packages("igraph")
install.packages("ggraph")
install.packages("visNetwork")
BiocManager::install("STRINGdb")

library(STRINGdb)
library(igraph)
