# =============================================================================
# UCEC — Histone Varyant GSEA Pathway Analizi
# Genler: H2AX (H2AFX), H2AZ1 (H2AFZ), MACROH2A1 (H2AFY), MACROH2A2 (H2AFY2)
# Yöntemler: GO-BP + KEGG + Reactome + Hallmark + C6 Oncogenic (GSEA) + GSVA
# Görselleştirme: Her varyant için ayrı dot plot + 4 varyant ortak bubble plot
#
# BRCA scriptinden farkı:
#   → Veri indirme yok: mevcut VST / dds / se_combat objeleri doğrudan yüklenir
#   → PATH_CONFIG bloğunda dosya yollarını kendi sisteminize göre ayarlayın
#   → UCEC'e özgü literatür referansları eklendi (bkz. aşağıdaki [L] blokları)
# =============================================================================
#
# LİTERATÜR BAĞLAMLARI:
#
#  [L1] H2AX / γH2AX — UCEC:
#       Yüksek H2AX ekspresyonu endometrial kanserde genomik instabilite ve
#       MSI (microsatellite instability) ile ilişkili; DNA onarım yolakları aktif.
#       Ref: TCGA UCEC kapsamlı genomik analizi (Cancer Cell 2013);
#            Chen et al. (2022) — γH2AX endometrial kanser prognoz belirteci.
#
#  [L2] H2AZ1 (H2A.Z) — UCEC:
#       H2A.Z promoter bölgelerinde aktif transkripsiyon için gerekli;
#       UCEC'de upregüle, proliferasyon ve hücre döngüsü pathway'leriyle ilişkili.
#       Ref: Li et al. (2024) s12964-024-01823 — RELA-HIF1A-EGFR ekseni;
#            Svotelis et al. (2010) PMID 20023423 — proliferasyon artışı.
#
#  [L3] MACROH2A1 (mH2A1) — UCEC — TEZ HİPOTEZİ:
#       MACROH2A1 UCEC'de anlamlı şekilde upregüle (DESeq2 p.adj < 0.05).
#       Tümör baskılayıcı ve onkojenik rol arasındaki bağlam bağımlı denge.
#       UCEC'de: EMT inhibisyonu, p53 yolağı, hormon sinyali ile ilişkili
#       pathway'lerde zenginleşme beklenmektedir.
#       Ref: Douet et al. (2022) PMC9016624; Gaspar-Maia (2023) PMC9950461;
#            Broggi et al. (2020) Frontiers Oncol.
#
#  [L4] MACROH2A2 (mH2A2) — UCEC — TEZ HİPOTEZİ:
#       MACROH2A2 UCEC'de DESeq2 ile non-significant ama Wilcoxon ile
#       anlamlı (imbalanced n: 35 normal vs 538 tümör); metodolojik dikkat.
#       Tümör baskılayıcı karakteri: enhancer baskılaması, onkojenik
#       program inhibisyonu. GSEA'da negatif NES beklentisi.
#       Ref: Filipescu et al. (2023) Nature Cell Biology;
#            Gaspar-Maia et al. (2023) PMC9950461.
#
#  [L5] UCEC'e özgü pathway beklentileri:
#       - PI3K/AKT/mTOR: UCEC'in en sık mutasyona uğrayan yolağı (%93 PIK3CA)
#       - POLE/MSI: hypermutator fenotip → DNA repair pathway zenginleşmesi
#       - Östrojen/progesteron sinyali: UCEC'de hormon bağımlı alt tipler
#       - p53 yolağı: serous UCEC alt tipinde sık mutasyon
#       Ref: TCGA UCEC (Cancer Cell 2013); Levine & Ellenson (2022) review.
#
# =============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(dplyr)
  library(ggplot2)
  library(matrixStats)
  library(ggrepel)
  library(RColorBrewer)
})

# =============================================================================
# PATH_CONFIG — Kendi sisteminize göre düzenleyin
# =============================================================================
# Aşağıdaki 4 satırdan sadece elinizde olan dosyanın yolunu doldurun.
# Script öncelik sırasına göre hangisi varsa onu kullanır:
#   Öncelik 1: vst_matrix_path  (en hızlı — direkt kullanır)
#   Öncelik 2: dds_path         (vst() yeniden hesaplar)
#   Öncelik 3: se_path          (DESeqDataSet + vst())
#   Öncelik 4: counts_path      (CSV/RDS — DESeqDataSet oluşturur + vst())

PATH_CONFIG <- list(
  vst_matrix_path = "results/UCEC/vst_counts.rds",
  dds_path        = "",
  se_path         = "",
  counts_path     = ""
)
# =============================================================================

cancer     <- "UCEC"
output_dir <- file.path("results", cancer, "gsea_pathway_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

TARGET_GENES <- list(
  H2AX      = "ENSG00000188486",   # H2AFX / γH2AX
  H2AZ1     = "ENSG00000164032",   # H2AFZ / H2A.Z
  MACROH2A1 = "ENSG00000134986",   # H2AFY  / mH2A1
  MACROH2A2 = "ENSG00000172264"    # H2AFY2 / mH2A2
)
histone_symbols <- names(TARGET_GENES)
ensembl_ids     <- unlist(TARGET_GENES)

VARIANT_COLORS <- c(
  H2AX      = "#F4A582",
  H2AZ1     = "#92C5DE",
  MACROH2A1 = "#4DAF4A",
  MACROH2A2 = "#F5C242"
)

cat("=== UCEC Histone Varyant GSEA Analizi ===\n")
cat("Hedef genler:", paste(histone_symbols, collapse = ", "), "\n\n")

# =============================================================================
# 1. VERİ YÜKLEME — Öncelik sırasına göre otomatik algılama
# =============================================================================
cat("─── Veri yükleniyor ───\n")

vst_matrix <- NULL

# ── Öncelik 1: VST matrisi ──────────────────────────────────────────────────
if (is.null(vst_matrix) && file.exists(PATH_CONFIG$vst_matrix_path)) {
  cat("VST matrisi yükleniyor:", PATH_CONFIG$vst_matrix_path, "\n")
  obj <- readRDS(PATH_CONFIG$vst_matrix_path)

  # DESeqTransform objesi mi yoksa düz matris mi?
  if (is(obj, "DESeqTransform") || is(obj, "SummarizedExperiment")) {
    vst_matrix <- assay(obj)
  } else if (is.matrix(obj)) {
    vst_matrix <- obj
  } else if (is.data.frame(obj)) {
    vst_matrix <- as.matrix(obj)
  }
  cat("  ✓ VST matrisi yüklendi:", nrow(vst_matrix), "gen ×",
      ncol(vst_matrix), "örnek\n\n")
}

# ── Öncelik 2: dds objesi ───────────────────────────────────────────────────
if (is.null(vst_matrix) && file.exists(PATH_CONFIG$dds_path)) {
  cat("dds objesi yükleniyor:", PATH_CONFIG$dds_path, "\n")
  dds <- readRDS(PATH_CONFIG$dds_path)

  # Sadece tümör örneklerini al (eğer karışık obje ise)
  if ("sample_type" %in% names(colData(dds))) {
    tumor_idx <- colData(dds)$sample_type %in%
      c("Primary Tumor", "Tumor", "01")
    if (sum(tumor_idx) > 0) dds <- dds[, tumor_idx]
    cat("  Tümör örnek sayısı:", ncol(dds), "\n")
  }

  # Düşük eksprese gen filtresi (zaten uygulanmışsa atlanır)
  if (max(rowSums(counts(dds))) > 10) {
    keep <- rowSums(counts(dds)) >= 10
    dds  <- dds[keep, ]
  }

  cat("  VST hesaplanıyor...\n")
  vst_obj    <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_obj)
  cat("  ✓ VST matrisi oluşturuldu:", nrow(vst_matrix), "gen ×",
      ncol(vst_matrix), "örnek\n\n")
}

# ── Öncelik 3: se_combat.RData ──────────────────────────────────────────────
if (is.null(vst_matrix) && file.exists(PATH_CONFIG$se_path)) {
  cat("SummarizedExperiment yükleniyor:", PATH_CONFIG$se_path, "\n")

  # .RData dosyası genellikle "se" adlı obje barındırır
  env <- new.env()
  load(PATH_CONFIG$se_path, envir = env)
  se_obj <- env[[ls(env)[1]]]   # ilk objeyi al

  # Tümör filtresi
  if ("sample_type" %in% names(colData(se_obj))) {
    tumor_idx <- colData(se_obj)$sample_type %in%
      c("Primary Tumor", "Tumor", "01")
    if (sum(tumor_idx) > 0) se_obj <- se_obj[, tumor_idx]
  }
  cat("  Örnek sayısı:", ncol(se_obj), "\n")

  dds <- DESeqDataSet(se_obj, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  keep <- rowSums(counts(dds)) >= 10
  dds  <- dds[keep, ]

  cat("  VST hesaplanıyor...\n")
  vst_obj    <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_obj)
  cat("  ✓ VST matrisi oluşturuldu:", nrow(vst_matrix), "gen ×",
      ncol(vst_matrix), "örnek\n\n")
}

# ── Öncelik 4: Normalize counts CSV/RDS ─────────────────────────────────────
if (is.null(vst_matrix) && file.exists(PATH_CONFIG$counts_path)) {
  cat("Normalize counts yükleniyor:", PATH_CONFIG$counts_path, "\n")

  ext <- tools::file_ext(PATH_CONFIG$counts_path)
  if (ext == "rds") {
    counts_mat <- readRDS(PATH_CONFIG$counts_path)
  } else {
    counts_mat <- read.csv(PATH_CONFIG$counts_path,
                           row.names = 1, check.names = FALSE)
  }
  counts_mat <- as.matrix(counts_mat)

  # Normalize counts'tan DESeq2 ile vst yapmak için integer gerekir.
  # Normalize counts zaten log-transformed ise direkt kullan.
  if (max(counts_mat, na.rm = TRUE) < 30) {
    # Muhtemelen log2-normalized — direkt VST yerine bu matrisi kullan
    cat("  Log-normalize matris algılandı (max =",
        round(max(counts_mat, na.rm = TRUE), 2), ") — direkt kullanılıyor.\n")
    vst_matrix <- counts_mat
  } else {
    # Ham normalize counts — vst() benzeri log dönüşüm uygula
    cat("  Count matrisi log2(x+1) dönüştürülüyor...\n")
    vst_matrix <- log2(counts_mat + 1)
  }
  cat("  ✓ Matris hazır:", nrow(vst_matrix), "gen ×",
      ncol(vst_matrix), "örnek\n\n")
}

# ── Hiçbiri bulunamazsa hata ver ────────────────────────────────────────────
if (is.null(vst_matrix)) {
  cat(
    "\n╔══════════════════════════════════════════════════════════════╗\n",
    "║  HATA: Hiçbir veri dosyası bulunamadı.                      ║\n",
    "║  PATH_CONFIG bloğundaki yolları kontrol edin.               ║\n",
    "║                                                              ║\n",
    "║  İpucu: Mevcut dosyanızın path'ini ilgili satıra yazın,     ║\n",
    "║  diğer satırları boş bırakabilirsiniz ('').                 ║\n",
    "╚══════════════════════════════════════════════════════════════╝\n"
  )
  stop("Veri bulunamadı — PATH_CONFIG'i güncelleyin.")
}

# ── Satır adı formatını kontrol et (versiyon numarasını kaldır) ─────────────
# TCGA Ensembl ID'leri bazen "ENSG00000188486.10" formatında gelir
if (any(grepl("\\.", rownames(vst_matrix)))) {
  cat("Not: Ensembl versiyon numaraları kaldırılıyor (ENSGXXX.N → ENSGXXX)\n")
  rownames(vst_matrix) <- gsub("\\..*", "", rownames(vst_matrix))
}

# ── Hedef genlerin varlığını kontrol et ─────────────────────────────────────
present <- intersect(ensembl_ids, rownames(vst_matrix))
missing <- setdiff(ensembl_ids, rownames(vst_matrix))

cat("Hedef gen durumu:\n")
for (g in histone_symbols) {
  eid   <- TARGET_GENES[[g]]
  found <- eid %in% rownames(vst_matrix)
  expr  <- if (found) round(mean(vst_matrix[eid, ]), 2) else NA
  cat(sprintf("  %-12s (%s): %s%s\n",
              g, eid,
              ifelse(found, "✓ BULUNDU", "✗ BULUNAMADI"),
              ifelse(found, paste0(" — ortalama expr = ", expr), "")))
}
cat("\n")

if (length(missing) > 0) {
  cat("UYARI:", length(missing), "gen matrisde yok, atlanacak.\n\n")
}

# =============================================================================
# 2. YARDIMCI FONKSİYONLAR
# =============================================================================
ensembl_to_entrez <- function(ensembl_vec) {
  mapIds(org.Hs.eg.db,
         keys      = gsub("\\..*", "", ensembl_vec),
         column    = "ENTREZID",
         keytype   = "ENSEMBL",
         multiVals = "first")
}

ensembl_to_symbol <- function(ensembl_vec) {
  mapIds(org.Hs.eg.db,
         keys      = gsub("\\..*", "", ensembl_vec),
         column    = "SYMBOL",
         keytype   = "ENSEMBL",
         multiVals = "first")
}

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

# =============================================================================
# 3. MSigDB GENE SET'LERİ HAZIRLA
# =============================================================================
cat("─── MSigDB gene set'leri yükleniyor ───\n")

# Hallmark
hallmark_entrez <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))
hallmark_list_entrez <- split(hallmark_entrez$entrez_gene, hallmark_entrez$gs_name)

hallmark_symbol      <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)
hallmark_list_symbol <- split(hallmark_symbol$gene_symbol, hallmark_symbol$gs_name)
cat("  Hallmark:", length(hallmark_list_entrez), "set\n")

# KEGG — versiyon-bağımsız subcategory tespiti
all_collections <- msigdbr_collections()
sub_col <- if ("gs_subcollection" %in% names(all_collections)) "gs_subcollection" else "gs_subcat"
cat_col  <- if ("gs_collection"   %in% names(all_collections)) "gs_collection"   else "gs_cat"

kegg_sub_all <- all_collections[[sub_col]][
  grepl("KEGG", all_collections[[sub_col]], ignore.case = TRUE) &
    all_collections[[cat_col]] == "C2"
]
kegg_sub_use <- if ("CP:KEGG_LEGACY" %in% kegg_sub_all) "CP:KEGG_LEGACY" else kegg_sub_all[1]
cat("  KEGG subcategory:", kegg_sub_use, "\n")

kegg_sets <- msigdbr(species = "Homo sapiens",
                     category = "C2", subcategory = kegg_sub_use) %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))
cat("  KEGG:", length(unique(kegg_sets$gs_name)), "set\n")

# Reactome
reactome_sets <- msigdbr(species = "Homo sapiens",
                         category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))
cat("  Reactome:", length(unique(reactome_sets$gs_name)), "set\n")

# C6 Oncogenic Signatures
# UCEC için özellikle kritik: PI3K/AKT/mTOR (%93 sıklık), KRAS, p53 [L5]
oncogenic_sets <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))
cat("  C6 Oncogenic:", length(unique(oncogenic_sets$gs_name)), "set\n\n")

# =============================================================================
# 4. GSEA ANALİZİ — Her varyant için GO-BP + KEGG + Reactome + Hallmark + C6
# =============================================================================
# Strateji: Korelasyon bazlı sıralı liste (Spearman r)
# Tüm örneklerdeki ekspresyon değerleri kullanılır — sadece tümör.
# UCEC için bağlam: MACROH2A1/2 paralog divergence hipotezi [L3, L4]
# =============================================================================
cat("─── GSEA analizleri başlıyor ───\n\n")

all_gsea_results <- list()

for (h_name in histone_symbols) {

  ens_id <- TARGET_GENES[[h_name]]

  if (!ens_id %in% rownames(vst_matrix)) {
    cat("ATLANADI:", h_name, "(matrisde yok)\n\n"); next
  }

  cat("═══════════════════════════════════\n")
  cat("GSEA:", h_name, "\n")
  cat("═══════════════════════════════════\n")

  # Spearman korelasyon — hedef gen vs. tüm genler
  target_expr  <- vst_matrix[ens_id, ]
  correlations <- apply(vst_matrix, 1, function(g) {
    tryCatch(
      cor(target_expr, g, method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    )
  })

  gene_list <- correlations[!is.na(correlations)]
  gene_list <- gene_list[abs(gene_list) < 0.9999]  # kendi korelasyonunu çıkar
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Ensembl → Entrez, duplikat temizle
  entrez_map <- ensembl_to_entrez(names(gene_list))
  gene_df <- data.frame(
    ensembl     = names(gene_list),
    entrez      = entrez_map,
    correlation = gene_list,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(entrez)) %>%
    group_by(entrez) %>%
    slice_max(abs(correlation), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene_list_final <- setNames(gene_df$correlation, gene_df$entrez)
  gene_list_final <- sort(gene_list_final, decreasing = TRUE)
  cat("  Korelasyon listesi:", length(gene_list_final), "gen\n")

  # Korelasyon dağılımı özeti (MACROH2A1/2 için paralog divergence kontrolü)
  cat(sprintf("  Korelasyon özeti: min=%.3f, Q25=%.3f, median=%.3f, Q75=%.3f, max=%.3f\n",
              min(gene_list_final), quantile(gene_list_final, 0.25),
              median(gene_list_final), quantile(gene_list_final, 0.75),
              max(gene_list_final)))

  # ── GO Biological Process ──────────────────────────────────────────────────
  cat("  GO-BP...\n")
  gsea_go <- tryCatch(
    gseGO(geneList      = gene_list_final,
          OrgDb         = org.Hs.eg.db,
          ont           = "BP",
          minGSSize     = 15,
          maxGSSize     = 500,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          verbose       = FALSE,
          seed          = 42),
    error = function(e) { cat("    HATA:", conditionMessage(e), "\n"); NULL }
  )
  n_go <- if (!is.null(gsea_go)) nrow(gsea_go@result) else 0
  cat("    GO-BP:", n_go, "pathway\n")

  # ── KEGG ──────────────────────────────────────────────────────────────────
  cat("  KEGG...\n")
  gsea_kegg <- tryCatch(
    GSEA(geneList      = gene_list_final,
         TERM2GENE     = kegg_sets %>% rename(term = gs_name, gene = entrez_gene),
         minGSSize     = 15,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         verbose       = FALSE,
         seed          = 42),
    error = function(e) { cat("    HATA:", conditionMessage(e), "\n"); NULL }
  )
  n_kegg <- if (!is.null(gsea_kegg)) nrow(gsea_kegg@result) else 0
  cat("    KEGG:", n_kegg, "pathway\n")

  # ── Reactome ──────────────────────────────────────────────────────────────
  cat("  Reactome...\n")
  gsea_reactome <- tryCatch(
    GSEA(geneList      = gene_list_final,
         TERM2GENE     = reactome_sets %>% rename(term = gs_name, gene = entrez_gene),
         minGSSize     = 15,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         verbose       = FALSE,
         seed          = 42),
    error = function(e) { cat("    HATA:", conditionMessage(e), "\n"); NULL }
  )
  n_reactome <- if (!is.null(gsea_reactome)) nrow(gsea_reactome@result) else 0
  cat("    Reactome:", n_reactome, "pathway\n")

  # ── Hallmark ──────────────────────────────────────────────────────────────
  cat("  Hallmark...\n")
  gsea_hallmark <- tryCatch(
    GSEA(geneList      = gene_list_final,
         TERM2GENE     = hallmark_entrez %>% rename(term = gs_name, gene = entrez_gene),
         minGSSize     = 15,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         verbose       = FALSE,
         seed          = 42),
    error = function(e) { cat("    HATA:", conditionMessage(e), "\n"); NULL }
  )
  n_hallmark <- if (!is.null(gsea_hallmark)) nrow(gsea_hallmark@result) else 0
  cat("    Hallmark:", n_hallmark, "pathway\n")

  # ── C6 Oncogenic — UCEC için PI3K/AKT, KRAS, p53 sinyalleri kritik [L5] ──
  cat("  C6 Oncogenic...\n")
  gsea_oncogenic <- tryCatch(
    GSEA(geneList      = gene_list_final,
         TERM2GENE     = oncogenic_sets %>% rename(term = gs_name, gene = entrez_gene),
         minGSSize     = 10,
         maxGSSize     = 500,
         pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         verbose       = FALSE,
         seed          = 42),
    error = function(e) { cat("    HATA:", conditionMessage(e), "\n"); NULL }
  )
  n_oncogenic <- if (!is.null(gsea_oncogenic)) nrow(gsea_oncogenic@result) else 0
  cat("    C6 Oncogenic:", n_oncogenic, "imza\n\n")

  # CSV kaydet
  save_gsea <- function(obj, n, tag) {
    if (!is.null(obj) && n > 0)
      write.csv(obj@result,
                file.path(output_dir, paste0(h_name, "_GSEA_", tag, ".csv")),
                row.names = FALSE)
  }
  save_gsea(gsea_go,        n_go,        "GO-BP")
  save_gsea(gsea_kegg,      n_kegg,      "KEGG")
  save_gsea(gsea_reactome,  n_reactome,  "Reactome")
  save_gsea(gsea_hallmark,  n_hallmark,  "Hallmark")
  save_gsea(gsea_oncogenic, n_oncogenic, "C6_Oncogenic")

  all_gsea_results[[h_name]] <- list(
    GO        = gsea_go,
    KEGG      = gsea_kegg,
    Reactome  = gsea_reactome,
    Hallmark  = gsea_hallmark,
    Oncogenic = gsea_oncogenic,
    gene_list = gene_list_final
  )
}

saveRDS(all_gsea_results,
        file.path(output_dir, "all_gsea_results_UCEC.rds"))
cat("✓ GSEA tamamlandı.\n\n")

# =============================================================================
# 5. GSVA — Hallmark pathway aktivasyon skorları
# =============================================================================
cat("─── GSVA Hallmark analizi ───\n")

row_symbols <- ensembl_to_symbol(rownames(vst_matrix))
symbol_df   <- data.frame(
  ensembl = rownames(vst_matrix),
  symbol  = row_symbols,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice(1) %>%
  ungroup()

vst_symbol <- vst_matrix[symbol_df$ensembl, ]
rownames(vst_symbol) <- symbol_df$symbol
cat("  GSVA matrisi:", nrow(vst_symbol), "gen ×", ncol(vst_symbol), "örnek\n")

gsva_param  <- gsvaParam(vst_symbol, hallmark_list_symbol)
gsva_scores <- gsva(gsva_param, verbose = FALSE)

write.csv(gsva_scores,
          file.path(output_dir, "GSVA_Hallmark_scores_UCEC.csv"))
saveRDS(gsva_scores,
        file.path(output_dir, "gsva_hallmark_scores_UCEC.rds"))
cat("  ✓ GSVA:", nrow(gsva_scores), "pathway ×", ncol(gsva_scores), "örnek\n\n")

# =============================================================================
# 6. GÖRSELLEŞTİRME
# =============================================================================
gsea_colors <- c(Activated = "#B2182B", Suppressed = "#2166AC")

# ── 6A. Her varyant için ayrı Hallmark dot plot ────────────────────────────
cat("─── Bireysel Hallmark dot plotlar ───\n")

for (h_name in histone_symbols) {

  obj <- all_gsea_results[[h_name]]$Hallmark
  if (is.null(obj) || nrow(obj@result) == 0) {
    cat("  Anlamlı sonuç yok:", h_name, "\n"); next
  }

  dot_df <- obj@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(
      Direction = ifelse(NES > 0, "Activated", "Suppressed"),
      log10padj = -log10(p.adjust),
      Pathway   = clean_pathway(ID),
      Pathway   = ifelse(nchar(Pathway) > 42,
                         paste0(substr(Pathway, 1, 40), ".."), Pathway),
      Pathway   = factor(Pathway, levels = Pathway[order(NES)])
    ) %>%
    arrange(NES)

  if (nrow(dot_df) == 0) next

  p_dot <- ggplot(dot_df,
                  aes(x = NES, y = Pathway, size = setSize, color = Direction)) +
    geom_point(alpha = 0.88) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.4) +
    scale_color_manual(values = gsea_colors) +
    scale_size_continuous(range = c(3, 9), name = "Gene Set Size",
                          guide = guide_legend(override.aes = list(color = "grey40"))) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title    = paste0("TCGA-UCEC — ", h_name, ": Hallmark GSEA"),
      subtitle = paste0("Spearman korelasyon bazlı | adj.p < 0.05, n = ",
                        nrow(dot_df)),
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 9),
      legend.position    = "right",
      panel.grid.major.y = element_line(color = "grey93"),
      panel.grid.minor   = element_blank()
    )

  out <- file.path(output_dir, paste0(h_name, "_GSEA_Hallmark_dotplot_UCEC.jpeg"))
  ggsave(out, p_dot,
         width  = 9,
         height = max(4, nrow(dot_df) * 0.38 + 2.5),
         dpi = 220, units = "in", limitsize = FALSE)
  cat("  →", basename(out), "\n")
}

# ── 6B. Ortak Hallmark bubble dot plot ────────────────────────────────────
cat("\n─── Ortak Hallmark bubble dot plot ───\n")

combined_df <- bind_rows(lapply(histone_symbols, function(h) {
  obj <- all_gsea_results[[h]]$Hallmark
  if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
  obj@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Gene = h, log10padj = -log10(p.adjust),
           Pathway = clean_pathway(ID)) %>%
    select(Pathway, NES, log10padj, setSize, Gene, p.adjust)
}))

if (!is.null(combined_df) && nrow(combined_df) > 0) {

  TOP_N_HALLMARK <- 30

  pathway_counts <- combined_df %>%
    group_by(Pathway) %>%
    summarise(n_variants = n_distinct(Gene), .groups = "drop")

  common_pathways <- pathway_counts %>%
    filter(n_variants >= 2) %>%
    left_join(
      combined_df %>%
        group_by(Pathway) %>%
        summarise(mean_abs_NES = mean(abs(NES)), .groups = "drop"),
      by = "Pathway"
    ) %>%
    arrange(desc(n_variants), desc(mean_abs_NES)) %>%
    slice_head(n = TOP_N_HALLMARK)

  cat("  Gösterilecek pathway:", nrow(common_pathways),
      "(>=2 varyant toplam:", nrow(pathway_counts %>% filter(n_variants >= 2)), ")\n")

  plot_df <- combined_df %>%
    filter(Pathway %in% common_pathways$Pathway)

  pathway_order <- plot_df %>%
    group_by(Pathway) %>%
    summarise(mean_NES = mean(NES), .groups = "drop") %>%
    arrange(mean_NES) %>% pull(Pathway)

  plot_df <- plot_df %>%
    mutate(Pathway = factor(Pathway, levels = pathway_order),
           Gene    = factor(Gene, levels = histone_symbols))

  n_paths <- length(pathway_order)
  nes_lim <- max(abs(plot_df$NES)) * 1.15

  p_common <- ggplot(plot_df,
                     aes(x = NES, y = Pathway, size = log10padj, color = Gene)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.45) +
    geom_point(alpha = 0.82, stroke = 0.3) +
    scale_color_manual(values = VARIANT_COLORS, name = "Histone Variant") +
    scale_size_continuous(range = c(2, 8),
                          name  = expression(-log[10](padj)),
                          breaks = c(1, 2, 3, 4)) +
    scale_x_continuous(limits = c(-nes_lim, nes_lim),
                       expand = expansion(mult = 0.05)) +
    labs(
      title    = "Common GSEA Pathways Across 4 Histone Variants",
      subtitle = paste0("TCGA-UCEC — top ", n_paths,
                        " pathways enriched in >=2 genes"),
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    guides(color = guide_legend(title = "Histone Variant",
                                override.aes = list(size = 4)),
           size  = guide_legend(title = expression(-log[10](padj)))) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 8.5),
      axis.text.x        = element_text(size = 9),
      axis.title.x       = element_text(size = 10),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8.5),
      panel.grid.major.y = element_line(color = "grey93"),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA)
    )

  out_h <- file.path(output_dir,
                     "UCEC_CommonPathways_4variants_bubble_dotplot.jpeg")
  ggsave(out_h, p_common,
         width = 9, height = n_paths * 0.28 + 2.8,
         dpi = 240, units = "in", limitsize = FALSE)
  cat("  →", basename(out_h), "\n")
}

# ── 6C. Ortak GO-BP bubble dot plot ────────────────────────────────────────
cat("\n─── Ortak GO-BP bubble dot plot ───\n")

combined_go_df <- bind_rows(lapply(histone_symbols, function(h) {
  obj <- all_gsea_results[[h]]$GO
  if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
  obj@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Gene = h, log10padj = -log10(p.adjust), Pathway = Description) %>%
    select(Pathway, NES, log10padj, setSize, Gene, p.adjust)
}))

if (!is.null(combined_go_df) && nrow(combined_go_df) > 0) {

  TOP_N_GO <- 25

  go_counts <- combined_go_df %>%
    group_by(Pathway) %>%
    summarise(n_variants = n_distinct(Gene), .groups = "drop")

  common_go <- go_counts %>%
    filter(n_variants >= 2) %>%
    left_join(
      combined_go_df %>%
        group_by(Pathway) %>%
        summarise(mean_abs_NES = mean(abs(NES)), .groups = "drop"),
      by = "Pathway"
    ) %>%
    arrange(desc(n_variants), desc(mean_abs_NES)) %>%
    slice_head(n = TOP_N_GO)

  cat("  GO-BP gösterilecek pathway:", nrow(common_go), "\n")

  go_plot_df <- combined_go_df %>%
    filter(Pathway %in% common_go$Pathway) %>%
    mutate(
      Pathway = ifelse(nchar(Pathway) > 48,
                       paste0(substr(Pathway, 1, 46), ".."), Pathway)
    )

  go_order <- go_plot_df %>%
    group_by(Pathway) %>%
    summarise(mean_NES = mean(NES), .groups = "drop") %>%
    arrange(mean_NES) %>% pull(Pathway)

  go_plot_df <- go_plot_df %>%
    mutate(Pathway = factor(Pathway, levels = go_order),
           Gene    = factor(Gene, levels = histone_symbols))

  n_go_paths  <- length(go_order)
  nes_lim_go  <- max(abs(go_plot_df$NES)) * 1.15

  p_go <- ggplot(go_plot_df,
                 aes(x = NES, y = Pathway, size = log10padj, color = Gene)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey55", linewidth = 0.45) +
    geom_point(alpha = 0.82, stroke = 0.3) +
    scale_color_manual(values = VARIANT_COLORS, name = "Histone Variant") +
    scale_size_continuous(range = c(2, 8),
                          name  = expression(-log[10](padj))) +
    scale_x_continuous(limits = c(-nes_lim_go, nes_lim_go),
                       expand = expansion(mult = 0.05)) +
    labs(
      title    = "Common GO-BP GSEA Pathways Across 4 Histone Variants",
      subtitle = paste0("TCGA-UCEC — top ", n_go_paths,
                        " pathways enriched in >=2 genes"),
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    guides(color = guide_legend(title = "Histone Variant",
                                override.aes = list(size = 4)),
           size  = guide_legend(title = expression(-log[10](padj)))) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 8.5),
      axis.text.x        = element_text(size = 9),
      axis.title.x       = element_text(size = 10),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 8.5),
      panel.grid.major.y = element_line(color = "grey93"),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA)
    )

  out_go <- file.path(output_dir,
                      "UCEC_CommonGO-BP_4variants_bubble_dotplot.jpeg")
  ggsave(out_go, p_go,
         width = 9, height = n_go_paths * 0.28 + 2.8,
         dpi = 240, units = "in", limitsize = FALSE)
  cat("  →", basename(out_go), "\n")
}

# ── 6D. GSVA bar plot — Her varyant için High vs Low ──────────────────────
cat("\n─── GSVA bar plotlar ───\n")

for (h_name in histone_symbols) {

  if (!h_name %in% rownames(vst_symbol)) {
    cat("  ATLANADI:", h_name, "\n"); next
  }

  h_expr   <- vst_symbol[h_name, ]
  high_idx <- which(h_expr > quantile(h_expr, 0.66))
  low_idx  <- which(h_expr < quantile(h_expr, 0.33))

  med_high <- rowMedians(gsva_scores[, high_idx])
  med_low  <- rowMedians(gsva_scores[, low_idx])
  names(med_high) <- rownames(gsva_scores)
  names(med_low)  <- rownames(gsva_scores)

  delta <- med_high - med_low
  top20 <- names(sort(abs(delta), decreasing = TRUE)[1:min(20, length(delta))])

  bar_df <- data.frame(
    Pathway    = rep(clean_pathway(top20), 2),
    Group      = rep(c(paste0(h_name, " High"), paste0(h_name, " Low")),
                     each = length(top20)),
    MedianGSVA = c(med_high[top20], med_low[top20])
  ) %>%
    mutate(
      Pathway = factor(Pathway,
                       levels = clean_pathway(top20[order(delta[top20])])),
      Group   = factor(Group,
                       levels = c(paste0(h_name, " High"), paste0(h_name, " Low")))
    )

  p_bar <- ggplot(bar_df,
                  aes(x = MedianGSVA, y = Pathway, fill = Group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.64, alpha = 0.90) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(
      values = setNames(c("#B2182B", "#4393C3"),
                        c(paste0(h_name, " High"), paste0(h_name, " Low"))),
      name = paste0(h_name, " Expression")
    ) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(
      title    = paste0("TCGA-UCEC — ", h_name, ": Hallmark GSVA Scores"),
      subtitle = "Median GSVA — High vs Low expression (top 20 by differential score)",
      x = "Median GSVA Score", y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(color = "grey40", size = 9),
      axis.text.y        = element_text(size = 9),
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.minor   = element_blank()
    )

  out_bar <- file.path(output_dir,
                       paste0(h_name, "_GSVA_Hallmark_barplot_UCEC.jpeg"))
  ggsave(out_bar, p_bar, width = 9, height = 7, dpi = 220, units = "in")
  cat("  →", basename(out_bar), "\n")
}

# =============================================================================
# 7. ÖZET TABLO
# =============================================================================
cat("\n═══════════════════════════════════════════════════\n")
cat("ÖZET — UCEC GSEA Sonuç Sayıları\n")
cat("═══════════════════════════════════════════════════\n")

summary_df <- bind_rows(lapply(histone_symbols, function(h) {
  res <- all_gsea_results[[h]]
  if (is.null(res)) return(NULL)
  data.frame(
    Gene         = h,
    GO_BP        = if (!is.null(res$GO))        nrow(res$GO@result)        else 0,
    KEGG         = if (!is.null(res$KEGG))       nrow(res$KEGG@result)      else 0,
    Reactome     = if (!is.null(res$Reactome))   nrow(res$Reactome@result)  else 0,
    Hallmark     = if (!is.null(res$Hallmark))   nrow(res$Hallmark@result)  else 0,
    C6_Oncogenic = if (!is.null(res$Oncogenic))  nrow(res$Oncogenic@result) else 0
  )
}))

print(summary_df)
write.csv(summary_df,
          file.path(output_dir, "UCEC_GSEA_summary.csv"),
          row.names = FALSE)

cat("\nÇıktı dizini:", output_dir, "\n")
cat("✓ Tüm analizler tamamlandı.\n")
