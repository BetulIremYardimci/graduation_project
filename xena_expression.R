library(data.table)

pheno <- fread("data/xena/TcgaTargetGTEX_phenotype.txt.gz", sep = "\t")

# Sütun isimlerini gör
print(colnames(pheno))


# GTEx tissue tiplerini listele
print(sort(unique(pheno[`_study` == "GTEX", `_primary_site`])))


# GTEx Ovary örneklerini filtrele
gtex_ovary <- pheno[`_study` == "GTEX" & `_primary_site` == "Ovary"]
message("GTEx Ovary normal örnek sayısı: ", nrow(gtex_ovary))
print(head(gtex_ovary))

# Sample ID'leri al
gtex_ids <- gtex_ovary$sample
message("Örnek IDs (ilk 5): ", paste(head(gtex_ids, 5), collapse = ", "))

# Toplam örnek sayısı
message("Toplam GTEx Ovary örnek sayısı: ", nrow(gtex_ovary))

# Count matrisinin header'ını oku (hangi kolonlar var?)
message("Count matrisi header okunuyor...")
header <- fread(
  "data/xena/TcgaTargetGtex_gene_expected_count.gz",
  nrows = 0,
  sep = "\t"
)
all_cols <- colnames(header)
message("Toplam kolon sayısı: ", length(all_cols))

# GTEx ovary ID'lerinin kaçı matris'te var?
found_ids <- intersect(all_cols, gtex_ids)
message("Matris'te bulunan GTEx Ovary örnek sayısı: ", length(found_ids))

# Sadece bu kolonları yükle (çok daha hızlı)
message("Sadece GTEx Ovary kolonları yükleniyor...")
gtex_counts_raw <- fread(
  "data/xena/TcgaTargetGtex_gene_expected_count.gz",
  select = c("sample", found_ids),
  sep = "\t"
)

message("Boyut: ", nrow(gtex_counts_raw), " gen x ", ncol(gtex_counts_raw) - 1, " örnek")


# Sadece GTEx Ovary kolonlarını yükle
message("Count matrisi yükleniyor (birkaç dakika sürebilir)...")
gtex_counts_raw <- fread(
  "data/xena/TcgaTargetGtex_gene_expected_count.gz",
  select = c("sample", found_ids),
  sep = "\t"
)

message("Ham matris boyutu: ", nrow(gtex_counts_raw), " x ", ncol(gtex_counts_raw))

# Gene ID'yi rowname yap
gene_ids <- gtex_counts_raw$sample
gtex_mat  <- as.matrix(gtex_counts_raw[, -1])
rownames(gtex_mat) <- gene_ids

# Xena log2(count+1) formatında gelir → ham count'a çevir
gtex_mat_counts <- round(2^gtex_mat - 1)

message("Negatif değer var mı: ", any(gtex_mat_counts < 0))
message("İlk 5 gene, ilk 3 örnek:")
print(gtex_mat_counts[1:5, 1:3])



# TCGA-OV count matrisini al
load("data/processed/OV_se.RData")  # se objesi
tcga_counts <- assay(se, "counts")

message("TCGA-OV boyutu: ", nrow(tcga_counts), " x ", ncol(tcga_counts))

# Gene ID versiyonlarını temizle
rownames(tcga_counts)     <- sub("\\..*", "", rownames(tcga_counts))
rownames(gtex_mat_counts) <- sub("\\..*", "", rownames(gtex_mat_counts))

# Ortak genler
common_genes <- intersect(rownames(tcga_counts), rownames(gtex_mat_counts))
message("Ortak gen sayısı: ", length(common_genes))

tcga_counts     <- tcga_counts[common_genes, ]
gtex_mat_counts <- gtex_mat_counts[common_genes, ]

# Metadata
tcga_meta <- data.frame(
  sample_id   = colnames(tcga_counts),
  sample_type = "Tumor",
  source      = "TCGA",
  row.names   = colnames(tcga_counts),
  stringsAsFactors = FALSE
)

gtex_meta <- data.frame(
  sample_id   = colnames(gtex_mat_counts),
  sample_type = "Normal",
  source      = "GTEx",
  row.names   = colnames(gtex_mat_counts),
  stringsAsFactors = FALSE
)

# Birleştir
combined_counts <- cbind(tcga_counts, gtex_mat_counts)
combined_meta   <- rbind(tcga_meta, gtex_meta)
combined_meta$sample_type <- factor(combined_meta$sample_type, levels = c("Normal", "Tumor"))

message("Birleşik matris: ", nrow(combined_counts), " gen x ", ncol(combined_counts), " örnek")
message("Normal: ", sum(combined_meta$sample_type == "Normal"),
        " | Tumor: ", sum(combined_meta$sample_type == "Tumor"))

# SE oluştur ve kaydet
se_combined <- SummarizedExperiment(
  assays  = list(counts = combined_counts),
  colData = DataFrame(combined_meta)
)

save(se_combined, file = "data/processed/OV_se_with_gtex.RData")
message("Kaydedildi.")




library(DESeq2)

# se_combat'tan direkt DESeq2
dds <- DESeqDataSet(se_combat, design = ~ source + sample_type)
dds$sample_type <- relevel(dds$sample_type, ref = "Normal")
dds$source      <- factor(dds$source)

# Düşük eksprese gen filtresi
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message("Filtre sonrası gen sayısı: ", nrow(dds))

# DESeq2
message("DESeq2 çalışıyor...")
dds <- DESeq(dds)

# Sonuçlar
res <- results(dds,
               contrast      = c("sample_type", "Tumor", "Normal"),
               alpha         = padj_threshold)

summary(res)

# Histone genleri
histone_res <- as.data.frame(res[ensembl_ids, ])
histone_res$gene_name <- gene_map[rownames(histone_res)]
print(histone_res[, c("gene_name", "log2FoldChange", "pvalue", "padj")])
