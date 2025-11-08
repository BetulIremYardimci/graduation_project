# graduation_project
Pan-cancer analysis of Histone Variants

This script automates:
- Download of TCGA-BRCA RNA-Seq HTSeq-Counts (public) via the GDC REST API
- Building a counts matrix
- Differential expression (Tumor vs Normal) for H2AFX, H2AFY, H2AFY2, H2AFZ
- Kaplan–Meier survival with log-rank tests (median split)
- CSV + PNG outputs in `out/`