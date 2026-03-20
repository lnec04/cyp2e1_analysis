suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(chromVAR)
  library(SummarizedExperiment)
  library(Matrix)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
  library(decoupleR)
  library(dorothea)
  library(tidyr)
  library(here)
})

# PARALLELISIERUNG
n_cores <- 8
 register(BiocParallel::MulticoreParam(
    workers  = n_cores,
    RNGseed  = 42,
    progressbar = TRUE
  ))


# PFADE
seurat_rds_path <- here("02_processed", "primary", "GSE281574_Liver_Multiome_Seurat_GEO.rds")
out_dir_tables  <- here("github", "cyp2e1", "01_tf_activity", "tables")
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# DATEN LADEN & HARMONISIEREN
obj <- readRDS(seurat_rds_path)

obj <- obj[, obj@meta.data$Condition == "Normal"]
obj <- obj[, !obj@meta.data$CellTypeSR %in% c("T cell", "Plasma Cell")]

obj@meta.data$CellTypeSR <- dplyr::case_when(
  obj@meta.data$CellTypeSR == "LSEC"         ~ "EC",
  obj@meta.data$CellTypeSR == "Kupffer Cell" ~ "Kupffer Cells",
  obj@meta.data$CellTypeSR == "HSC"          ~ "Stromal Cells",
  TRUE ~ obj@meta.data$CellTypeSR
)

# CHROMVAR
DefaultAssay(obj) <- "ATAC"
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

motif_map <- data.frame(
  motif_id  = names(pfm),
  tf_primary = vapply(pfm, function(x) {
    strsplit(name(x), split = "::")[[1]][1]
  }, character(1))
)

motif_matrix <- motifMatches(
  matchMotifs(
    pwms    = pfm,
    subject = granges(obj),
    genome  = BSgenome.Hsapiens.UCSC.hg38
  )
)
rownames(motif_matrix) <- rownames(obj)
Motifs(obj) <- CreateMotifObject(data = motif_matrix, pwm = pfm)

# chromVAR direkt — umgeht RunChromVAR Seurat-Wrapper (slot/layer Konflikt)
counts_mat  <- GetAssayData(obj, assay = "ATAC", layer = "counts")
peaks_keep  <- rowSums(counts_mat) > 0
counts_filt <- counts_mat[peaks_keep, ]
peak_ranges <- granges(obj[["ATAC"]])[peaks_keep]
motif_matrix_filt <- motif_matrix[peaks_keep, ]

se <- SummarizedExperiment::SummarizedExperiment(
  assays    = list(counts = counts_filt),
  rowRanges = peak_ranges
)
SummarizedExperiment::colData(se)$depth <- Matrix::colSums(counts_filt)
se  <- chromVAR::addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg  <- chromVAR::getBackgroundPeaks(se, niterations = 200, w = 0.1, bs = 50)
dev <- chromVAR::computeDeviations(
  object           = se,
  annotations      = motif_matrix_filt,
  background_peaks = bg
)

chromvar_scores <- chromVAR::deviationScores(dev)
colnames(chromvar_scores) <- colnames(obj)
obj[["chromvar"]] <- SeuratObject::CreateAssayObject(data = chromvar_scores)

# CHROMVAR PSEUDOBULK — alle Zelltypen (für cross-lineage accessibility)
chromvar_data  <- GetAssayData(obj, assay = "chromvar", layer = "data")
cell_types_vec <- obj$CellTypeSR

chrom_pb_all <- as.data.frame(
  sapply(sort(unique(cell_types_vec)), function(ct)
    rowMeans(chromvar_data[, cell_types_vec == ct, drop = FALSE])
  )
)
chrom_pb_all$motif_id <- rownames(chrom_pb_all)
rownames(chrom_pb_all) <- NULL

# Exportiere vollständige Matrix für cross-lineage heatmap
chromvar_tf_all <- chrom_pb_all %>%
  left_join(motif_map, by = "motif_id") %>%
  filter(!is.na(tf_primary)) %>%
  group_by(TF = tf_primary) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

write.csv(
  chromvar_tf_all,
  file.path(out_dir_tables, "chromvar_pseudobulk_all_celltypes.csv"),
  row.names = FALSE
)

# Z-SCORE NUR FÜR HEPATOZYTEN (CYP2E1 hepatozyten-spezifisch)
z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

chromvar_z <- chromvar_tf_all %>%
  mutate(across(-TF, z_score)) %>%
  select(TF, chromvar_z = Hepatocyte)

# DOROTHEA
rna_pb <- AggregateExpression(
  object      = obj,
  assays      = "RNA",
  features    = rownames(obj[["RNA"]]),
  group.by    = "CellTypeSR",
  return.seurat = FALSE
)$RNA

rna_counts <- as.data.frame(rna_pb)
lib_sizes  <- colSums(rna_counts)
rna_pb     <- log1p(sweep(rna_counts, 2, lib_sizes / 1e6, "/"))

regulons <- dorothea::dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  dplyr::select(source = tf, target = target, mor = mor)

tf_activities <- run_ulm(
  mat     = as.matrix(rna_pb),
  network = regulons,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5,
  center  = TRUE
)

tf_wide <- tf_activities %>%
  filter(statistic == "ulm") %>%
  select(TF = source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score)

# Z-Score nur Hepatozyten
hep_col <- grep("^Hepatocyte", colnames(tf_wide), value = TRUE)[1]

dorothea_z <- tf_wide %>%
  mutate(across(-TF, z_score)) %>%
  select(TF, dorothea_z = !!sym(hep_col))

# INTEGRATION & LINEAGE DRIVERS
comparison_all <- inner_join(chromvar_z, dorothea_z, by = "TF") %>%
  mutate(combined_score = (chromvar_z + dorothea_z) / 2)

write.csv(
  comparison_all,
  file.path(out_dir_tables, "tf_activity_scores_comparison.csv"),
  row.names = FALSE
)

# LINEAGE DRIVERS
hep_col_rna <- grep("^Hepatocyte", colnames(rna_pb), value = TRUE)[1]
avg_expr_hep <- rna_pb[intersect(comparison_all$TF, rownames(rna_pb)), hep_col_rna]

lineage_drivers <- comparison_all %>%
  filter(chromvar_z > 1.0 & dorothea_z > 1.0) %>%
  mutate(avg_expr = avg_expr_hep[TF]) %>%
  arrange(desc(combined_score))

write.csv(
  lineage_drivers,
  file.path(out_dir_tables, "top_hepatocyte_lineage_drivers.csv"),
  row.names = FALSE
)