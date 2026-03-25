suppressPackageStartupMessages({
  library(here)
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(motifmatchr)
  library(dplyr)
  library(readr)
  library(tidyr)
})

in_peaks_csv <- here("export", "02_cyp2e1_hub_peaks_metadata.csv")

peaks_df <- read_csv(in_peaks_csv, show_col_types = FALSE)
pfm_all  <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", species = 9606, all_versions = FALSE))

# Motif Matching
gr_fg <- StringToGRanges(peaks_df$Peak)
seq_fg <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_fg)
mat_fg <- motifMatches(matchMotifs(pfm_all, seq_fg, genome = BSgenome.Hsapiens.UCSC.hg38, out = "matches"))

hits_list <- lapply(seq_len(nrow(mat_fg)), function(i) {
  idx <- which(mat_fg[i, ])
  if (length(idx) > 0) {
    data.frame(Peak_ID = peaks_df$Peak_ID[i], TF = toupper(name(pfm_all)[idx]), Found = 1)
  }
})

occupancy_matrix <- bind_rows(hits_list) %>%
  distinct() %>%
  pivot_wider(names_from = TF, values_from = Found, values_fill = 0)

write_csv(occupancy_matrix, here("export", "05_motif_occupancy_matrix_p1_p8.csv"))

# Fisher's exact test: driver TF motif enrichment in P1-P8 vs background
obj        <- readRDS(here("data", "GSE281574_Liver_Multiome_Seurat_GEO.rds"))
obj        <- obj[, obj@meta.data$Condition == "Normal"]
obj        <- obj[, obj@meta.data$CellTypeSR == "Hepatocyte"]
DefaultAssay(obj) <- "ATAC"
all_peaks  <- granges(obj)

driver_tfs <- c("HNF4A", "HNF4G", "HNF1A", "SREBF2", "PPARA")
tf_names   <- sapply(pfm_all, function(x) x@name)
pfm_driver <- pfm_all[tf_names %in% driver_tfs]

bg_matrix   <- motifMatches(matchMotifs(pfm_driver, all_peaks,  genome = BSgenome.Hsapiens.UCSC.hg38, out = "matches"))
p1p8_matrix <- motifMatches(matchMotifs(pfm_driver, gr_fg,      genome = BSgenome.Hsapiens.UCSC.hg38, out = "matches"))

fisher_results <- do.call(rbind, Filter(Negate(is.null), lapply(seq_along(colnames(p1p8_matrix)), function(i) {
  tf      <- colnames(p1p8_matrix)[i]
  col_idx <- which(colnames(bg_matrix) == tf)
  if (length(col_idx) == 0) return(NULL)

  p1p8_with    <- sum(p1p8_matrix[, i])
  p1p8_without <- nrow(p1p8_matrix) - p1p8_with
  bg_with      <- sum(bg_matrix[, col_idx])
  bg_without   <- nrow(bg_matrix) - bg_with

  ft <- fisher.test(
    matrix(c(p1p8_with, p1p8_without, bg_with, bg_without), nrow = 2),
    alternative = "greater"
  )

  data.frame(
    TF              = tf,
    P1P8_with_motif = p1p8_with,
    P1P8_total      = nrow(p1p8_matrix),
    BG_with_motif   = bg_with,
    BG_total        = nrow(bg_matrix),
    Pct_P1P8        = round(p1p8_with / nrow(p1p8_matrix) * 100, 1),
    Pct_BG          = round(bg_with   / nrow(bg_matrix)   * 100, 1),
    OddsRatio       = round(ft$estimate, 2),
    P_value         = ft$p.value
  )
})))

write_csv(fisher_results, here("export", "05_motif_enrichment_fisher.csv"))