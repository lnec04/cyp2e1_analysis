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