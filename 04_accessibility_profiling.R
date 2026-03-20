suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(here)
})

in_peaks  <- here("github", "cyp2e1", "02_peak_linking", "tables", "cyp2e1_hub_peaks_metadata.csv")
in_seurat <- here("02_processed", "primary", "GSE281574_Liver_Multiome_Seurat_GEO.rds")
out_dir   <- here("github", "cyp2e1", "04_accessibility", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_seurat)
obj <- obj[, obj@meta.data$Condition == "Normal"]
obj <- obj[, !obj@meta.data$CellTypeSR %in% c("T cell", "Plasma Cell")]
obj@meta.data$CellTypeSR <- dplyr::case_when(
  obj@meta.data$CellTypeSR == "LSEC" ~ "EC",
  obj@meta.data$CellTypeSR == "Kupffer Cell" ~ "Kupffer Cells",
  obj@meta.data$CellTypeSR == "HSC" ~ "Stromal Cells",
  TRUE ~ obj@meta.data$CellTypeSR
)

peak_map <- read_csv(in_peaks, show_col_types = FALSE) %>% select(peak, Peak_ID)
atac_data <- GetAssayData(obj, assay = "ATAC", layer = "data")
cell_types <- unique(obj@meta.data$CellTypeSR)

stats_df <- bind_rows(lapply(cell_types, function(ct) {
  data.frame(
    peak = rownames(atac_data),
    Accessibility = rowMeans(atac_data[, obj@meta.data$CellTypeSR == ct, drop = FALSE]),
    CellType = ct
  )
})) %>%
  filter(peak %in% peak_map$peak) %>%
  left_join(peak_map, by = "peak")

write_csv(stats_df, file.path(out_dir, "peak_accessibility_by_celltype.csv"))