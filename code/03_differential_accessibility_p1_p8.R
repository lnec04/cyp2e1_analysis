suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(here)
})

seurat_path <- here("data", "GSE281574_Liver_Multiome_Seurat_GEO.rds")
peaks_path  <- here("export", "02_cyp2e1_hub_peaks_metadata.csv")
out_path    <- here("export", "03_da_peaks_hepatocyte_enriched_stats.csv")

peak_map      <- read_csv(peaks_path, show_col_types = FALSE)
peaks_to_test <- peak_map$peak

obj <- readRDS(seurat_path)
obj <- obj[, obj@meta.data$Condition == "Normal"]
obj <- obj[, !grepl("plasma|^T[ _]|^T$|T.cell|Tcell", obj@meta.data$CellTypeSR, ignore.case = TRUE)]
obj@meta.data$CellTypeSR <- case_when(
  obj@meta.data$CellTypeSR == "LSEC"         ~ "EC",
  obj@meta.data$CellTypeSR == "Kupffer Cell" ~ "Kupffer Cells",
  obj@meta.data$CellTypeSR == "HSC"          ~ "Stromal Cells",
  TRUE ~ obj@meta.data$CellTypeSR
)

DefaultAssay(obj) <- "ATAC"
Idents(obj) <- "CellTypeSR"

latent_var <- if ("nCount_ATAC" %in% colnames(obj@meta.data)) "nCount_ATAC" else "nCount_peaks"

da_results <- FindMarkers(
  object          = obj,
  ident.1         = "Hepatocyte",
  features        = peaks_to_test,
  test.use        = "LR",
  latent.vars     = latent_var,
  logfc.threshold = 0,
  min.pct         = 0
)

da_results$peak <- rownames(da_results)

final_da_stats <- da_results %>%
  left_join(peak_map %>% select(peak, Peak_ID), by = "peak") %>%
  mutate(
    is_hepatocyte_enriched = p_val_adj < 0.05 & avg_log2FC > 0
  ) %>%
  select(Peak_ID, peak, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, is_hepatocyte_enriched) %>%
  arrange(desc(avg_log2FC))

write_csv(final_da_stats, out_path)
