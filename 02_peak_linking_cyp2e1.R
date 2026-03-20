suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(here)
})

# Paths
seurat_rds_path <- here("02_processed", "primary", "GSE281574_Liver_Multiome_Seurat_GEO.rds")
out_dir_tables  <- here("github", "cyp2e1", "02_peak_linking", "tables")
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)

# Load and Harmonize
obj <- readRDS(seurat_rds_path)
obj <- obj[, obj@meta.data$Condition == "Normal"]
obj <- obj[, !obj@meta.data$CellTypeSR %in% c("T cell", "Plasma Cell")]
obj@meta.data$CellTypeSR <- dplyr::case_when(
  obj@meta.data$CellTypeSR == "LSEC" ~ "EC",
  obj@meta.data$CellTypeSR == "Kupffer Cell" ~ "Kupffer Cells",
  obj@meta.data$CellTypeSR == "HSC" ~ "Stromal Cells",
  TRUE ~ obj@meta.data$CellTypeSR
)

DefaultAssay(obj) <- "ATAC"
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# LinkPeaks
obj_hep <- obj[, obj@meta.data$CellTypeSR == "Hepatocyte"]
obj_hep <- LinkPeaks(obj_hep, peak.assay = "ATAC", expression.assay = "RNA", genes.use = "CYP2E1", distance = 1e6)

# Extract P1-P8 Hub
hub_peaks <- as.data.frame(Links(obj_hep)) %>%
  filter(score > 0 & pvalue < 0.05) %>%
  arrange(desc(score)) %>%
  mutate(Peak_ID = paste0("P", row_number()), Peak = gsub(":", "-", peak))

coords <- do.call(rbind, lapply(strsplit(hub_peaks$Peak, "-"), function(x) data.frame(chr = x[1], start = as.numeric(x[2]), end = as.numeric(x[3]))))
hub_peaks <- cbind(hub_peaks, coords)

write_csv(hub_peaks, file.path(out_dir_tables, "cyp2e1_hub_peaks_metadata.csv"))