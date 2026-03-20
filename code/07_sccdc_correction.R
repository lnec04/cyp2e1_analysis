suppressPackageStartupMessages({
  library(scCDC)
  library(Seurat)
  library(here)
})

# paths
seurat_rds_path  <- here("data", "GSE281574_Liver_Multiome_Seurat_GEO.rds")
out_intermediate <- here("export")
out_final        <- here("export", "GSE281574_Liver_scCDC_Corrected.rds")

dir.create(out_intermediate, recursive = TRUE, showWarnings = FALSE)

# load
seurat_obj           <- readRDS(seurat_rds_path)
DefaultAssay(seurat_obj) <- "RNA"

pools     <- unique(seurat_obj$sample)
cond_meta <- setNames(seurat_obj$Condition, colnames(seurat_obj))

# pool-wise correction
corrected_files <- c()

for (p_id in pools) {
  pool_out <- file.path(out_intermediate, p_id)
  dir.create(pool_out, recursive = TRUE, showWarnings = FALSE)

  seurat_sub          <- subset(seurat_obj, subset = sample == p_id)
  Idents(seurat_sub)  <- seurat_sub$CellTypeSR
  seurat_sub[["RNA"]] <- as(seurat_sub[["RNA"]], "Assay5")
  seurat_sub[["RNA"]] <- JoinLayers(seurat_sub[["RNA"]])
  seurat_sub          <- NormalizeData(seurat_sub, verbose = FALSE)

  GCGs <- ContaminationDetection(seurat_sub)
  saveRDS(GCGs, file.path(pool_out, paste0("GCGs_", p_id, ".rds")))

  ContaminationQuantification(seurat_sub, rownames(GCGs))

  seurat_corrected               <- ContaminationCorrection(seurat_sub, rownames(GCGs))
  DefaultAssay(seurat_corrected) <- "Corrected"

  atac_assays <- names(seurat_corrected@assays)[
    sapply(names(seurat_corrected@assays), function(a)
      inherits(seurat_corrected[[a]], "ChromatinAssay"))
  ]
  for (a in atac_assays) seurat_corrected[[a]] <- NULL

  out_path        <- file.path(pool_out, paste0("seurat_corrected_", p_id, ".rds"))
  saveRDS(seurat_corrected, out_path)
  corrected_files <- c(corrected_files, out_path)

  rm(seurat_sub, GCGs, seurat_corrected)
  gc()
}

rm(seurat_obj)
gc()

# merge
seurat_merged                <- readRDS(corrected_files[1])
merge_list                   <- lapply(corrected_files[-1], readRDS)
seurat_merged                <- merge(x = seurat_merged, y = merge_list)
DefaultAssay(seurat_merged)  <- "Corrected"
seurat_merged[["Corrected"]] <- JoinLayers(seurat_merged[["Corrected"]])
seurat_merged$Condition      <- cond_meta[colnames(seurat_merged)]

saveRDS(seurat_merged, out_final)
