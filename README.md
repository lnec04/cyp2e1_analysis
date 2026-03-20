# cyp2e1_figure3_analysis

Bioinformatics pipeline accompanying the manuscript: *"Integrated single-cell and epigenomic analysis reveals cellular specificity of CYP2E1 in human alcohol-related liver disease"*. Reconstructs the transcriptional and epigenetic regulatory architecture governing CYP2E1 expression.

## Data availability

The pre-processed Seurat Multiome object (`GSE281574_Liver_Multiome_Seurat_GEO.rds`, ~2 GB) is **not included** in this repository. Download it from GEO accession [GSE281574](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281574), provided by Chembazhi et al. (2025, *Nature Communications*, https://doi.org/10.1038/s41467-025-63251-2).

After download, place the file at:

```
data/GSE281574_Liver_Multiome_Seurat_GEO.rds
```

## Pipeline

### Figure 3 analysis (`01 → 06`)

Scripts must be executed sequentially. All tabular outputs are written to `export/`.

| Script                                    | Function                                                                                                                                    | Output                                                                                                                                                 |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `01_tf_activity_inference.R`            | Integrates chromVAR (motif accessibility) and DoRothEA (TF target expression via ULM) to identify hepatocyte lineage drivers                | `export/01_chromvar_pseudobulk_all_celltypes.csv<br>``export/01_tf_activity_scores_comparison.csv<br>``export/01_top_hepatocyte_lineage_drivers.csv` |
| `02_peak_linking_cyp2e1.R`              | Identifies the CYP2E1 regulatory enhancer hub in untreated hepatocytes via Signac LinkPeaks (1 Mb window, score > 0, FDR-adjusted p < 0.05) | `export/02_cyp2e1_hub_peaks_metadata.csv`                                                                                                            |
| `03_differential_accessibility_p1_p8.R` | Validates hepatocyte-enriched accessibility of P1–P8 via logistic regression adjusting for ATAC library size bias (adjusted p < 0.05)      | `export/03_da_peaks_hepatocyte_enriched_stats.csv`                                                                                                   |
| `04_motif_analysis.R`                   | Maps JASPAR2020 Homo sapiens motifs to P1–P8 elements; constructs binary occupancy matrix                                                  | `export/04_motif_occupancy_matrix_p1_p8.csv`                                                                                                         |
| `05_accessibility_profiling.R`          | Quantifies baseline ATAC accessibility (prevalence and magnitude) of P1–P8 across all cell types                                           | `export/05_peak_accessibility_by_celltype.csv`                                                                                                       |
| `06_pipeline_synthesis.R`               | Intersects TF lineage drivers with motif occupancy to determine combinatorial binding architecture at the CYP2E1 locus                      | `export/06_final_tf_peak_binding_synthesis.csv`                                                                                                      |

### Ambient RNA decontamination

Independent analysis. Can be run separately from the Figure 3 pipeline.

| Script                    | Function                                                                                                                                                                               | Output                                          |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| `07_sccdc_correction.R` | Pool-wise ambient RNA correction via scCDC (ContaminationDetection → Quantification → Correction); corrected pools merged into single Seurat object with `Corrected` assay as default | `export/GSE281574_Liver_scCDC_Corrected.rds`  |

## Reproducibility

- **OS:** Linux
- **Environment:** Restore exact package versions via `renv::restore()` using the provided `renv.lock`
- **Paths:** Managed via the `here` package; all paths are relative to the repository root
