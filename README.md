# cyp2e1_figure3_analysis

This repository contains the bioinformatics pipeline accompanying the manuscript: "Integrated single-cell and epigenomic analysis reveals cellular specificity of CYP2E1 in human alcohol-related liver disease". It details the computational reconstruction of the epigenetic and transcriptional regulatory architecture governing CYP2E1 expression.
Pipeline Architecture

## Workflow

    01_tf_activity_inference.R: Infers putative lineage drivers through the integration of chromatin accessibility (chromVAR) and transcription factor target gene expression (decoupleR/DoRothEA via a univariate linear model).

    02_peak_linking_cyp2e1.R: Identifies the CYP2E1 regulatory enhancer hub (P1–P8 peaks) in pericentral hepatocytes utilizing the LinkPeaks framework to correlate distal chromatin accessibility with local gene expression.

    03_motif_analysis.R: Executes iterative motif mapping of the JASPAR core vertebrate database to the delineated CYP2E1 enhancer elements.

    04_accessibility_profiling.R: Quantifies cross-lineage chromatin accessibility of the P1–P8 elements across parenchymal and non-parenchymal cellular compartments.

    05_pipeline_synthesis.R: Integrates transcription factor activity scores with motif occupancy data to determine the combinatorial binding architecture at the CYP2E1 locus.

## Computational Environment and Reproducibility

    Operating System: The pipeline was developed and validated under Linux.

    Environment Management: To ensure strict computational reproducibility, an renv.lock file is provided. Execute renv::restore() within the R console to replicate the exact package ecosystem and dependencies utilized in this study.

    Data Provisioning: Path management is handled via the here package. The pre-processed Seurat Multiome object (GSE281574_Liver_Multiome_Seurat_GEO.rds) must be placed in the 02_processed/primary/ directory prior to execution. All tabular outputs are systematically generated within their respective tables/ subdirectories.
