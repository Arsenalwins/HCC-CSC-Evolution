# HCC-CSC-Evolution

Code repository for the manuscript:

**Clonal evolution and epigenetic memory of stem/progenitor cells in hepatocellular carcinoma: a single-cell transcriptomic study**

Yikai Hu, Kunjiang Tan, Na Feng, Jing Li, Lu Chen, Yuxuan Lin, Hongyang Wang*, Yufei He*

## Overview

This repository contains the analysis code for integrating single-cell RNA-seq data from 109 samples to construct a single-cell-resolution somatic mutation landscape in hepatocellular carcinoma (HCC), combined with targeted methylation sequencing (EPIC-seq) to trace the evolutionary trajectory from bipotent progenitors (BPs) to cancer stem cells (CSCs).

## Repository structure（Coming Soon......）
```
├── 01_data_processing/       # Cell Ranger, quality control, batch integration
├── 02_cell_annotation/        # Malignant cell identification, CSC/BP classification
├── 03_RNA_velocity/           # scVelo, veloVI, Monocle3, scTour pseudotime
├── 04_mutation_analysis/      # Monopogen SNV calling, mutational signatures, dN/dS
├── 05_CNV_phylogenetics/      # CopyKAT, MEDALT, phylogenetic tree construction
├── 06_cell_of_origin/         # TCR-based cell-of-origin inference
├── 07_epigenetic_analysis/    # EPIC-seq methylation, erosion rate, core imprint score
├── 08_WGCNA/                  # Weighted gene co-expression network analysis
├── 09_survival_analysis/      # Kaplan–Meier, Cox regression, timeROC
├── 10_tissue_microarray/      # Multiplex immunofluorescence quantification
└── utils/                     # Helper functions and shared utilities
```

## Key dependencies

- Python: scanpy, scvi-tools, scvelo, velovi, pyscenic, cytotrace2
- R: Seurat, CopyKAT, Monocle3, MutationalPatterns, methylKit, ComplexHeatmap, clusterProfiler, fgsea, timeROC

## Data availability

- Public scRNA-seq datasets: GSE149614, GSE156625, GSE282701, GSE242889 (GEO); SRP318499 (SRA)
- EpCAM-enriched scRNA-seq data: GSA [accession pending]

## Citation

[Manuscript under review]

## License

This project is licensed under the MIT License.
