# Clonal Evolution and Epigenetic Memory of Stem/Progenitor Cells in Hepatocellular Carcinoma

This repository contains the analysis code for the manuscript:

> **Clonal evolution and epigenetic memory of stem/progenitor cells in hepatocellular carcinoma**
>
> YiKai Hu _et al._

## Overview

We integrated single-cell transcriptomic data from 109 samples (44 patients; 410,608 cells) to construct a somatic mutation landscape at single-cell resolution (384,867 SNVs across 31,908 cells from 20 patients). Combined with epigenetic profiling, we traced the evolutionary trajectory from bipotent progenitors (BPs) to cancer stem cells (CSCs), and validated the model using EpCAM-enriched scRNA-seq from two additional patients. Dynamic methylation experiments revealed that ~78% of CSC-specific methylation alterations exhibit stable epigenetic memory.

## Repository Structure

```
HCC_CSC_Evolution/
├── README.md
└── code/
    ├── Fig1_Atlas/              # Single-cell atlas and CSC/BP identification
    │   └── Fig1.py              # Fig 1b–f
    │
    ├── Fig2_Mutation/           # Somatic mutations inferred from scRNA-seq
    │   └── Fig2.R               # Fig 2a–f
    │
    ├── Fig3_CNV_COO/            # CNV and cell-of-origin analyses
    │   ├── Fig3_cnv.py          # Fig 3a–c
    │   └── Fig3.R               # Fig 3d–g
    │
    ├── Fig4_MACS_Validation/    # EpCAM-enriched sample validation
    │   ├── Fig4.py              # Fig 4a, 4c (Python)
    │   └── Fig4c.R              # Fig 4c (Monocle3)
    │
    ├── Fig5_WGCNA_DNMT/         # DNMT and CSC stemness maintenance
    │   ├── Fig5.R               # Fig 5a–j
    │   └── Fig5_WGCNA.py        # Fig 5c (WGCNA plotting)
    │
    ├── Fig6_Methylation/        # Dynamic methylation and epigenetic memory
    │   └── Fig6.R               # Fig 6b–k
    │
    ├── Fig7_TMA/                # Tissue microarray immunofluorescence
    │   └── Fig7.R               # Fig 7h–j
    │
    └── Supplementary/           # Extended data figures
        ├── FigS1.py             # Fig S1a–d
        ├── FigS2.py             # Fig S2a–c, e–h
        ├── FigS2d.R             # Fig S2d
        ├── FigS3_S4.R           # Fig S3a, S4a–c
        └── FigS6_S7.R           # Fig S6, S7a–d
```

## Figure–Script Mapping

| Figure | Script | Description |
|--------|--------|-------------|
| Fig 1b–f | `Fig1_Atlas/Fig1.py` | scANVI integration, CytoTRACE2 KDE, AUCell heatmap, stacked violin, RNA velocity |
| Fig 2a | `Fig2_Mutation/Fig2.R` | Single-cell OncoPrint + dN/dS lollipop |
| Fig 2b | `Fig2_Mutation/Fig2.R` | Variant classification distribution |
| Fig 2c | `Fig2_Mutation/Fig2.R` | De novo mutational signatures (NMF) |
| Fig 2d–f | `Fig2_Mutation/Fig2.R` | Shared-site O/E mutation burden (LMM + Wilcoxon) |
| Fig 3a–c | `Fig3_CNV_COO/Fig3_cnv.py` | CNV score paired comparisons |
| Fig 3d,e | `Fig3_CNV_COO/Fig3.R` | MEDALT evolutionary distance |
| Fig 3f,g | `Fig3_CNV_COO/Fig3.R` | Cell-of-origin inference (TCR-based) |
| Fig 4a | `Fig4_MACS_Validation/Fig4.py` | EpCAM-enriched RNA velocity |
| Fig 4c | `Fig4_MACS_Validation/Fig4.py` + `Fig4c.R` | SNP-transcriptome pseudotime |
| Fig 5a–j | `Fig5_WGCNA_DNMT/Fig5.R` + `Fig5_WGCNA.py` | WGCNA, DNMT regulatory chain, KM survival |
| Fig 6b–k | `Fig6_Methylation/Fig6.R` | Methylation entropy, ChromHMM, memory classification, AUC |
| Fig 7h–j | `Fig7_TMA/Fig7.R` | MIF quantification and prognosis |
| Fig S1a–d | `Supplementary/FigS1.py` | Atlas marker dotplot, NMF clustering, MP analysis |
| Fig S2a–h | `Supplementary/FigS2.py` + `FigS2d.R` | BP gating, velocity, scTour, Monocle3, CytoTRACE2 |
| Fig S3, S4 | `Supplementary/FigS3_S4.R` | SNV lollipop, CNV phylogenetic trees |
| Fig S6, S7 | `Supplementary/FigS6_S7.R` | GO enrichment, entropy pipeline, compartment, epigenetics |
| Fig S8 | `Fig7_TMA/Fig7.R` | Adjacent tissue IF and prognosis |

> **Note:** Fig 1a, 6a, 8 are schematic diagrams created manually. Fig 4b (CNV tree for MACS patients) uses the same pipeline as Fig S4 with patient-specific inputs. Fig 7a–g and S8a–g are microscopy images from multiplex immunofluorescence.

## Software Dependencies

### Python
- scanpy (≥1.10.3), scvi-tools (≥1.1.6), scvelo (≥0.3.3)
- sctour (≥1.0.0), velovi (≥0.3.1), cellrank
- numpy, pandas, scipy, scikit-learn, matplotlib, seaborn

### R
- Seurat (≥5.1.0), monocle3 (≥1.4.26)
- ComplexHeatmap, MutationalPatterns, maftools
- survival, survminer, timeROC
- methylKit, clusterProfiler, fgsea
- lmerTest, data.table, ggplot2, ggpubr
- phangorn, ggtree, ape
- compartmap, GenomicRanges, rtracklayer

## Data Availability

- **Public scRNA-seq**: GSE149614, GSE156625, GSE282701, GSE242889 (GEO); SRP318499 (SRA)
- **EpCAM-enriched data**: deposited in GSA database
- **TCGA-LIHC / ICGC LIRI-JP / GSE14520**: downloaded via [HCCDB](http://lifeome.net/database/hccdb/home.html)
- **CLCA WGS**: Genome Sequence Archive (PRJCA002666)

## License

This project is licensed under the MIT License.

## Contact

- Hu Yikai: huyikai2001@gmail.com
