#!/usr/bin/env Rscript
# Figure 5 — DNMT upregulation is associated with CSC stemness maintenance
# 5a: GO/GSEA enrichment
# 5b: DNMT/TET expression dotplot
# 5c: WGCNA module Hallmark enrichment
# 5d: KM survival (BP-M1 in adjacent)
# 5e: BP-M1 vs Hoshida poor-prognosis correlation
# 5f-j: DNMT → promoter methylation → gene silencing regulatory chain

#!/usr/bin/env Rscript
#
#
#   Fig1A  ME Cell-Type Correlation Heatmap
# Fig1B  Module Hallmark Enrichment ()
#   Fig2A  KM Survival: BP-M1 score in Adjacent
#   Fig2B  Scatter: BP score vs Hoshida score (ρ=0.71)
#   Fig3A  Volcano: CSC vs Other Tumor
#   Fig3B  Scatter: DNMT score vs CSC hub score (ρ=0.559)
# Fig3C  4-panel: DNMT→DiffBeta→DiffExpr
#   Fig3D  KM ×3: DNMT alone / CSC alone / DNMT×CSC 4-group
#   Fig3E  Forest plot: HR + 95% CI
# Fig4A  DotPlot: 14  across CSC/BP/Hepatocytes
# Fig4B  Boxplot + Heatmap
#   Fig5A  Scatter: EpiSCORE CSC fraction vs CSC-M1 (Tumor)
#   Fig5B  Scatter: EpiSCORE BP fraction vs BP-M1 (Adjacent)
#
# Python  (celltype_specific_wgcna / similarity_enrichment)
#

###############################################################################
# S0.                                #
###############################################################################

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c(
  # Bioconductor
  "TCGAbiolinks", "SummarizedExperiment",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "msigdbr",
  # Seurat ()
  "Seurat",
  "survival", "survminer",
  "dplyr", "tidyr", "tibble", "janitor",
  "ggplot2", "ggpubr", "ggrepel", "ggVennDiagram",
  "pheatmap", "RColorBrewer", "corrplot",
  "gridExtra", "grid"
)

for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE),
      error = function(e) install.packages(pkg)
    )
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

###############################################################################
# S1.  ( )                         #
###############################################################################

WGCNA_DIR    <- "/data1/WGCNA/CellType_Specific_Results/"
OUTPUT_DIR   <- "/data1/WGCNA/CellType_Specific_Results/"
FIG_DIR      <- "/home/download/csc_article/fig4/wgcna/"
TCGA_DIR     <- "~/JWD_LIHC/"
EPISCORE_DIR <- file.path(OUTPUT_DIR, "Plan2_EpiSCORE/")
RDATA_PATH   <- "/data2/adata/atlas.RData"          # Seurat

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(TCGA_DIR)

# --- WGCNA hub gene  ---
N_HUB        <- 30
CUTOFF_YEARS <- 5
cutoff_days  <- CUTOFF_YEARS * 365.25

FOCUS_MODULES <- list(
  darkturquoise = list(label = "CSC-M1",  group = "CSC",      color = "#D62728"),
  lightcyan     = list(label = "CSC-M2",  group = "CSC",      color = "#FF7F0E"),
  saddlebrown   = list(label = "BP-M1",   group = "Bipotent", color = "#1F77B4"),
  darkorange    = list(label = "BP-M2",   group = "Bipotent", color = "#2CA02C")
)

DNMT_GENES <- c("DNMT1", "DNMT3A", "DNMT3B")

EPI_GENES <- c("TET1","TET2","TET3","DNMT1","DNMT3A","DNMT3B",
               "UHRF1","TDG","GADD45A","GADD45B","HDAC1","IDH1","IDH2","EZH2")

# --- Hoshida 2008 NEJM  ---
HOSHIDA_POOR <- c(
  "FSHB","SH3GL2","RBM34","NCAPH","EGF","TRIO","COL6A3","ABLIM1","ITGA9",
  "NTS","SERPINB2","MMP7","PRKG2","EDG4","NOS2A","EPHA4","SP100","FMO1",
  "GPX2","ATP4B","IER3","WIPF1","IGFBP6","CTNND2","FBN1","GBA","ICK",
  "CHERP","HDAC9","NOL7","IQGAP1","ADD3","ANXA3","HMG20B","SLC12A2",
  "COL4A1","CPA3","KRT7","SERPINB8","NFKB2","AEBP1","TGFB1I1","EMP2",
  "BCL2","PSMB9","ACTR2","DDR1","SLC7A1","PODXL","COL16A1","IFI30",
  "EPM2AIP1","ANXA1","CCL21","CHSY1","AP1B1","TEAD4","ELOVL2","TCF4",
  "TSC22D2","CCDC6","CD48","TNK2","DAB2","LOXL2","RNASE1","LPP","CXCR4",
  "SLIT3","FILIP1L","CCL19","AKAP13"
)

HOSHIDA_GOOD <- c(
  "ALDH9A1","TTR","RLF","IMPA1","PFKFB1","ACSM3","ADRA2B","PTPN2","PSMB3",
  "PPP1R1A","TMEM97","PKLR","RPS6KA5","CYB5A","SCG5","TXN2","PLG","SC5DL",
  "AR","IGF1","SUCLG1","HAAO","C9","TAF1C","CPOX","XPA","HABP2","GHR","PCK1",
  "AKR1D1","ADH5","AARS","C8B","MGC29506","ATP6AP2","DOCK4","PROS1","ZBTB17",
  "DAD1","TIMM8A","HMGCL","C4BPB","TRAF6","EIF2B1","LIPC","PIGK","WDR23",
  "RFC2","GRM5","SDHC","ERCC5","F9","ANKRD46","ART1","CTBS","SLC37A4",
  "ALAS1","VPS41","GCGR","CCT8","BRP44","GRK4","HSPE1","NARS2","DST","ATP2C1",
  "AKR1A1","EMD","CALCR","DLGAP4","RRM1","NENF","SNX10","PMM1","TDO2","GSTM1",
  "SREBF2","PTPN18","ASAHL","PLCG2","KCNJ3","PCYT2","GJB1","TM7SF2",
  "SELENBP1","AOX1","ZER1","ADH6","MSH6","SLC4A4","USP14","C5","RAD52",
  "FAM129A","BAIAP2","SSFA2","PON3","GCKR","CREB1","CUTL2","SFRS2","HMGCR",
  "GGCX","CYP2B6","ZNF185","ARF4","ACOT2","ATP5D","CPN1","PLCB3","INSM1",
  "POLRMT","HRASLS3"
)

###############################################################################
# S2.  WGCNA  +  Hub Gene Lists                      #
###############################################################################

message("\n", strrep("=", 70))
message(strrep("=", 70))

all_kme      <- read.csv(file.path(WGCNA_DIR, "Table_All_Gene_kME.csv"))
gene_modules <- read.csv(file.path(WGCNA_DIR, "gene_modules.csv"), row.names = 1)

# hub genes  ( pseudogene / lncRNA)
get_hub_genes <- function(mod_name, n = N_HUB) {
  all_kme %>%
    filter(module == mod_name) %>%
    filter(!grepl("^(LINC|AC[0-9]|AL[0-9]|AP[0-9]|RN7SL|CTD-|CTC-|RP[0-9]|KB-)", gene)) %>%
    arrange(desc(kME)) %>%
    head(n) %>%
    pull(gene)
}

hub_gene_lists <- lapply(setNames(names(FOCUS_MODULES), names(FOCUS_MODULES)), get_hub_genes)
csc_hub_genes  <- unique(c(hub_gene_lists[["darkturquoise"]], hub_gene_lists[["lightcyan"]]))
bp_hub_genes   <- unique(c(hub_gene_lists[["saddlebrown"]],  hub_gene_lists[["darkorange"]]))

for (m in names(FOCUS_MODULES)) {
  info <- FOCUS_MODULES[[m]]
  message(sprintf("  %s (%s): %d hub genes — top5: %s",
                  info$label, m, length(hub_gene_lists[[m]]),
                  paste(head(hub_gene_lists[[m]], 5), collapse = ", ")))
}

bp_m1_full  <- rownames(gene_modules)[gene_modules$moduleColors == "saddlebrown"]
bp_m2_full  <- rownames(gene_modules)[gene_modules$moduleColors == "darkorange"]
bp_module_full  <- unique(c(bp_m1_full, bp_m2_full))
csc_m1_full <- rownames(gene_modules)[gene_modules$moduleColors == "darkturquoise"]
csc_m2_full <- rownames(gene_modules)[gene_modules$moduleColors == "lightcyan"]
csc_module_full <- unique(c(csc_m1_full, csc_m2_full))

###############################################################################
# S3.  TCGA-LIHC  +  +                        #
###############################################################################

message("\n", strrep("=", 70))
message(strrep("=", 70))

# --- 3.1  ---
tcga_rdata <- file.path(OUTPUT_DIR, "TCGA_LIHC_full.RData")

if (file.exists(tcga_rdata)) {
  load(tcga_rdata)
} else {
  query <- GDCquery(
    project = "TCGA-LIHC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  data_obj <- GDCprepare(query)

  tcga_expr <- assay(data_obj, "fpkm_unstrand")
  rownames(tcga_expr) <- as.data.frame(rowData(data_obj))$gene_name

  clinical_raw <- as.data.frame(colData(data_obj))

  clinical_clean <- clinical_raw
  list_cols <- sapply(clinical_clean, is.list)
  clinical_clean[list_cols] <- lapply(clinical_clean[list_cols], function(x)
    sapply(x, function(y) paste(y, collapse = "; ")))

  save(tcga_expr, clinical_raw, clinical_clean, file = tcga_rdata)
}

list_cols <- sapply(clinical_raw, is.list)
if (any(list_cols)) {
  clinical_raw[list_cols] <- lapply(clinical_raw[list_cols], function(x)
    sapply(x, function(y) {
      if (is.null(y) || length(y) == 0) NA_character_
      else paste(y, collapse = "; ")
    }))
}

# clinical_raw  S4 DataFrame,  data.frame
if (!is.data.frame(clinical_raw) || is(clinical_raw, "DFrame")) {
  clinical_raw <- as.data.frame(clinical_raw, stringsAsFactors = FALSE)
}

if (exists("clinical_clean") && is.data.frame(clinical_clean) && ncol(clinical_clean) >= ncol(clinical_raw)) {
  clinical_raw <- clinical_clean
}

message(sprintf("  tcga_expr: %d genes × %d samples", nrow(tcga_expr), ncol(tcga_expr)))

# --- 3.2 OS  ---
safe_col <- function(df, col) {
  if (col %in% colnames(df)) df[[col]] else rep(NA, nrow(df))
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) df[[hit[1]]] else rep(NA, nrow(df))
}

df_os_raw <- data.frame(
  patient_id             = substr(as.character(find_col(clinical_raw, c("barcode","bcr_patient_barcode","patient"))), 1, 12),
  barcode                = as.character(find_col(clinical_raw, c("barcode","bcr_patient_barcode"))),
  sample_type            = as.character(find_col(clinical_raw, c("sample_type","shortLetterCode","definition"))),
  vital_status           = as.character(find_col(clinical_raw, c("vital_status"))),
  days_to_death          = as.numeric(find_col(clinical_raw, c("days_to_death"))),
  days_to_last_follow_up = as.numeric(find_col(clinical_raw, c("days_to_last_follow_up","days_to_last_followup"))),
  age                    = as.numeric(find_col(clinical_raw, c("age_at_index","age_at_diagnosis","age_at_initial_pathologic_diagnosis"))),
  gender                 = as.character(find_col(clinical_raw, c("gender","sex"))),
  stage                  = as.character(find_col(clinical_raw, c("ajcc_pathologic_stage","tumor_stage.diagnoses","ajcc_pathologic_t"))),
  grade                  = as.character(find_col(clinical_raw, c("paper_Histological.Grade","histological_grade","neoplasm_histologic_grade"))),
  stringsAsFactors       = FALSE
)

print(table(df_os_raw$sample_type, useNA = "ifany"))

df_os_raw <- df_os_raw %>%
  filter(grepl("Primary Tumor|Tumor|TP|01", sample_type, ignore.case = TRUE)) %>%
  mutate(
    os_status = ifelse(grepl("Dead|dead|1", vital_status), 1, 0),
    os_time   = ifelse(os_status == 1, days_to_death, days_to_last_follow_up),
    stage_simple = case_when(
      grepl("^Stage I$|^Stage IA|^Stage IB", stage)  ~ "I",
      grepl("^Stage II$|^Stage IIA|^Stage IIB", stage) ~ "II",
      grepl("^Stage III", stage) ~ "III",
      grepl("^Stage IV", stage)  ~ "IV",
      TRUE ~ NA_character_),
    grade_simple = case_when(
      grepl("G1", grade) ~ "G1", grepl("G2", grade) ~ "G2",
      grepl("G3", grade) ~ "G3", grepl("G4", grade) ~ "G4",
      TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(os_time) & os_time > 0) %>%
  distinct(patient_id, .keep_all = TRUE)

# --- 3.3 RFS (cBioPortal) ---
cbio_file <- "lihc_tcga_pan_can_atlas_2018_clinical_data.tsv"
cbio_clinical <- NULL
if (file.exists(cbio_file)) {
  cbio_clinical <- read.delim(cbio_file, check.names = FALSE) %>%
    clean_names() %>%
    select(patient_id, disease_free_months, disease_free_status) %>%
    mutate(
      rfs_status = case_when(
        grepl("1:Recurred", disease_free_status)   ~ 1,
        grepl("0:DiseaseFree", disease_free_status) ~ 0,
        TRUE ~ NA_real_),
      rfs_time = as.numeric(disease_free_months) * 30.44
    ) %>%
    filter(!is.na(rfs_time) & rfs_time > 0 & !is.na(rfs_status))
}

# --- 3.4 Methylation beta matrix ---
if (!exists("tcga_beta") || is.null(tcga_beta)) {
  tcga_beta <- NULL; probe_gene_map <- NULL
  beta_paths <- c(
    "/home/download/JWD_LIHC/tcga_lihc_methylation_final.rds",
    file.path(OUTPUT_DIR, "TCGA_beta.RData")
  )
  for (bp in beta_paths) {
    if (!file.exists(bp)) next
    message(sprintf("  Loading methylation: %s", bp))
    if (grepl("\\.rds$", bp, ignore.case = TRUE)) {
      obj <- readRDS(bp)
      if (inherits(obj, "SummarizedExperiment") || inherits(obj, "RangedSummarizedExperiment")) {
        message("  RDS is SummarizedExperiment, extracting assay...")
        tcga_beta <- SummarizedExperiment::assay(obj)
      } else if (is.list(obj) && "beta" %in% names(obj)) {
        tcga_beta <- obj$beta; probe_gene_map <- obj$probe_gene_map
      } else if (is.matrix(obj) || is.data.frame(obj)) {
        tcga_beta <- as.matrix(obj)
      } else {
        tcga_beta <- tryCatch(SummarizedExperiment::assay(obj), error = function(e) NULL)
      }
      rm(obj); gc(verbose = FALSE)
    } else {
      load(bp)
    }
    if (!is.null(tcga_beta)) break
  }

  # Build probe_gene_map from HM450 annotation (same logic as wc2methyl.R)
  if (!is.null(tcga_beta) && is.null(probe_gene_map)) {
    anno_file <- "/data1/WGCNA/CellType_Specific_Results/Hoshida_Integration/GSE157973_RAW_extracted/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
    if (file.exists(anno_file)) {
      message("  Building probe_gene_map from HM450 annotation...")
      # GPL13534 CSV has comment lines at top; find real header with IlmnID
      con <- gzfile(anno_file, "r")
      skip_n <- 0
      while (TRUE) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) break
        skip_n <- skip_n + 1
        if (grepl("^IlmnID", line)) { skip_n <- skip_n - 1; break }
      }
      close(con)
      message(sprintf("  Skipping %d header lines", skip_n))
      anno_raw <- data.table::fread(anno_file, sep = ",", header = TRUE,
        skip = skip_n, fill = TRUE,
        select = c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group"))
      promoter_anno <- anno_raw[grepl("TSS200|TSS1500|5'UTR|1stExon", anno_raw$UCSC_RefGene_Group), ]
      promoter_anno <- promoter_anno[!is.na(promoter_anno$UCSC_RefGene_Name) &
                                      promoter_anno$UCSC_RefGene_Name != "", ]
      probe_gene_map <- promoter_anno %>%
        as.data.frame() %>%
        tidyr::separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep = ";") %>%
        dplyr::mutate(probe_id = Name, gene = trimws(UCSC_RefGene_Name)) %>%
        dplyr::select(probe_id, gene) %>%
        dplyr::distinct()
      rm(anno_raw, promoter_anno); gc(verbose = FALSE)
      message(sprintf("  probe_gene_map: %d probes -> %d genes",
                      length(unique(probe_gene_map$probe_id)),
                      length(unique(probe_gene_map$gene))))
    } else {
      message("  WARNING: HM450 annotation not found at: ", anno_file)
    }
  }

  if (!is.null(tcga_beta)) {
    message(sprintf("  tcga_beta: %d probes x %d samples", nrow(tcga_beta), ncol(tcga_beta)))
  } else {
    message("  WARNING: methylation not loaded")
  }
}

# --- 3.5 GSE14520  (Hoshida ) ---
GSE14520_DIR <- file.path(OUTPUT_DIR, "Hoshida_Integration/")
gse_meta_file <- file.path(GSE14520_DIR, "GSE14520_metadata.txt")
gse_expr_file <- file.path(GSE14520_DIR, "GSE14520_expr.RData")

gse_expr   <- NULL
gse_clin   <- NULL
gse_probe_gene <- NULL

if (file.exists(gse_meta_file)) {
  gse_meta <- read.delim(gse_meta_file, stringsAsFactors = FALSE)
  message(sprintf("  GSE14520 metadata: %d samples", nrow(gse_meta)))

  if (file.exists(gse_expr_file)) {
    gse_env <- new.env()
    tryCatch({
      load(gse_expr_file, envir = gse_env)
      obj_names <- ls(gse_env)
      message(sprintf("  GSE14520 RData contains: %s", paste(obj_names, collapse = ", ")))
      # Try known variable names
      if ("expr_all" %in% obj_names) gse_expr <- gse_env$expr_all
      else if ("gse_expr" %in% obj_names) gse_expr <- gse_env$gse_expr
      else if (length(obj_names) >= 1) gse_expr <- get(obj_names[1], envir = gse_env)
      if ("probe_gene" %in% obj_names) gse_probe_gene <- gse_env$probe_gene
      else if ("gse_probe_gene" %in% obj_names) gse_probe_gene <- gse_env$gse_probe_gene
      else if (length(obj_names) >= 2) gse_probe_gene <- get(obj_names[2], envir = gse_env)
      message(sprintf("  GSE14520 expr: %d probes x %d samples", nrow(gse_expr), ncol(gse_expr)))
    }, error = function(e) {
      message(sprintf("  WARNING: GSE14520 RData load failed: %s", e$message))
      gse_expr <<- NULL
    })
    rm(gse_env)
  } else {
    if (requireNamespace("GEOquery", quietly = TRUE)) {
      tryCatch({
        library(GEOquery)
        gse <- getGEO("GSE14520", GSEMatrix = TRUE, getGPL = TRUE)
        gse_obj <- gse[[1]]
        gse_expr <- exprs(gse_obj)
        feat <- fData(gse_obj)
        gene_col <- intersect(c("Gene Symbol","Gene.Symbol","GENE_SYMBOL","gene_assignment"), colnames(feat))[1]
        if (!is.na(gene_col)) {
          gse_probe_gene <- data.frame(
            probe = rownames(feat),
            gene_clean = as.character(feat[[gene_col]]),
            stringsAsFactors = FALSE
          ) %>% filter(gene_clean != "" & !is.na(gene_clean))
        }
        save(expr_all = gse_expr, probe_gene = gse_probe_gene, file = gse_expr_file)
      }, error = function(e) message(sprintf("   GSE14520 download failed: %s", e$message)))
    }
  }

  if (!is.null(gse_expr)) {
    gse_clin <- gse_meta %>%
      filter(Tissue.Type == "Tumor" & Affy_GSM != "") %>%
      dplyr::select(Affy_GSM, Survival.status, Survival.months, Recurr.status, Recurr.months) %>%
      mutate(
        sample_id  = Affy_GSM,
        os_status  = as.numeric(Survival.status),
        os_time    = as.numeric(Survival.months),
        rfs_status = as.numeric(Recurr.status),
        rfs_time   = as.numeric(Recurr.months)
      ) %>%
      filter(!is.na(os_time) & os_time > 0)
    rownames(gse_clin) <- gse_clin$sample_id
    message(sprintf("  GSE14520 tumor samples: %d", nrow(gse_clin)))
  }
} else {
}

###############################################################################
# #
###############################################################################

# ---- z-score signature score ----
calc_sig_score <- function(genes, expr_mat, samples = NULL) {
  available <- intersect(genes, rownames(expr_mat))
  if (length(available) < 3) return(rep(NA, ncol(expr_mat)))
  sub <- expr_mat[available, , drop = FALSE]
  if (!is.null(samples)) sub <- sub[, samples, drop = FALSE]
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

plot_km <- function(df, time_col, status_col, group_col = "group",
                    title = "", palette = c("#4DBBD5", "#E64B35"),
                    legend_labs = c("Low", "High"), x_max = NULL) {
  df$time   <- df[[time_col]]
  df$status <- df[[status_col]]
  df$grp    <- df[[group_col]]

  fit <- survfit(Surv(time, status) ~ grp, data = df)
  p <- ggsurvplot(
    fit, data = df, pval = TRUE, pval.method = TRUE,
    conf.int = TRUE, risk.table = TRUE, risk.table.col = "strata",
    palette = palette, title = title, legend.labs = legend_labs,
    legend.title = "Group", xlab = "Time", ylab = "Survival Probability",
    ggtheme = theme_bw(base_size = 11), risk.table.height = 0.25,
    surv.median.line = "hv"
  )

  if (!is.null(x_max)) {
    brk <- seq(0, x_max, 365.25); lab <- 0:CUTOFF_YEARS
    p$plot  <- p$plot  + scale_x_continuous(breaks = brk, labels = lab, limits = c(0, x_max)) + xlab("Time (Years)")
    p$table <- p$table + scale_x_continuous(breaks = brk, labels = lab, limits = c(0, x_max))
  }

  # Cox HR
  cox <- coxph(Surv(time, status) ~ grp, data = df)
  hr <- exp(coef(cox)); ci <- exp(confint(cox)); cp <- summary(cox)$coefficients[, "Pr(>|z|)"]
  p$plot <- p$plot +
    annotate("text", x = 0, y = 0.05, hjust = 0, size = 3.2, color = "grey30",
             label = sprintf("HR=%.2f (%.2f-%.2f), Cox p=%.2e", hr, ci[1], ci[2], cp))
  list(plot = p, hr = hr, ci = ci, cox_p = cp)
}

# ---- Spearman  ----
plot_scatter <- function(df, xvar, yvar, xlabel, ylabel,
                         color = "#D62728", title = "") {
  ok <- is.finite(df[[xvar]]) & is.finite(df[[yvar]])
  ct <- cor.test(df[[xvar]][ok], df[[yvar]][ok], method = "spearman")
  ggplot(df[ok, ], aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(alpha = .3, size = 1.5, color = color) +
    geom_smooth(method = "lm", color = "black", linewidth = .8) +
    labs(title = title,
         subtitle = sprintf("\u03c1=%.3f  P=%.2e  (n=%d)", ct$estimate, ct$p.value, sum(ok)),
         x = xlabel, y = ylabel) +
    theme_bw(base_size = 12)
}

###############################################################################
#                                                                             #
# ██  : WGCNA   ██                           #
#                                                                             #
###############################################################################

message("\n\n", strrep("█", 60))
message(strrep("█", 60))

# Fig1A: ME Cell-Type Correlation Heatmap
message("\n--- Fig1A: ME Cell-Type Correlation Heatmap ---")

# datME / MEs
me_data <- NULL
for (vn in c("datME", "MEs", "ME_df", "moduleEigengenes")) {
  if (exists(vn)) {
    me_data <- get(vn)
    break
  }
}

# metacell metadata
mc_meta <- NULL
for (vn in c("metacell_meta", "mc_meta", "datTraits", "sample_info", "mc_info")) {
  if (exists(vn)) {
    tmp <- get(vn)
    ct_col <- intersect(c("cell_type","final_type","celltype","CellType","type"), colnames(tmp))
    if (length(ct_col) > 0) {
      mc_meta <- tmp; mc_meta$cell_type <- mc_meta[[ct_col[1]]]
      break
    }
  }
}

if (!is.null(me_data) && !is.null(mc_meta)) {
  common <- intersect(rownames(me_data), rownames(mc_meta))
  if (length(common) == 0 && nrow(me_data) == nrow(mc_meta)) {
    common <- rownames(me_data); rownames(mc_meta) <- common
  }
  if (length(common) > 10) {
    me_cols <- grep("^ME", colnames(me_data), value = TRUE)
    if (length(me_cols) == 0) me_cols <- colnames(me_data)

    me_by_ct <- aggregate(me_data[common, me_cols], by = list(CellType = mc_meta[common, "cell_type"]),
                          FUN = mean, na.rm = TRUE)
    rownames(me_by_ct) <- me_by_ct$CellType; me_by_ct$CellType <- NULL
    ct_cor <- cor(t(as.matrix(me_by_ct)), method = "pearson")

    keep <- rownames(ct_cor)[grepl("CSC|Bipotent|Hepatocyte|Hep", rownames(ct_cor), ignore.case = TRUE)]
    ct_cor_plot <- if (length(keep) >= 3) ct_cor[keep, keep] else ct_cor

    pdf(file.path(FIG_DIR, "Fig1A_ME_CellType_Correlation.pdf"), width = 10, height = 9)
    pheatmap(ct_cor_plot,
             color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
             breaks = seq(-1, 1, length.out = 101),
             display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
             main = "Cell-Type Similarity (Pearson Correlation of Mean MEs)",
             fontsize = 11)
    dev.off()
    write.csv(ct_cor_plot, file.path(FIG_DIR, "Table_Fig1A_Correlation.csv"))
    message(sprintf("   Fig1A: %d cell types", nrow(ct_cor_plot)))
  }
} else {
}

# Fig1B: Module Hallmark Enrichment (, 4  top8)
message("\n--- Fig1B: Module Hallmark Enrichment ---")

# Hallmark  (msigdbr)
if (requireNamespace("msigdbr", quietly = TRUE)) {
  library(msigdbr)
  hm_db <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_gene_list <- split(hm_db$gene_symbol, hm_db$gs_name)
} else {
  gmt_candidates <- c("/data2/reference/h.all.v2023.2.Hs.symbols.gmt",
                       "/data2/reference/h.all.v7.5.1.symbols.gmt")
  gmt_file <- gmt_candidates[file.exists(gmt_candidates)][1]
  hallmark_gene_list <- NULL
  if (!is.na(gmt_file)) {
    gmt_lines <- readLines(gmt_file)
    hallmark_gene_list <- lapply(gmt_lines, function(line) {
      parts <- strsplit(line, "\t")[[1]]; parts[-(1:2)]
    })
    names(hallmark_gene_list) <- sapply(strsplit(gmt_lines, "\t"), `[`, 1)
  }
}

if (!is.null(hallmark_gene_list)) {
  moduleColors <- setNames(gene_modules$moduleColors, rownames(gene_modules))
  all_genes_bg <- names(moduleColors)
  enrich_records <- list()

  for (mod_name in names(FOCUS_MODULES)) {
    info <- FOCUS_MODULES[[mod_name]]
    mod_gene_set <- names(moduleColors)[moduleColors == mod_name]
    if (length(mod_gene_set) < 5) next

    for (pw_name in names(hallmark_gene_list)) {
      pw_genes <- hallmark_gene_list[[pw_name]]
      overlap  <- intersect(mod_gene_set, pw_genes)
      if (length(overlap) < 2) next

      a <- length(overlap)
      b <- length(setdiff(mod_gene_set, pw_genes))
      cc <- length(setdiff(pw_genes, mod_gene_set))
      d <- length(setdiff(all_genes_bg, union(mod_gene_set, pw_genes)))
      ft <- fisher.test(matrix(c(a, b, cc, d), 2, 2), alternative = "greater")

      enrich_records[[paste(mod_name, pw_name)]] <- data.frame(
        Module = info$label, Pathway = gsub("HALLMARK_", "", pw_name),
        Overlap = a, ModSize = length(mod_gene_set), PwSize = length(pw_genes),
        OR = ft$estimate, P = ft$p.value, stringsAsFactors = FALSE)
    }
  }

  if (length(enrich_records) > 0) {
    enrich_df <- bind_rows(enrich_records) %>%
      group_by(Module) %>% mutate(P_adj = p.adjust(P, "BH")) %>% ungroup() %>%
      filter(P_adj < 0.05) %>% arrange(Module, P_adj)
    write.csv(enrich_df, file.path(FIG_DIR, "Table_Fig1B_Hallmark.csv"), row.names = FALSE)

    top_pw <- enrich_df %>% group_by(Module) %>% slice_min(P_adj, n = 8) %>% ungroup()
    top_pw$neg_log10P <- -log10(top_pw$P_adj)

    p1b <- ggplot(top_pw, aes(neg_log10P, reorder(Pathway, neg_log10P))) +
      geom_col(aes(fill = Module), alpha = .85) +
      geom_text(aes(label = sprintf("n=%d", Overlap)), hjust = -0.1, size = 3) +
      facet_wrap(~Module, scales = "free", ncol = 2) +
      labs(title = "Module Hallmark Enrichment (Fisher's Exact)",
           x = "-log10(adj.P)", y = "") +
      theme_bw(base_size = 10) + theme(legend.position = "none")

    ggsave(file.path(FIG_DIR, "Fig1B_Module_Hallmark_Enrichment.pdf"), p1b,
           width = 12, height = 8)
    message(sprintf("   Fig1B: %d modules, %d sig pathways", length(unique(top_pw$Module)), nrow(enrich_df)))
  }
} else {
}

###############################################################################
#                                                                             #

# ██  : BP  + Hoshida  ██                        #
#                                                                             #
###############################################################################

message("\n\n", strrep("█", 60))
message(strrep("█", 60))

sample_type_code <- substr(colnames(tcga_expr), 14, 15)
tumor_bc   <- colnames(tcga_expr)[sample_type_code == "01"]
adj_bc     <- colnames(tcga_expr)[sample_type_code == "11"]
message(sprintf("  TCGA samples: Tumor=%d, Adjacent=%d", length(tumor_bc), length(adj_bc)))

# signature score ()
sig_scores <- data.frame(
  barcode    = colnames(tcga_expr),
  patient_id = substr(colnames(tcga_expr), 1, 12),
  tissue     = ifelse(sample_type_code == "01", "Tumor", "Adjacent"),
  stringsAsFactors = FALSE
)

for (mod_name in names(FOCUS_MODULES)) {
  info <- FOCUS_MODULES[[mod_name]]
  sig_scores[[info$label]] <- calc_sig_score(hub_gene_lists[[mod_name]], tcga_expr)
}

# DNMT / CSC_hub / BP_hub score ()
sig_scores$DNMT_score    <- calc_sig_score(DNMT_GENES, tcga_expr)
sig_scores$CSC_hub_score <- calc_sig_score(csc_hub_genes, tcga_expr)
sig_scores$BP_hub_score  <- calc_sig_score(bp_hub_genes, tcga_expr)

# Fig2A: KM Survival — BP-M1 score in Adjacent Tissue

adj_scores <- sig_scores %>% filter(tissue == "Adjacent") %>%
  distinct(patient_id, .keep_all = TRUE)

surv_adj <- merge(adj_scores, df_os_raw[, c("patient_id","os_time","os_status")], by = "patient_id")
surv_adj <- surv_adj %>% filter(!is.na(os_time) & os_time > 0)

if (nrow(surv_adj) > 10) {
  surv_adj$group <- factor(
    ifelse(surv_adj$`BP-M1` > median(surv_adj$`BP-M1`, na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High"))

  res2a <- plot_km(surv_adj, "os_time", "os_status",
                   title = sprintf("BP-M1 Score in Adjacent Tissue (n=%d)", nrow(surv_adj)),
                   palette = c("#1F77B4", "#FF7F0E"))
  pdf(file.path(FIG_DIR, "Fig2A_KM_BP_Adjacent.pdf"), width = 3.8, height = 4.5)
  print(res2a$plot); dev.off()
  message(sprintf("   Fig2A: HR=%.2f, p=%.2e", res2a$hr, res2a$cox_p))
}

# Fig2B: Scatter — BP score vs Hoshida score
message("\n--- Fig2B: BP × Hoshida Signature ---")

# --- GSE14520  BP-M1  Hoshida  ---
calc_gse_score <- function(genes, gse_mat, probe_map) {
  if (is.null(probe_map)) {
    valid <- intersect(genes, rownames(gse_mat))
    if (length(valid) < 3) return(NULL)
    sub <- gse_mat[valid, , drop = FALSE]
  } else {
    mapped <- probe_map %>% filter(gene_clean %in% genes)
    probes <- intersect(mapped$probe, rownames(gse_mat))
    if (length(probes) < 3) return(NULL)
    sub <- gse_mat[probes, , drop = FALSE]
    sub_df <- as.data.frame(sub)
    sub_df$gene <- probe_map$gene_clean[match(rownames(sub_df), probe_map$probe)]
    sub <- sub_df %>% group_by(gene) %>% summarise(across(everything(), mean, na.rm = TRUE)) %>%
      tibble::column_to_rownames("gene") %>% as.matrix()
  }
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

if (!is.null(gse_expr) && !is.null(gse_clin)) {

  common_smp <- intersect(colnames(gse_expr), gse_clin$sample_id)
  gse_clin_sub <- gse_clin[common_smp, ]

  bp_m1_gse_score   <- calc_gse_score(hub_gene_lists[["saddlebrown"]], gse_expr[, common_smp], gse_probe_gene)
  hoshida_poor_score <- calc_gse_score(HOSHIDA_POOR, gse_expr[, common_smp], gse_probe_gene)
  hoshida_good_score <- calc_gse_score(HOSHIDA_GOOD, gse_expr[, common_smp], gse_probe_gene)

  if (!is.null(bp_m1_gse_score) && !is.null(hoshida_poor_score)) {
    gse_scores <- data.frame(
      sample_id    = common_smp,
      BP_M1        = bp_m1_gse_score[common_smp],
      Hoshida_Poor = hoshida_poor_score[common_smp],
      stringsAsFactors = FALSE
    )
    if (!is.null(hoshida_good_score)) gse_scores$Hoshida_Good <- hoshida_good_score[common_smp]

    # Fig2B: BP-M1 vs Hoshida Poor (GSE14520)
    p2b_gse <- plot_scatter(gse_scores, "BP_M1", "Hoshida_Poor",
                             "BP-M1 Signature Score", "Hoshida Poor-Prognosis Score",
                             color = "#2CA02C",
                             title = "BP-M1 vs Hoshida Poor (GSE14520, Hoshida Original Cohort)")
    ggsave(file.path(FIG_DIR, "Fig2B_BP_vs_Hoshida_GSE14520.pdf"), p2b_gse, width = 6, height = 5.5)

    gse_scores <- merge(gse_scores, gse_clin_sub[, c("sample_id","os_time","os_status")], by = "sample_id")
    gse_scores <- gse_scores %>% filter(!is.na(os_time) & os_time > 0)
    if (nrow(gse_scores) > 20) {
      gse_scores$group <- factor(
        ifelse(gse_scores$BP_M1 > median(gse_scores$BP_M1, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High"))
      res_gse_km <- plot_km(gse_scores, "os_time", "os_status",
                             title = sprintf("BP-M1 in GSE14520 Tumor (n=%d)", nrow(gse_scores)),
                             palette = c("#1F77B4", "#FF7F0E"))
      pdf(file.path(FIG_DIR, "Fig2B_KM_BP_GSE14520.pdf"), width = 3.8, height = 4.5)
      print(res_gse_km$plot); dev.off()
      message(sprintf("   GSE14520 KM: HR=%.2f, p=%.2e", res_gse_km$hr, res_gse_km$cox_p))
    }
  }
} else {
}

sig_scores$Hoshida_Poor <- calc_sig_score(HOSHIDA_POOR, tcga_expr)
sig_scores$Hoshida_Good <- calc_sig_score(HOSHIDA_GOOD, tcga_expr)

tumor_scores <- sig_scores %>% filter(tissue == "Tumor") %>% distinct(patient_id, .keep_all = TRUE)
adj_scores_h <- sig_scores %>% filter(tissue == "Adjacent") %>% distinct(patient_id, .keep_all = TRUE)

# Fig2B: BP-M1 vs Hoshida  Adjacent  ()
p2b_adj <- plot_scatter(adj_scores_h, "BP-M1", "Hoshida_Poor",
                         "BP-M1 Signature Score", "Hoshida Poor-Prognosis Score",
                         color = "#2CA02C", title = "BP-M1 vs Hoshida Poor (TCGA Adjacent)")
ggsave(file.path(FIG_DIR, "Fig2B_BP_vs_Hoshida_Adjacent.pdf"), p2b_adj, width = 6, height = 5.5)

# Tumor  ( supplementary)
p2b_tcga <- plot_scatter(tumor_scores, "BP-M1", "Hoshida_Poor",
                          "BP-M1 Signature Score", "Hoshida Poor-Prognosis Score",
                          color = "#2CA02C", title = "BP-M1 vs Hoshida Poor (TCGA Tumor)")
ggsave(file.path(FIG_DIR, "FigS_BP_vs_Hoshida_Tumor.pdf"), p2b_tcga, width = 6, height = 5.5)
message("  -> FigS_BP_vs_Hoshida_Tumor.pdf (supplementary)")

# ---- : Hoshida × BP Fisher  ----
overlap_poor <- intersect(bp_module_full, HOSHIDA_POOR)
message(sprintf("  BP Module ∩ Hoshida Poor = %d genes: %s",
                length(overlap_poor), paste(overlap_poor, collapse = ", ")))

# Venn
p_venn <- ggVennDiagram(
  list(`BP Module` = bp_module_full, `Hoshida Poor` = HOSHIDA_POOR, `Hoshida Good` = HOSHIDA_GOOD),
  label_alpha = 0, set_size = 4, label_size = 4) +
  scale_fill_gradient(low = "#F7F7F7", high = "#E41A1C", name = "Count") +
  ggtitle(sprintf("BP Module (%d genes) vs Hoshida", length(bp_module_full)))
ggsave(file.path(FIG_DIR, "Fig2B_Venn_BP_Hoshida.pdf"), p_venn, width = 8, height = 6)

###############################################################################

#                                                                             #
# ██  : DNMT × CSC   ██                     #
#                                                                             #
###############################################################################

message("\n\n", strrep("█", 60))
message(strrep("█", 60))

# Fig3A: Volcano — CSC vs Other Tumor ()
message("\n--- Fig3A: Volcano (CSC vs Other Tumor) ---")

CSC_DIFF_GENES <- NULL

if (!exists("seu_obj") && file.exists(RDATA_PATH)) {
  load(RDATA_PATH)
}

if (exists("seu_obj")) {
  DefaultAssay(seu_obj) <- "RNA"

  # Seurat v5:  RNA assay  split  layer,  JoinLayers
  if (packageVersion("Seurat") >= "5.0.0") {
    tryCatch({
      seu_obj[["RNA"]] <- JoinLayers(seu_obj[["RNA"]])
    }, error = function(e) message(sprintf("  JoinLayers skipped: %s", e$message)))
  }

  # data layer ( counts,  Normalize)
  has_data <- tryCatch({
    if (packageVersion("Seurat") >= "5.0.0") {
      d <- LayerData(seu_obj, assay = "RNA", layer = "data")
      !is.null(d) && ncol(d) > 0 && max(d[1:min(10,nrow(d)), 1:min(10,ncol(d))]) > 0
    } else {
      d <- GetAssayData(seu_obj, slot = "data")
      !is.null(d) && ncol(d) > 0
    }
  }, error = function(e) FALSE)

  if (!has_data) {
    seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
  }

  Idents(seu_obj) <- "final_type"
  all_types <- sort(unique(Idents(seu_obj)))
  csc_idents <- all_types[grepl("^CSC", all_types, ignore.case = TRUE)]

  seu_tumor_sc <- tryCatch(subset(seu_obj, malignant == "tumor"), error = function(e) NULL)
  if (is.null(seu_tumor_sc)) seu_tumor_sc <- seu_obj   # fallback: use all cells

  # Seurat v5: subset  JoinLayers
  if (packageVersion("Seurat") >= "5.0.0") {
    tryCatch({ seu_tumor_sc[["RNA"]] <- JoinLayers(seu_tumor_sc[["RNA"]]) }, error = function(e) NULL)
  }

  if (length(csc_idents) > 0) {
    Idents(seu_tumor_sc) <- "final_type"
    markers_csc <- FindMarkers(seu_tumor_sc, ident.1 = csc_idents, ident.2 = NULL,
                                logfc.threshold = 0, min.pct = 0.05, verbose = FALSE)
    markers_csc$gene <- rownames(markers_csc)

    csc_down <- markers_csc %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>% arrange(avg_log2FC)
    csc_up   <- markers_csc %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5)
    CSC_DIFF_GENES <- head(csc_down$gene, min(200, nrow(csc_down)))
    write.csv(csc_down, file.path(FIG_DIR, "Table_Fig3A_DiffGenes.csv"), row.names = FALSE)

    markers_csc$sig <- case_when(
      markers_csc$p_val_adj < 0.05 & markers_csc$avg_log2FC > 0.5  ~ "CSC Up",
      markers_csc$p_val_adj < 0.05 & markers_csc$avg_log2FC < -0.5 ~ "Differentiation",
      TRUE ~ "NS")
    markers_csc$label <- ifelse(markers_csc$gene %in% c(head(CSC_DIFF_GENES, 8), DNMT_GENES,
                                                         head(csc_up$gene, 5)),
                                markers_csc$gene, NA)

    p3a <- ggplot(markers_csc, aes(avg_log2FC, -log10(p_val_adj), color = sig)) +
      geom_point(alpha = .3, size = .8) +
      geom_point(data = markers_csc %>% filter(!is.na(label)), size = 2) +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
      scale_color_manual(values = c("CSC Up" = "#D62728", "Differentiation" = "#1F77B4", "NS" = "grey70")) +
      geom_vline(xintercept = c(-0.5, 0.5), lty = 2, color = "grey50") +
      geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey50") +
      labs(title = "CSC vs Other Tumor Cells",
           subtitle = sprintf("%d up, %d down (Differentiation genes)", nrow(csc_up), nrow(csc_down)),
           x = "log2(FC)", y = "-log10(adj.P)") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
    ggsave(file.path(FIG_DIR, "Fig3A_Volcano_CSC.pdf"), p3a, width = 8, height = 6.5)
    message(sprintf("   Fig3A: %d up, %d down, %d diff genes", nrow(csc_up), nrow(csc_down), length(CSC_DIFF_GENES)))
  }
}

if (is.null(CSC_DIFF_GENES) || length(CSC_DIFF_GENES) == 0) {
  CSC_DIFF_GENES <- c("ALB","HNF4A","CYP3A4","CYP2E1","APOB","APOA1","TF","TAT",
                       "HPD","SERPINA1","TTR","FGB","FGA","AHSG","RBP4","PCK1","G6PC")
}

# Fig3B: Scatter — DNMT score vs CSC hub score (ρ=0.559)
message("\n--- Fig3B: DNMT vs CSC Score ---")

p3b <- plot_scatter(tumor_scores, "DNMT_score", "CSC_hub_score",
                    "DNMT Score (Z-score mean)", "CSC Hub Score",
                    color = "#D62728", title = "DNMT Activity vs CSC Hub Score (Tumor)")
ggsave(file.path(FIG_DIR, "Fig3B_DNMT_vs_CSC.pdf"), p3b, width = 6, height = 5.5)

# Fig3C:  4-panel (DNMT→DiffBeta→DiffExpr, CSC→DiffBeta)
message("\n--- Fig3C: Regulation Chain ---")

diff_in <- intersect(CSC_DIFF_GENES, rownames(tcga_expr))
tumor_scores$Diff_gene_expr <- colMeans(tcga_expr[diff_in, tumor_scores$barcode, drop = FALSE], na.rm = TRUE)

if (!is.null(tcga_beta) && !is.null(probe_gene_map)) {
  diff_probes <- probe_gene_map$probe_id[probe_gene_map$gene %in% CSC_DIFF_GENES]
  diff_probes <- diff_probes[diff_probes %in% rownames(tcga_beta)]

  if (length(diff_probes) > 0) {
    expr2short  <- data.frame(expr_bc = tumor_scores$barcode, short = substr(tumor_scores$barcode, 1, 16))
    methyl2short <- data.frame(methyl_bc = colnames(tcga_beta), short = substr(colnames(tcga_beta), 1, 16))
    bc_map <- merge(expr2short, methyl2short, by = "short")

    diff_beta_vec <- colMeans(tcga_beta[diff_probes, bc_map$methyl_bc, drop = FALSE], na.rm = TRUE)
    methyl_df <- data.frame(barcode = bc_map$expr_bc, Diff_gene_beta = diff_beta_vec)
    tumor_scores <- merge(tumor_scores, methyl_df, by = "barcode", all.x = TRUE)
  }
}

chain_pairs <- list(
  list("DNMT_score",    "Diff_gene_beta", "DNMT Score",    "Diff Gene \u03b2",     "# D62728", "DNMT\u2191\u2192\u2191")
  list("Diff_gene_beta","Diff_gene_expr", "Diff Gene \u03b2",   "Diff Gene Expr",  "# 1F77B4", "\u2191\u2192\u2193")
  list("CSC_hub_score", "Diff_gene_beta", "CSC Hub Score", "Diff Gene \u03b2",     "# FF7F0E", "CSC\u2191\u2192\u2191")
  list("DNMT_score",    "Diff_gene_expr", "DNMT Score",    "Diff Gene Expr",  "# 9467BD", ": DNMT\u2191\u2192\u2193")
)

chain_plots <- list()
chain_stats <- list()

for (pair in chain_pairs) {
  xv <- pair[[1]]; yv <- pair[[2]]
  if (!xv %in% colnames(tumor_scores) || !yv %in% colnames(tumor_scores)) next
  ok <- is.finite(tumor_scores[[xv]]) & is.finite(tumor_scores[[yv]])
  if (sum(ok) < 10) next

  ct <- cor.test(tumor_scores[[xv]][ok], tumor_scores[[yv]][ok], method = "spearman")
  tag <- paste(pair[[3]], "vs", pair[[4]])
  chain_plots[[tag]] <- ggplot(tumor_scores[ok, ], aes(.data[[xv]], .data[[yv]])) +
    geom_point(alpha = .3, size = 1.5, color = pair[[5]]) +
    geom_smooth(method = "lm", color = "black", linewidth = .7) +
    labs(title = pair[[6]],
         subtitle = sprintf("\u03c1=%.3f  P=%.2e  (n=%d)", ct$estimate, ct$p.value, sum(ok)),
         x = pair[[3]], y = pair[[4]]) +
    theme_bw(base_size = 11)
  chain_stats[[tag]] <- data.frame(X = pair[[3]], Y = pair[[4]], rho = ct$estimate,
                                    p = ct$p.value, n = sum(ok), Label = pair[[6]])
  message(sprintf("  %s: ρ=%.3f P=%.2e", pair[[6]], ct$estimate, ct$p.value))
}

if (length(chain_plots) > 0) {
  pdf(file.path(FIG_DIR, "Fig3C_Regulation_Chain.pdf"), width = 12, height = 10)
  do.call(grid.arrange, c(chain_plots, ncol = 2)); dev.off()
  write.csv(bind_rows(chain_stats), file.path(FIG_DIR, "Table_Fig3C_Chain.csv"), row.names = FALSE)
}

# Fig3D: KM ×3 — DNMT alone / CSC alone / DNMT×CSC 4-group
message("\n--- Fig3D: Survival (DNMT × CSC) ---")

message(sprintf("  tumor_scores: %d rows, patient_id e.g.: %s",
                nrow(tumor_scores), head(tumor_scores$patient_id, 2)[1]))
message(sprintf("  df_os_raw: %d rows, patient_id e.g.: %s",
                nrow(df_os_raw), head(df_os_raw$patient_id, 2)[1]))
message(sprintf("  Overlap: %d patients",
                length(intersect(tumor_scores$patient_id, df_os_raw$patient_id))))

# If patient_id doesn't match (different barcode lengths), try truncating
if (length(intersect(tumor_scores$patient_id, df_os_raw$patient_id)) < 10) {
  message("  WARNING: patient_id mismatch, trying barcode truncation fix...")
  # tumor_scores patient_id should be 12-char TCGA ID
  # df_os_raw might have full barcodes
  df_os_raw$patient_id <- substr(df_os_raw$patient_id, 1, 12)
  tumor_scores$patient_id <- substr(tumor_scores$patient_id, 1, 12)
  message(sprintf("  After fix overlap: %d",
                  length(intersect(tumor_scores$patient_id, df_os_raw$patient_id))))
}

surv_df <- merge(tumor_scores, df_os_raw[, c("patient_id","os_time","os_status","age","gender",
                                               "stage_simple","grade_simple")],
                 by = "patient_id") %>%
  filter(!is.na(os_time) & os_time > 0 & !is.na(os_status)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  mutate(
    DNMT_group = ifelse(DNMT_score > median(DNMT_score, na.rm = TRUE), "DNMT_High", "DNMT_Low"),
    CSC_group  = ifelse(CSC_hub_score > median(CSC_hub_score, na.rm = TRUE), "CSC_High", "CSC_Low"),
    Combined   = factor(paste(DNMT_group, CSC_group, sep = "/"),
                        levels = c("DNMT_Low/CSC_Low","DNMT_Low/CSC_High",
                                   "DNMT_High/CSC_Low","DNMT_High/CSC_High"))
  )

if (nrow(surv_df) >= 20) {
# (i) DNMT alone
fit_d <- survfit(Surv(os_time, os_status) ~ DNMT_group, data = surv_df)
p3d1 <- ggsurvplot(fit_d, data = surv_df, pval = TRUE, risk.table = TRUE,
                    palette = c("#D62728","#1F77B4"), title = "DNMT Score (OS)",
                    ggtheme = theme_bw(base_size = 11))

# (ii) CSC alone
fit_c <- survfit(Surv(os_time, os_status) ~ CSC_group, data = surv_df)
p3d2 <- ggsurvplot(fit_c, data = surv_df, pval = TRUE, risk.table = TRUE,
                    palette = c("#D62728","#1F77B4"), title = "CSC Hub Score (OS)",
                    ggtheme = theme_bw(base_size = 11))

# (iii) 4-group
fit_4 <- survfit(Surv(os_time, os_status) ~ Combined, data = surv_df)
p3d3 <- ggsurvplot(fit_4, data = surv_df, pval = TRUE, risk.table = TRUE,
                    palette = c("#1F77B4","#FF7F0E","#2CA02C","#D62728"),
                    title = "DNMT × CSC Combined (OS)", legend.title = "Group",
                    ggtheme = theme_bw(base_size = 11))

pdf(file.path(FIG_DIR, "Fig3D_KM_DNMT_alone.pdf"), width = 3.8, height = 4.5); print(p3d1); dev.off()
pdf(file.path(FIG_DIR, "Fig3D_KM_CSC_alone.pdf"),  width = 3.8, height = 4.5); print(p3d2); dev.off()
pdf(file.path(FIG_DIR, "Fig3D_KM_Combined.pdf"),    width = 4.5, height = 5); print(p3d3); dev.off()

# Fig3E: Cox Forest Plot (DNMT + CSC + )
message("\n--- Fig3E: Cox Forest Plot ---")

tryCatch({
  # Identify which covariates actually have non-NA data
  cox_vars <- c("DNMT_score", "CSC_hub_score")
  optional_vars <- c("age", "gender", "stage_simple", "grade_simple")
  
  for (v in optional_vars) {
    if (v %in% colnames(surv_df)) {
      non_na <- sum(!is.na(surv_df[[v]]))
      n_levels <- if (is.character(surv_df[[v]]) || is.factor(surv_df[[v]]))
        length(unique(na.omit(surv_df[[v]]))) else 2
      if (non_na >= 20 && n_levels >= 2) {
        cox_vars <- c(cox_vars, v)
      } else {
        message(sprintf("  Dropping %s from Cox: %d non-NA, %d levels", v, non_na, n_levels))
      }
    }
  }
  
  # Build formula dynamically
  rhs <- paste(sapply(cox_vars, function(v) {
    if (v %in% c("stage_simple", "grade_simple")) sprintf("factor(%s)", v) else v
  }), collapse = " + ")
  cox_formula <- as.formula(paste("Surv(os_time, os_status) ~", rhs))
  message(sprintf("  Cox formula: %s", deparse(cox_formula)))
  
  # Remove rows with NA in used covariates
  cox_data <- surv_df
  for (v in cox_vars) cox_data <- cox_data[!is.na(cox_data[[v]]), ]
  message(sprintf("  Cox samples: %d (from %d)", nrow(cox_data), nrow(surv_df)))
  
  if (nrow(cox_data) < 30) stop("Too few samples for multivariate Cox")
  
  multi_cox <- coxph(cox_formula, data = cox_data)
  s <- summary(multi_cox)

  cox_df <- data.frame(
    Variable = rownames(s$conf.int),
    HR       = s$conf.int[, "exp(coef)"],
    Lower    = s$conf.int[, "lower .95"],
    Upper    = s$conf.int[, "upper .95"],
    P        = s$coefficients[, "Pr(>|z|)"]
  )
  cox_df$Variable <- gsub("factor\\(stage_simple\\)", "Stage ", cox_df$Variable)
  cox_df$Variable <- gsub("factor\\(grade_simple\\)", "Grade ", cox_df$Variable)
  cox_df$Variable <- gsub("gendermale", "Gender (Male)", cox_df$Variable)
  cox_df$Variable <- factor(cox_df$Variable, levels = rev(cox_df$Variable))
  cox_df$sig <- ifelse(cox_df$P < 0.001, "***", ifelse(cox_df$P < 0.01, "**",
                 ifelse(cox_df$P < 0.05, "*", "")))

  p3e <- ggplot(cox_df, aes(HR, Variable)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = .25, linewidth = .8) +
    geom_vline(xintercept = 1, lty = 2, color = "grey50") +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f) %s", HR, Lower, Upper, sig)),
              hjust = -0.1, size = 3) +
    scale_x_log10() +
    labs(title = sprintf("Multivariate Cox (n=%d)", nrow(cox_data)),
         x = "Hazard Ratio (log scale)", y = "") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG_DIR, "Fig3E_ForestPlot.pdf"), p3e, width = 9, height = 5)
  write.csv(cox_df, file.path(FIG_DIR, "Table_Fig3E_Cox.csv"), row.names = FALSE)
}, error = function(e) message(sprintf("   Cox error: %s", e$message)))

} else {
}

###############################################################################
#                                                                             #

###############################################################################
#                                                                             #
#    ██████  Fig4A: DotPlot — 14 epigenetic enzymes  ██████                   #
#                                                                             #
###############################################################################
message("\n", strrep("=", 70))
message(strrep("=", 70))

if (exists("seu_obj")) {
  DefaultAssay(seu_obj) <- "RNA"
  Idents(seu_obj) <- "final_type"
  tryCatch({ seu_obj[["RNA"]] <- JoinLayers(seu_obj[["RNA"]]) },
           error = function(e) NULL)
  # Fix empty data layer
  data_empty <- tryCatch({
    d <- LayerData(seu_obj, assay = "RNA", layer = "data")
    is.null(d) || all(dim(d) == 0)
  }, error = function(e) TRUE)
  if (data_empty) seu_obj <- NormalizeData(seu_obj, verbose = FALSE)

  all_types <- sort(unique(Idents(seu_obj)))
  csc_ids <- all_types[grepl("^CSC", all_types)]
  bp_ids  <- all_types[grepl("^Bipotent", all_types)]
  hep_ids <- all_types[grepl("^Hepatocyte", all_types)]
  key_types <- c(bp_ids, hep_ids, csc_ids)

  seu_sub <- subset(seu_obj, idents = key_types)
  tryCatch({ seu_sub[["RNA"]] <- JoinLayers(seu_sub[["RNA"]]) }, error = function(e) NULL)
  seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
  epi_in <- intersect(EPI_GENES, rownames(seu_sub))

  p4a <- DotPlot(seu_sub, features = epi_in, group.by = "final_type") +
    RotatedAxis() +
    labs(title = "Epigenetic Enzymes: CSC vs Bipotent vs Hepatocytes") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))
  ggsave(file.path(FIG_DIR, "Fig4A_DotPlot_EpiEnzymes.pdf"), p4a,
         width = max(8, length(epi_in) * 0.7), height = max(5, length(key_types) * 0.5))
  message("  -> Fig4A_DotPlot_EpiEnzymes.pdf")

  # --- FigS8: Heatmap ---
  fd <- tryCatch(FetchData(seu_sub, vars = epi_in, layer = "data"),
                 error = function(e) tryCatch(FetchData(seu_sub, vars = epi_in, layer = "counts"),
                                               error = function(e2) NULL))
  if (!is.null(fd)) {
    fd$cell_type <- as.character(Idents(seu_sub))
    hm_data <- fd %>% group_by(cell_type) %>%
      summarise(across(all_of(epi_in), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      tibble::column_to_rownames("cell_type")
    hm_z <- scale(as.matrix(hm_data))
    hm_z[is.na(hm_z)] <- 0; hm_z[is.infinite(hm_z)] <- 0
    hm_z <- pmax(pmin(hm_z, 2), -2)

    gene_anno <- data.frame(
      Category = case_when(
        colnames(hm_z) %in% c("DNMT1","DNMT3A","DNMT3B","UHRF1") ~ "Methylation",
        colnames(hm_z) %in% c("TET1","TET2","TET3","TDG","GADD45A","GADD45B") ~ "Demethylation",
        TRUE ~ "Chromatin"), row.names = colnames(hm_z))

    pdf(file.path(FIG_DIR, "FigS8_Heatmap_EpiEnzymes.pdf"), width = 10, height = 6)
    pheatmap(hm_z, color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
             annotation_col = gene_anno,
             annotation_colors = list(Category = c(Methylation="#D62728",
                                                    Demethylation="#2CA02C", Chromatin="#9467BD")),
             main = "Epigenetic Enzymes (Z-score)", fontsize = 10)
    dev.off()
    message("  -> FigS8_Heatmap_EpiEnzymes.pdf")

    # --- FigS9: Boxplot ---
    fd$Epi_total <- rowSums(fd[, epi_in, drop = FALSE])
    type_order <- fd %>% group_by(cell_type) %>%
      summarise(med = median(Epi_total), .groups = "drop") %>%
      arrange(med) %>% pull(cell_type)
    fd$cell_type <- factor(fd$cell_type, levels = type_order)
    type_colors <- setNames(
      c(rep("#1F77B4", length(bp_ids)), rep("#2CA02C", length(hep_ids)), rep("#D62728", length(csc_ids))),
      c(as.character(bp_ids), as.character(hep_ids), as.character(csc_ids)))

    pS9 <- ggplot(fd, aes(cell_type, Epi_total, fill = cell_type)) +
      geom_boxplot(outlier.shape = NA, alpha = .8) +
      stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
      labs(title = "Total Epigenetic Enzyme Activity", y = "Sum of 14 genes", x = "") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      scale_fill_manual(values = type_colors)
    ggsave(file.path(FIG_DIR, "FigS9_Boxplot_EpiActivity.pdf"), pS9, width = 10, height = 6)
    message("  -> FigS9_Boxplot_EpiActivity.pdf")
  }

  # --- BP vs Hep DEG (for Fig4C later) ---
  message("\n  Computing BP vs Hep DEG...")
  bp_type <- "Bipotent Progenitor"
  if (!bp_type %in% levels(Idents(seu_obj))) {
    bp_cand <- grep("^Bipotent|^BP", levels(Idents(seu_obj)), value = TRUE, ignore.case = TRUE)
    if (length(bp_cand) > 0) bp_type <- bp_cand[1]
  }
  hep_types <- grep("^Hepatocytes", levels(Idents(seu_obj)), value = TRUE)

  bp_down_genes <- NULL
  if (bp_type %in% levels(Idents(seu_obj)) && length(hep_types) > 0) {
    sub2 <- subset(seu_obj, idents = c(bp_type, hep_types))
    tryCatch({ sub2[["RNA"]] <- JoinLayers(sub2[["RNA"]]) }, error = function(e) NULL)
    sub2 <- NormalizeData(sub2, verbose = FALSE)
    sub2$cmp <- ifelse(sub2$final_type == bp_type, "BP", "Hep")
    Idents(sub2) <- "cmp"
    message(sprintf("  BP=%d Hep=%d", sum(Idents(sub2)=="BP"), sum(Idents(sub2)=="Hep")))

    mk <- FindMarkers(sub2, ident.1 = "BP", ident.2 = "Hep",
                       min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox", only.pos = FALSE)
    mk$gene <- rownames(mk)
    bp_down_genes <- mk$gene[mk$avg_log2FC < -0.25 & mk$p_val_adj < 0.05]
    message(sprintf("  BP-down genes: %d", length(bp_down_genes)))
    write.csv(mk, file.path(FIG_DIR, "Table_BP_vs_Hep_DEG.csv"), row.names = FALSE)
    rm(sub2); gc(verbose = FALSE)
  }
  rm(seu_sub); gc(verbose = FALSE)
} else {
  message("  WARNING: seu_obj not available")
  bp_down_genes <- module_gene_lists[["BP-M1"]]
}

###############################################################################
#                                                                             #
# ██████  Fig4B: BP-M1 TF Enrichment ( DNMT1)  ██████                 #
#                                                                             #
###############################################################################
message("\n", strrep("=", 70))
# Helper: label to color name
if (!exists("label_to_color")) {
  label_to_color <- function(lab)
    names(FOCUS_MODULES)[sapply(FOCUS_MODULES, function(x) x$label == lab)]
}

# Build module_gene_lists if not yet defined
if (!exists("module_gene_lists")) {
  module_gene_lists <- list()
  for (m in names(FOCUS_MODULES))
    module_gene_lists[[ FOCUS_MODULES[[m]]$label ]] <-
      rownames(gene_modules)[gene_modules$moduleColors == m]
  all_background <- unique(rownames(gene_modules))
  message(sprintf("  Built module_gene_lists: %s",
    paste(sapply(names(module_gene_lists), function(x) sprintf("%s=%d", x, length(module_gene_lists[[x]]))), collapse=", ")))
}

message("  Fig4B: BP-M1 TF Enrichment")
message(strrep("=", 70))

# MSigDB C3 TFT (GTRD + Legacy) 4TF
c3_gtrd   <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
  dplyr::select(gs_name, gene_symbol)
c3_legacy <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>%
  dplyr::select(gs_name, gene_symbol)
c3_all <- bind_rows(c3_gtrd, c3_legacy) %>% distinct()

enr_all <- list()
for (ml in names(module_gene_lists)) {
  enr <- tryCatch(enricher(gene = module_gene_lists[[ml]], TERM2GENE = c3_all,
                            pvalueCutoff = 0.2, qvalueCutoff = 0.5,
                            minGSSize = 10, maxGSSize = 3000), error = function(e) NULL)
  if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
    df <- as.data.frame(enr); df$Module <- ml
    n_dnmt <- sum(grepl("DNMT", df$ID, ignore.case = TRUE))
    message(sprintf("  %s: %d sig, %d DNMT-related", ml, sum(df$qvalue < 0.05), n_dnmt))
    enr_all[[ml]] <- df
  } else {
    message(sprintf("  %s: no results", ml))
  }
}
if (length(enr_all) > 0)
  write.csv(bind_rows(enr_all), file.path(FIG_DIR, "Table_C3TFT_AllModules.csv"), row.names = FALSE)

# Fig4B: BP-M1  TF  barplot, DNMT
if ("BP-M1" %in% names(enr_all)) {
  top_bp <- enr_all[["BP-M1"]] %>% arrange(qvalue) %>% head(15)
  top_bp$is_dnmt <- grepl("DNMT", top_bp$ID, ignore.case = TRUE)
  top_bp$Term <- gsub("^.*?_", "", top_bp$ID)
  top_bp$Term <- make.unique(top_bp$Term, sep = " ")
  top_bp$Term <- factor(top_bp$Term, levels = rev(top_bp$Term))

  p4b <- ggplot(top_bp, aes(x = -log10(qvalue), y = Term)) +
    geom_col(aes(fill = is_dnmt), width = 0.7, show.legend = FALSE) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("FALSE" = "#4DBBD5", "TRUE" = "#D62728")) +
    labs(title = "BP-M1 Upstream Regulators (C3 TFT)",
         subtitle = "Red = DNMT-related | Dashed = q = 0.05",
         x = "-log10(q-value)", y = "") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG_DIR, "Fig4B_TF_Enrichment_BP_M1.pdf"), p4b, width = 9, height = 5)
  message("  -> Fig4B_TF_Enrichment_BP_M1.pdf")
}

# --- FigS4: CSC-M1  C3 TFT barplot ---
if ("CSC-M1" %in% names(enr_all)) {
  top_csc <- enr_all[["CSC-M1"]] %>% arrange(qvalue) %>% head(20)
  top_csc$is_dnmt <- grepl("DNMT", top_csc$ID, ignore.case = TRUE)
  top_csc$Term <- gsub("^.*?_", "", top_csc$ID)
  top_csc$Term <- make.unique(top_csc$Term, sep = " ")
  top_csc$Term <- factor(top_csc$Term, levels = rev(top_csc$Term))

  pS4 <- ggplot(top_csc, aes(x = -log10(qvalue), y = Term)) +
    geom_col(aes(fill = is_dnmt), width = 0.7, show.legend = FALSE) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("FALSE" = "#4DBBD5", "TRUE" = "#D62728")) +
    labs(title = "CSC-M1 Upstream Regulators (C3 TFT)",
         subtitle = "Red = DNMT-related", x = "-log10(q-value)", y = "") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG_DIR, "FigS4_C3TFT_CSC_M1.pdf"), pS4, width = 9, height = 6)
  message("  -> FigS4_C3TFT_CSC_M1.pdf")
}

# Ensure tumor/adjacent log-expr matrices exist for Fig4C/D and supplementary
if (!exists("exp_tumor")) {
  if (!exists("tumor_bc")) {
    stype <- substr(colnames(tcga_expr), 14, 15)
    tumor_bc <- colnames(tcga_expr)[stype == "01"]
    adj_bc   <- colnames(tcga_expr)[stype == "11"]
  }
  exp_tumor   <- log2(tcga_expr[, tumor_bc, drop = FALSE] + 1)
  adj_exp_log <- log2(tcga_expr[, adj_bc,   drop = FALSE] + 1)
  message(sprintf("  exp_tumor: %d x %d, adj_exp_log: %d x %d",
                  nrow(exp_tumor), ncol(exp_tumor), nrow(adj_exp_log), ncol(adj_exp_log)))
}
if (!exists("adj_exp_log")) {
  adj_exp_log <- log2(tcga_expr[, adj_bc, drop = FALSE] + 1)
}

# Validate probe_gene_map column names
if (!is.null(probe_gene_map)) {
  if (!"probe_id" %in% colnames(probe_gene_map)) {
    # Try common alternatives
    cand_probe <- intersect(c("IlmnID", "ID", "Name", "probe", "Probe_ID"), colnames(probe_gene_map))
    cand_gene  <- intersect(c("Gene_Symbol", "gene_symbol", "GENE", "Gene", "UCSC_RefGene_Name"), colnames(probe_gene_map))
    if (length(cand_probe) > 0) {
      message(sprintf("  Renaming probe_gene_map: %s -> probe_id", cand_probe[1]))
      colnames(probe_gene_map)[colnames(probe_gene_map) == cand_probe[1]] <- "probe_id"
    } else if (ncol(probe_gene_map) >= 2) {
      message("  Using first two columns of probe_gene_map as probe_id and gene")
      colnames(probe_gene_map)[1:2] <- c("probe_id", "gene")
    }
    if (length(cand_gene) > 0 && !"gene" %in% colnames(probe_gene_map)) {
      colnames(probe_gene_map)[colnames(probe_gene_map) == cand_gene[1]] <- "gene"
    }
  }
  message(sprintf("  probe_gene_map columns: %s", paste(colnames(probe_gene_map), collapse = ", ")))
}

###############################################################################
#                                                                             #
#    ██████  Fig4C: BP-down beta-Expr Histogram  ██████                       #
#                                                                             #
###############################################################################
message("\n", strrep("=", 70))
message("  Fig4C: BP-down genes beta-Expr (Adjacent)")
message(strrep("=", 70))

message(sprintf("  DEBUG: bp_down_genes=%s, tcga_beta=%s, probe_gene_map=%s",
  if(is.null(bp_down_genes)) "NULL" else as.character(length(bp_down_genes)),
  if(is.null(tcga_beta)) "NULL" else paste(dim(tcga_beta), collapse="x"),
  if(is.null(probe_gene_map)) "NULL" else paste(dim(probe_gene_map), collapse="x")))

if (!is.null(bp_down_genes) && length(bp_down_genes) > 5 &&
    !is.null(tcga_beta) && !is.null(probe_gene_map)) {

  adj_bc_m   <- colnames(tcga_beta)[substr(colnames(tcga_beta), 14, 15) == "11"]
  common_adj <- intersect(substr(adj_bc_m, 1, 12), substr(adj_bc, 1, 12))
  message(sprintf("  Paired adjacent: %d", length(common_adj)))

  bp_cor <- list()
  for (gene in bp_down_genes) {
    if (!gene %in% rownames(tcga_expr)) next
    pa <- probe_gene_map$probe_id[probe_gene_map$gene == gene]
    pa <- pa[pa %in% rownames(tcga_beta)]; if (length(pa) == 0) next
    mb <- adj_bc_m[substr(adj_bc_m, 1, 12) %in% common_adj]; if (length(mb) < 10) next
    mv <- if (length(pa) == 1) as.numeric(tcga_beta[pa, mb]) else
      colMeans(tcga_beta[pa, mb, drop = FALSE], na.rm = TRUE)
    dm <- data.frame(pid = substr(mb, 1, 12), beta = mv) %>% distinct(pid, .keep_all = TRUE)
    eb <- adj_bc[substr(adj_bc, 1, 12) %in% common_adj]; if (length(eb) < 10) next
    de <- data.frame(pid = substr(eb, 1, 12), expr = as.numeric(tcga_expr[gene, eb])) %>%
      distinct(pid, .keep_all = TRUE)
    dp <- inner_join(dm, de, by = "pid") %>% filter(is.finite(beta), is.finite(expr))
    if (nrow(dp) < 10) next
    ct <- tryCatch(cor.test(dp$beta, dp$expr, method = "spearman"), error = function(e) NULL)
    if (!is.null(ct))
      bp_cor[[gene]] <- data.frame(Gene = gene, rho = ct$estimate, p = ct$p.value, n = nrow(dp),
                                    mean_beta = mean(dp$beta), mean_expr = mean(dp$expr))
  }

  if (length(bp_cor) > 0) {
    df_cor <- bind_rows(bp_cor) %>% mutate(fdr = p.adjust(p, "BH")) %>% arrange(rho)
    n_neg      <- sum(df_cor$rho < 0 & df_cor$fdr < 0.05)
    median_rho <- median(df_cor$rho)
    overall_t  <- t.test(df_cor$rho, mu = 0, alternative = "less")
    message(sprintf("  %d tested, %d sig.neg, median rho=%.3f P=%.2e",
                    nrow(df_cor), n_neg, median_rho, overall_t$p.value))
    write.csv(df_cor, file.path(FIG_DIR, "Table_BP_Down_BetaExprCor.csv"), row.names = FALSE)

    p4c <- ggplot(df_cor, aes(x = rho)) +
      geom_histogram(bins = 30, fill = "#9467BD", alpha = 0.75, color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = median_rho, color = "#D62728", linewidth = 0.8) +
      annotate("text", x = median_rho, y = Inf, vjust = 2,
               label = sprintf("median rho=%.3f\nP=%.2e", median_rho, overall_t$p.value),
               color = "#D62728", fontface = "bold", size = 3.5) +
      labs(title = "BP-Downregulated Genes: beta vs Expr in Adjacent",
           subtitle = sprintf("n=%d | %d sig. negative (FDR<0.05)", nrow(df_cor), n_neg),
           x = "Spearman rho (promoter beta vs expression)", y = "Count") +
      theme_bw(base_size = 11)
    ggsave(file.path(FIG_DIR, "Fig4C_BP_Down_BetaExpr_Histogram.pdf"), p4c, width = 6, height = 4.5)
    message("  -> Fig4C_BP_Down_BetaExpr_Histogram.pdf")

    # --- FigS6: Top6 scatter ---
    top6 <- df_cor %>% filter(rho < 0) %>% arrange(p) %>% head(6)
    if (nrow(top6) > 0) {
      slist <- list()
      for (i in seq_len(nrow(top6))) {
        g <- top6$Gene[i]
        pa <- probe_gene_map$probe_id[probe_gene_map$gene == g]
        pa <- pa[pa %in% rownames(tcga_beta)]
        mb2 <- adj_bc_m[substr(adj_bc_m, 1, 12) %in% common_adj]
        mv2 <- if (length(pa)==1) as.numeric(tcga_beta[pa,mb2]) else
          colMeans(tcga_beta[pa,mb2,drop=FALSE],na.rm=TRUE)
        dm2 <- data.frame(pid=substr(mb2,1,12),beta=mv2) %>% distinct(pid,.keep_all=TRUE)
        eb2 <- adj_bc[substr(adj_bc,1,12) %in% common_adj]
        de2 <- data.frame(pid=substr(eb2,1,12),expr=as.numeric(tcga_expr[g,eb2])) %>%
          distinct(pid,.keep_all=TRUE)
        dp2 <- inner_join(dm2,de2,by="pid") %>% filter(is.finite(beta),is.finite(expr))
        ct2 <- tryCatch(cor.test(dp2$beta,dp2$expr,method="spearman"),error=function(e)NULL)
        if (!is.null(ct2))
          slist[[g]] <- ggplot(dp2,aes(beta,expr)) +
            geom_point(alpha=0.5,size=2,color="#9467BD") +
            geom_smooth(method="lm",color="black",linewidth=0.7,formula=y~x) +
            labs(title=g,subtitle=sprintf("rho=%.3f P=%.2e",ct2$estimate,ct2$p.value),
                 x="Promoter beta",y="Expression") + theme_bw(base_size=9)
      }
      if (length(slist)>0) {
        nc <- min(3,length(slist)); nr <- ceiling(length(slist)/nc)
        pdf(file.path(FIG_DIR,"FigS6_BP_Down_TopScatters.pdf"),width=4*nc,height=4*nr)
        do.call(grid.arrange, c(slist, ncol=nc)); dev.off()
        message("  -> FigS6_BP_Down_TopScatters.pdf")
      }
    }
  }
} else {
  message("  SKIPPED: bp_down_genes or tcga_beta missing")
}

###############################################################################
#                                                                             #
#    ██████  Fig4D: Dual Pathway Dot Plot  ██████                             #
#                                                                             #
###############################################################################
message("\n", strrep("=", 70))
message(strrep("=", 70))

dnmt_avail <- DNMT_GENES[DNMT_GENES %in% rownames(exp_tumor)]
dual <- list()

for (ml in c("CSC-M1","BP-M1")) {
  cn <- label_to_color(ml)
  for (tissue in c("Tumor","Adjacent")) {
    em <- if (tissue == "Tumor") exp_tumor else adj_exp_log
    mg <- intersect(hub_gene_lists[[cn]], rownames(em))
    if (length(mg) < 5) next
    ms <- colMeans(t(scale(t(em[mg,,drop=FALSE]))), na.rm = TRUE)
    for (d in dnmt_avail) {
      dv <- as.numeric(em[d,]); ok <- is.finite(ms) & is.finite(dv)
      if (sum(ok) < 15) next
      ct <- cor.test(ms[ok], dv[ok], method = "spearman")
      dual[[paste(ml,tissue,d)]] <- data.frame(
        Module=ml, Tissue=tissue, DNMT=d, rho=ct$estimate, p=ct$p.value, n=sum(ok))
      sig <- if(ct$p.value<0.001)"***" else if(ct$p.value<0.05)"*" else ""
      message(sprintf("    %s x %s (%s): rho=%.3f P=%.2e %s",
                      ml, d, tissue, ct$estimate, ct$p.value, sig))
    }
  }
}

if (length(dual) > 0) {
  df_dual <- bind_rows(dual) %>%
    mutate(sig = p < 0.05, Panel = paste(Module, Tissue, sep = "\n"))
  df_dual$Panel <- factor(df_dual$Panel,
    levels = c("CSC-M1\nTumor","CSC-M1\nAdjacent","BP-M1\nTumor","BP-M1\nAdjacent"))
  write.csv(df_dual, file.path(FIG_DIR, "Table_DualPathway_Compare.csv"), row.names = FALSE)

  p4d <- ggplot(df_dual, aes(x = rho, y = Panel, color = DNMT, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 4, position = position_dodge(0.5)) +
    scale_color_manual(values = c(DNMT1="#D62728", DNMT3A="#FF7F0E", DNMT3B="#9467BD")) +
    scale_shape_manual(values = c("TRUE"=16, "FALSE"=1), guide = "none") +
    labs(title = "DNMT Dual Pathway: CSC (Tumor) vs BP (Adjacent)",
         subtitle = "Filled = P<0.05 | Open = NS",
         x = "Spearman rho (DNMT x Module activity)", y = "") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, "Fig4D_DualPathway_Compare.pdf"), p4d, width = 7, height = 4.5)
  message("  -> Fig4D_DualPathway_Compare.pdf")
}

###############################################################################
#                                                                             #
#    ██████  Supplementary: FigS1-S5, S7  ██████                             #
#                                                                             #
###############################################################################
message("\n", strrep("=", 70))
message("  Supplementary Figures")
message(strrep("=", 70))

# --- FigS1: DNMT x CSC-M1 Hub heatmap ---
message("\n--- FigS1: DNMT x Hub heatmap ---")
hub_avail <- hub_gene_lists[["darkturquoise"]]
hub_avail <- hub_avail[hub_avail %in% rownames(exp_tumor)]

cor_mat  <- matrix(NA, length(hub_avail), length(dnmt_avail),
                   dimnames = list(hub_avail, dnmt_avail))
pval_mat <- cor_mat
for (d in dnmt_avail) for (h in hub_avail) {
  ct <- cor.test(as.numeric(exp_tumor[d,]), as.numeric(exp_tumor[h,]), method="spearman")
  cor_mat[h,d] <- ct$estimate; pval_mat[h,d] <- ct$p.value
}
ss <- matrix("", nrow(pval_mat), ncol(pval_mat))
ss[pval_mat < 0.001] <- "***"; ss[pval_mat >= 0.001 & pval_mat < 0.01] <- "**"
ss[pval_mat >= 0.01 & pval_mat < 0.05] <- "*"
ord <- order(cor_mat[,"DNMT1"], decreasing = TRUE)
cor_mat <- cor_mat[ord,,drop=FALSE]; ss <- ss[ord,,drop=FALSE]

pdf(file.path(FIG_DIR, "FigS1_Heatmap_DNMT_vs_Hub.pdf"), width = 4.5, height = 8)
pheatmap(cor_mat, display_numbers = ss,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         breaks = seq(-0.6, 0.6, length.out = 101),
         cluster_cols = FALSE, cluster_rows = FALSE,
         main = "DNMT x CSC-M1 Hub Genes (Spearman rho)",
         fontsize = 9, fontsize_number = 11, angle_col = 0,
         cellwidth = 40, cellheight = 14)
dev.off()
write.csv(cor_mat, file.path(FIG_DIR, "Table_DNMT_Hub_Cor.csv"))
message("  -> FigS1_Heatmap_DNMT_vs_Hub.pdf")

# --- FigS2: Fisher DNMT enrichment ---
message("\n--- FigS2: Fisher test ---")
background_avail <- intersect(all_background, rownames(exp_tumor))
fisher_res <- list()
for (d in dnmt_avail) {
  dv <- as.numeric(exp_tumor[d,])
  cors <- sapply(background_avail, function(g) {
    ct <- tryCatch(cor.test(dv, as.numeric(exp_tumor[g,]), method="spearman"), error=function(e) NULL)
    if (is.null(ct)) return(c(NA, NA)); c(unname(ct$estimate), ct$p.value)
  })
  dfc <- data.frame(gene=background_avail, r=cors[1,], p=cors[2,]) %>%
    filter(is.finite(r)) %>% mutate(fdr=p.adjust(p,"BH"))
  dpos <- dfc$gene[dfc$r > 0.3 & dfc$fdr < 0.05]
  message(sprintf("  %s pos-cor genes: %d", d, length(dpos)))
  for (ml in names(module_gene_lists)) {
    mg <- intersect(module_gene_lists[[ml]], background_avail)
    a <- length(intersect(mg,dpos)); b <- length(setdiff(mg,dpos))
    cv <- length(setdiff(dpos,mg)); dv2 <- length(background_avail)-a-b-cv
    ft <- fisher.test(matrix(c(a,b,cv,dv2),2), alternative="greater")
    fisher_res[[paste(ml,d)]] <- data.frame(
      Module=ml, DNMT=d, Overlap=a,
      Expected=round(length(mg)*length(dpos)/length(background_avail),1),
      FoldEnrich=round(a/max(1,length(mg)*length(dpos)/length(background_avail)),2),
      OddsRatio=round(ft$estimate,3), P_value=ft$p.value)
  }
}
dff <- bind_rows(fisher_res)
write.csv(dff, file.path(FIG_DIR, "Table_Fisher_DNMT_Overlap.csv"), row.names = FALSE)
dff$Module <- factor(dff$Module, levels=c("CSC-M1","CSC-M2","BP-M1","BP-M2"))
dff$sig <- dff$P_value < 0.05

pS2 <- ggplot(dff, aes(x=OddsRatio, y=Module, fill=DNMT)) +
  geom_col(position=position_dodge(0.7), width=0.6, alpha=0.85) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_text(aes(label=ifelse(sig, sprintf("%.1fx",FoldEnrich),""), x=OddsRatio+0.1),
            position=position_dodge(0.7), size=3, hjust=0) +
  scale_fill_manual(values=c(DNMT1="#D62728",DNMT3A="#FF7F0E",DNMT3B="#9467BD")) +
  labs(title="DNMT-Correlated Genes Enrichment in Modules",
       subtitle="Fisher test | Labels=fold (P<0.05)", x="Odds Ratio", y="") +
  theme_bw(base_size=12)
ggsave(file.path(FIG_DIR, "FigS2_Fisher_DNMT_Enrichment.pdf"), pS2, width=8, height=4)
message("  -> FigS2_Fisher_DNMT_Enrichment.pdf")

# --- FigS3: 4-module x DNMT dot plot ---
message("\n--- FigS3: 4-module x DNMT ---")
mod_res <- list()
for (ml in names(module_gene_lists)) {
  g <- intersect(module_gene_lists[[ml]], rownames(exp_tumor))
  if (length(g)<3) next
  ms <- colMeans(t(scale(t(exp_tumor[g,,drop=FALSE]))), na.rm=TRUE)
  for (d in dnmt_avail) {
    dv <- as.numeric(exp_tumor[d,]); ok <- is.finite(ms)&is.finite(dv)
    ct <- cor.test(ms[ok],dv[ok],method="spearman")
    mod_res[[paste(ml,d)]] <- data.frame(Module=ml,DNMT=d,rho=ct$estimate,p=ct$p.value)
  }
}
dfm <- bind_rows(mod_res)
write.csv(dfm, file.path(FIG_DIR, "Table_Module_DNMT_Cor.csv"), row.names=FALSE)
dfm$Module <- factor(dfm$Module, levels=c("CSC-M1","CSC-M2","BP-M1","BP-M2"))
dfm$sig <- dfm$p < 0.05

pS3 <- ggplot(dfm, aes(x=rho, y=Module)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(color=DNMT, shape=sig), size=4, position=position_dodge(0.5)) +
  scale_color_manual(values=c(DNMT1="#D62728",DNMT3A="#FF7F0E",DNMT3B="#9467BD")) +
  scale_shape_manual(values=c("TRUE"=16,"FALSE"=1), guide="none") +
  labs(title="DNMT x Module Activity (TCGA Tumor)",
       subtitle="Filled=P<0.05", x="Spearman rho", y="") +
  theme_bw(base_size=12) + theme(legend.position="right")
ggsave(file.path(FIG_DIR, "FigS3_DNMT_4Module_Compare.pdf"), pS3, width=7, height=3.5)
message("  -> FigS3_DNMT_4Module_Compare.pdf")

# --- FigS5: DNMT x CSC-M1 scatter trio ---
message("\n--- FigS5: DNMT x CSC-M1 scatter ---")
sp <- list()
for (d in dnmt_avail) {
  g <- intersect(module_gene_lists[["CSC-M1"]], rownames(exp_tumor))
  ms <- colMeans(t(scale(t(exp_tumor[g,,drop=FALSE]))), na.rm=TRUE)
  dv <- as.numeric(exp_tumor[d,]); ok <- is.finite(ms)&is.finite(dv)
  ct <- cor.test(ms[ok],dv[ok],method="spearman")
  sdf <- data.frame(DNMT=dv[ok], Score=ms[ok])
  sp[[d]] <- ggplot(sdf, aes(DNMT,Score)) +
    geom_point(alpha=0.25, size=1.2, color="#D62728") +
    geom_smooth(method="lm", color="black", linewidth=0.7, formula=y~x) +
    labs(title=sprintf("%s vs CSC-M1",d),
         subtitle=sprintf("rho=%.3f P=%.2e",ct$estimate,ct$p.value),
         x=sprintf("log2(%s+1)",d), y="CSC-M1 Score") +
    theme_bw(base_size=10)
}
if (length(sp)>0) {
  pdf(file.path(FIG_DIR,"FigS5_Scatter_DNMT_vs_CSC_M1.pdf"),width=4*length(sp),height=4)
  do.call(grid.arrange, c(sp, nrow=1)); dev.off()
  message("  -> FigS5_Scatter_DNMT_vs_CSC_M1.pdf")
}

# --- FigS7: DNMT1 x BP-M1 Adjacent scatter ---
message("\n--- FigS7: DNMT1 x BP-M1 Adjacent ---")
bp_mg <- intersect(hub_gene_lists[["saddlebrown"]], rownames(adj_exp_log))
if (length(bp_mg)>=5 && "DNMT1" %in% rownames(adj_exp_log)) {
  z <- t(scale(t(adj_exp_log[bp_mg,,drop=FALSE])))
  bps <- colMeans(z,na.rm=TRUE); d1 <- as.numeric(adj_exp_log["DNMT1",])
  ok <- is.finite(bps)&is.finite(d1)
  ct <- cor.test(bps[ok],d1[ok],method="spearman")
  sdf <- data.frame(DNMT1=d1[ok], BP_M1=bps[ok])
  pS7 <- ggplot(sdf,aes(DNMT1,BP_M1)) +
    geom_point(alpha=0.5,size=2,color="#1F77B4") +
    geom_smooth(method="lm",color="black",linewidth=0.7,formula=y~x) +
    labs(title="DNMT1 x BP-M1 Score (Adjacent)",
         subtitle=sprintf("rho=%.3f P=%.2e (n=%d)",ct$estimate,ct$p.value,sum(ok)),
         x="log2(DNMT1+1)", y="BP-M1 Module Score") +
    theme_bw(base_size=11)
  ggsave(file.path(FIG_DIR,"FigS7_DNMT1_vs_BP_M1_Adjacent.pdf"),pS7,width=5.5,height=5)
  message("  -> FigS7_DNMT1_vs_BP_M1_Adjacent.pdf")
}

# --- BP-M1 Hub beta-Expr ---
if (!is.null(tcga_beta) && !is.null(probe_gene_map)) {
  message("\n--- Table: BP-M1 Hub beta-Expr ---")
  hcor <- list()
  abm <- colnames(tcga_beta)[substr(colnames(tcga_beta),14,15)=="11"]
  cadj <- intersect(substr(abm,1,12), substr(adj_bc,1,12))
  for (gene in hub_gene_lists[["saddlebrown"]]) {
    if (!gene %in% rownames(tcga_expr)) next
    pa <- probe_gene_map$probe_id[probe_gene_map$gene==gene]
    pa <- pa[pa %in% rownames(tcga_beta)]; if(length(pa)==0) next
    mb <- abm[substr(abm,1,12) %in% cadj]; if(length(mb)<10) next
    mv <- if(length(pa)==1) as.numeric(tcga_beta[pa,mb]) else colMeans(tcga_beta[pa,mb,drop=FALSE],na.rm=TRUE)
    dm <- data.frame(pid=substr(mb,1,12),beta=mv) %>% distinct(pid,.keep_all=TRUE)
    eb <- adj_bc[substr(adj_bc,1,12) %in% cadj]; if(length(eb)<10) next
    de <- data.frame(pid=substr(eb,1,12),expr=as.numeric(tcga_expr[gene,eb])) %>% distinct(pid,.keep_all=TRUE)
    dp <- inner_join(dm,de,by="pid") %>% filter(is.finite(beta),is.finite(expr))
    if(nrow(dp)<10) next
    ct <- tryCatch(cor.test(dp$beta,dp$expr,method="spearman"),error=function(e)NULL)
    if(!is.null(ct)) hcor[[gene]] <- data.frame(Gene=gene,rho=ct$estimate,p=ct$p.value,n=nrow(dp))
  }
  if(length(hcor)>0) {
    dfh <- bind_rows(hcor) %>% mutate(fdr=p.adjust(p,"BH")) %>% arrange(rho)
    write.csv(dfh, file.path(FIG_DIR,"Table_BP_Hub_BetaExprCor.csv"), row.names=FALSE)
    message(sprintf("  BP-M1 Hub: %d genes, median rho=%.3f", nrow(dfh), median(dfh$rho)))
  }
}

###############################################################################
#                         DONE                                                #
###############################################################################
message("\n", strrep("=", 70))
message("  ALL DONE!")
message(strrep("=", 70))
message(sprintf("\nOutput: %s\n", FIG_DIR))
message("
=== Main figures (Paper Figure 4) ===
  Fig4A_DotPlot_EpiEnzymes.pdf        14 epigenetic enzymes DotPlot
  Fig4B_TF_Enrichment_BP_M1.pdf       BP-M1 TF enrichment (DNMT1 highlighted)
  Fig4C_BP_Down_BetaExpr_Histogram.pdf BP downregulated genes beta-Expr histogram
  Fig4D_DualPathway_Compare.pdf        dual pathway comparison dotplot

=== Supplementary figures ===
  FigS1_Heatmap_DNMT_vs_Hub.pdf        DNMT x 30 CSC-M1 hub heatmap
  FigS2_Fisher_DNMT_Enrichment.pdf     Fisher OR barplot
  FigS3_DNMT_4Module_Compare.pdf       4-module x DNMT dot
  FigS4_C3TFT_CSC_M1.pdf              CSC-M1 TF enrichment
  FigS5_Scatter_DNMT_vs_CSC_M1.pdf    DNMT x CSC-M1 scatter
  FigS6_BP_Down_TopScatters.pdf        Top6 beta-Expr scatter
  FigS7_DNMT1_vs_BP_M1_Adjacent.pdf   DNMT1 x BP-M1 Adjacent
  FigS8_Heatmap_EpiEnzymes.pdf        epi-enzyme Z-score heatmap
  FigS9_Boxplot_EpiActivity.pdf        total epi-activity boxplot

=== Tables ===
  Hub_CSC_M1_BP_M1.csv
  Table_BP_vs_Hep_DEG.csv
  Table_C3TFT_AllModules.csv
  Table_BP_Down_BetaExprCor.csv
  Table_DualPathway_Compare.csv
  Table_DNMT_Hub_Cor.csv
  Table_Fisher_DNMT_Overlap.csv
  Table_Module_DNMT_Cor.csv
  Table_BP_Hub_BetaExprCor.csv
")

# ██  : EpiSCORE   ██                          #
#                                                                             #
###############################################################################

message("\n\n", strrep("█", 60))
message(strrep("█", 60))

# EpiSCORE
if (!exists("frac_tumor_all")) {
  f <- file.path(EPISCORE_DIR, "Table_Fraction_Tumor.csv")
  if (file.exists(f)) { frac_tumor_all <- read.csv(f); message("  restored frac_tumor_all") }
}
if (!exists("frac_normal_all")) {
  f <- file.path(EPISCORE_DIR, "Table_Fraction_Normal.csv")
  if (file.exists(f)) { frac_normal_all <- read.csv(f); message("  restored frac_normal_all") }
}

fig5_plots <- list()
fig5_stats <- list()

# Fig5A: CSC fraction vs CSC-M1 (Tumor only)
message("\n--- Fig5A: CSC EpiSCORE (Tumor) ---")

if (exists("frac_tumor_all")) {
  if (!"patient_id" %in% colnames(frac_tumor_all))
    frac_tumor_all$patient_id <- substr(frac_tumor_all$barcode, 1, 12)

  sig_tumor <- sig_scores %>% filter(tissue == "Tumor") %>% distinct(patient_id, .keep_all = TRUE)

  ct_cols <- setdiff(colnames(frac_tumor_all), c("barcode","patient_id","tissue"))
  frac_t  <- frac_tumor_all %>% distinct(patient_id, .keep_all = TRUE)
  val_t   <- inner_join(frac_t, sig_tumor, by = "patient_id")

  for (ct in ct_cols) for (mod in names(FOCUS_MODULES)) {
    ml <- FOCUS_MODULES[[mod]]$label
    if (!ct %in% colnames(val_t) || !ml %in% colnames(val_t)) next
    ok <- is.finite(val_t[[ct]]) & is.finite(val_t[[ml]]); if (sum(ok) < 10) next
    ct_test <- tryCatch(cor.test(val_t[[ct]][ok], val_t[[ml]][ok], method = "spearman"), error = function(e) NULL)
    if (is.null(ct_test)) next

    tag <- paste("Tumor", ct, mod)
    fig5_plots[[tag]] <- plot_scatter(val_t[ok, ], ct, ml,
                                      sprintf("%s fraction", ct), ml,
                                      color = "#D62728", title = sprintf("%s vs %s (Tumor)", ct, ml))
    fig5_stats[[tag]] <- data.frame(Panel = "Tumor", Fraction = ct, Module = ml,
                                     rho = ct_test$estimate, p = ct_test$p.value, n = sum(ok))
    message(sprintf("  %s vs %s: ρ=%.3f", ct, ml, ct_test$estimate))
  }
}

# Fig5B: BP fraction vs BP-M1 (Adjacent only)
message("\n--- Fig5B: Bipotent EpiSCORE (Adjacent) ---")

if (exists("frac_normal_all")) {
  if (!"patient_id" %in% colnames(frac_normal_all))
    frac_normal_all$patient_id <- substr(frac_normal_all$barcode, 1, 12)

  sig_adj  <- sig_scores %>% filter(tissue == "Adjacent") %>% distinct(patient_id, .keep_all = TRUE)
  ct_cols  <- setdiff(colnames(frac_normal_all), c("barcode","patient_id","tissue"))
  frac_n   <- frac_normal_all %>% distinct(patient_id, .keep_all = TRUE)
  val_n    <- inner_join(frac_n, sig_adj, by = "patient_id")

  for (ct in ct_cols) for (mod in names(FOCUS_MODULES)) {
    ml <- FOCUS_MODULES[[mod]]$label
    if (!ct %in% colnames(val_n) || !ml %in% colnames(val_n)) next
    ok <- is.finite(val_n[[ct]]) & is.finite(val_n[[ml]]); if (sum(ok) < 10) next
    ct_test <- tryCatch(cor.test(val_n[[ct]][ok], val_n[[ml]][ok], method = "spearman"), error = function(e) NULL)
    if (is.null(ct_test)) next

    tag <- paste("Adj", ct, mod)
    fig5_plots[[tag]] <- plot_scatter(val_n[ok, ], ct, ml,
                                      sprintf("%s fraction", ct), ml,
                                      color = "#1F77B4", title = sprintf("%s vs %s (Adjacent)", ct, ml))
    fig5_stats[[tag]] <- data.frame(Panel = "Adjacent", Fraction = ct, Module = ml,
                                     rho = ct_test$estimate, p = ct_test$p.value, n = sum(ok))
    message(sprintf("  %s vs %s: ρ=%.3f", ct, ml, ct_test$estimate))
  }
}

if (length(fig5_plots) > 0) {
  nc <- min(3, length(fig5_plots))
  pdf(file.path(FIG_DIR, "Fig5_EpiSCORE_Validation.pdf"),
      width = 5 * nc, height = 5 * ceiling(length(fig5_plots) / nc))
  do.call(grid.arrange, c(fig5_plots, ncol = nc)); dev.off()
  write.csv(bind_rows(fig5_stats), file.path(FIG_DIR, "Table_Fig5_EpiSCORE.csv"), row.names = FALSE)
}

###############################################################################

###############################################################################
#                         COMPLETE SUMMARY                                    #
###############################################################################
message("\n", strrep("=", 70))
message("  ALL FIGURES COMPLETE!")
message(strrep("=", 70))
message(sprintf("\nOutput: %s\n", FIG_DIR))
message("
=== Main Figures ===
Fig1A  ME Cell-Type Correlation Heatmap
Fig1B  Module Hallmark Enrichment
Fig2A  KM: BP-M1 in Adjacent
Fig2B  Scatter: BP-M1 vs Hoshida
Fig3A  Volcano: CSC vs Other Tumor
Fig3B  Scatter: DNMT score vs CSC hub score
Fig3C  4-panel regulation chain
Fig3D  KM x3 (DNMT, CSC, DNMTxCSC)
Fig3E  Cox Forest Plot
Fig4A  DotPlot: 14 epigenetic enzymes
Fig4B  BP-M1 TF Enrichment (DNMT1 highlighted)
Fig4C  BP-down genes beta-Expr histogram
Fig4D  Dual Pathway dot plot
Fig5A  EpiSCORE CSC fraction (Tumor)
Fig5B  EpiSCORE BP fraction (Adjacent)

=== Supplementary ===
FigS1  DNMT x CSC-M1 Hub heatmap
FigS2  Fisher DNMT enrichment
FigS3  4-module x DNMT dot plot
FigS4  CSC-M1 C3 TFT enrichment
FigS5  DNMT x CSC-M1 scatter trio
FigS6  BP-down Top6 beta-Expr scatter
FigS7  DNMT1 x BP-M1 Adjacent scatter
FigS8  Epigenetic enzymes Z-score heatmap
FigS9  Total epigenetic activity boxplot
")
