#!/usr/bin/env Rscript
# Figure 2 — Somatic mutations inferred from scRNA-seq

# --- 2a: Single-cell OncoPrint (dN/dS significant genes) ---

# !/usr/bin/env Rscript
# /home/download/csc_article/fig3/main/panel_f_oncoprint_percell.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

MAF_FILES <- c(
  "/data2/Final_24_SingleCell.maf",
  "/data2/Final_GSE156625_2.maf",
  "/data2/Final_GSE156625_SingleCell.maf",
  "/data2/Final_GSE282701_SingleCell.maf",
  "/data2/Final_SRP_SingleCell.maf"
)

META_FILE      <- "/home/download/Monopogen/monopogen_snv/cell_annotation.csv"
SIG_GENES_FILE <- "/data2/HCC_SNV_Landscape/dNdS_significant_genes.csv"

OUT_DIR <- "/home/download/csc_article/fig3/main"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

maf_list <- lapply(MAF_FILES, function(f) {
  if (file.exists(f)) fread(f, fill = TRUE) else NULL
})
maf_raw <- rbindlist(maf_list, fill = TRUE)

meta <- fread(META_FILE)
if (!"Tumor_Sample_Barcode" %in% names(meta)) setnames(meta, 1, "Tumor_Sample_Barcode")

maf_m <- merge(maf_raw, meta[, .(Tumor_Sample_Barcode, cell_type, patient_id)],
               by = "Tumor_Sample_Barcode", all.x = FALSE)

message(sprintf("  %d mutations | %d cells", nrow(maf_m), uniqueN(maf_m$Tumor_Sample_Barcode)))

if (!file.exists(SIG_GENES_FILE)) stop("SIG_GENES_FILE not found: ", SIG_GENES_FILE)

sig_genes_dt <- fread(SIG_GENES_FILE)
gene_col <- intersect(names(sig_genes_dt), c("gene_name", "Hugo_Symbol", "gene"))
if (length(gene_col) == 0) stop("Gene name column not found")
sig_genes <- unique(sig_genes_dt[[gene_col[1]]])

maf_sig <- maf_m[Hugo_Symbol %in% sig_genes]

classify_variant <- function(vc) {
  vc <- tolower(vc)
  ifelse(grepl("missense", vc), "Missense",
  ifelse(grepl("nonsense|stop_gained", vc), "Nonsense",
  ifelse(grepl("splice", vc), "Splice",
  ifelse(grepl("synonymous|silent", vc), "Synonymous",
  ifelse(grepl("frame_shift|frameshift", vc), "Frameshift",
  ifelse(grepl("in_frame", vc), "InFrame", "Other"))))))
}
maf_sig[, MutType := classify_variant(Variant_Classification)]

print(maf_sig[, .N, by = MutType][order(-N)])

cells_with_mut <- unique(maf_sig$Tumor_Sample_Barcode)
valid_cells <- unique(maf_m$Tumor_Sample_Barcode)
cells_no_mut <- setdiff(valid_cells, cells_with_mut)

set.seed(42)
MAX_BG <- 500; MAX_MUT <- 2000
cells_no_mut_sample <- if (length(cells_no_mut) > MAX_BG) sample(cells_no_mut, MAX_BG) else cells_no_mut
cells_with_mut_plot <- if (length(cells_with_mut) > MAX_MUT) sample(cells_with_mut, MAX_MUT) else cells_with_mut
all_plot_cells <- c(cells_with_mut_plot, cells_no_mut_sample)

                length(all_plot_cells), length(cells_with_mut_plot), length(cells_no_mut_sample)))

mut_types <- unique(maf_sig$MutType)
mat_list <- list()

for (mt in mut_types) {
  sub <- maf_sig[MutType == mt & Tumor_Sample_Barcode %in% all_plot_cells]
  m <- matrix(FALSE, nrow = length(sig_genes), ncol = length(all_plot_cells),
              dimnames = list(sig_genes, all_plot_cells))
  if (nrow(sub) > 0) {
    for (i in seq_len(nrow(sub))) {
      g <- sub$Hugo_Symbol[i]; bc <- sub$Tumor_Sample_Barcode[i]
      if (g %in% sig_genes && bc %in% all_plot_cells) m[g, bc] <- TRUE
    }
  }
  mat_list[[mt]] <- m
}

# ========================  &  ========================
anno_dt <- meta[Tumor_Sample_Barcode %in% all_plot_cells]
anno_dt <- anno_dt[match(all_plot_cells, Tumor_Sample_Barcode)]

has_patient <- "patient_id" %in% names(anno_dt)
if (has_patient) {
  order_idx <- order(anno_dt$cell_type, anno_dt$patient_id)
} else {
  order_idx <- order(anno_dt$cell_type)
}
all_plot_cells_sorted <- all_plot_cells[order_idx]
for (mt in names(mat_list)) mat_list[[mt]] <- mat_list[[mt]][, all_plot_cells_sorted, drop = FALSE]
anno_dt <- anno_dt[order_idx]

ct_levels <- sort(unique(anno_dt$cell_type))
ct_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
                "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7","#dbdb8d","#9edae5")
ct_colors <- setNames(ct_palette[seq_along(ct_levels)], ct_levels)

top_anno_args <- list(
  Cell_Type = anno_dt$cell_type,
  col = list(Cell_Type = ct_colors),
  show_legend = TRUE, annotation_name_side = "left", show_annotation_name = TRUE
)
if (has_patient) {
  pat_levels <- sort(unique(anno_dt$patient_id))
  pat_palette <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                   "#FFFF33","#A65628","#F781BF","#999999","#66C2A5",
                   "#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F",
                   "#E5C494","#B3B3B3","#80B1D3","#FDB462","#FCCDE5")
  pat_colors <- setNames(pat_palette[seq_along(pat_levels)], pat_levels)
  top_anno_args$Patient <- anno_dt$patient_id
  top_anno_args$col$Patient <- pat_colors
}
ha_top <- do.call(HeatmapAnnotation, top_anno_args)

# ========================  & alter_fun ========================
mut_colors <- c(
  "Missense" = "#26A269", "Nonsense" = "#E74C3C", "Splice" = "#F39C12",
  "Synonymous" = "#3498DB", "Frameshift" = "#8E44AD", "InFrame" = "#E91E63", "Other" = "#95A5A6"
)
mut_colors <- mut_colors[names(mut_colors) %in% mut_types]

alter_fun_list <- list(
  background = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#F5F5F5", col = NA))
)
for (mt in names(mut_colors)) {
  local({
    col <- mut_colors[[mt]]
    alter_fun_list[[mt]] <<- function(x, y, w, h) {
      grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = col, col = NA))
    }
  })
}

n_cells <- length(all_plot_cells_sorted)
fig_w <- min(60, max(16, n_cells * 0.012 + 6))
fig_h <- min(20, max(8, length(sig_genes) * 0.5 + 4))

out_file <- file.path(OUT_DIR, "panel_f_oncoprint_percell.pdf")

pdf(out_file, width = fig_w, height = fig_h)

ht <- oncoPrint(
  mat_list,
  alter_fun = alter_fun_list,
  alter_fun_is_vectorized = FALSE,
  col = mut_colors,
  top_annotation = ha_top,
  show_column_names = FALSE,
  show_pct = FALSE,
  remove_empty_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = sprintf("Single-Cell OncoPrint: %d Genes × %d Cells",
                         length(sig_genes), n_cells),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Mutation Type",
    at = names(mut_colors),
    labels = names(mut_colors),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

draw(ht)
dev.off()

message(sprintf("   %s", out_file))

# --- 2a supplement: dN/dS analysis + lollipop plots ---

# !/usr/bin/env Rscript
# 03_dNdS_lollipop.R —  dN/dS  + Lollipop

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(maftools)
})

MAF_FILES <- c(
  "/data2/Final_24_SingleCell.maf",
  "/data2/Final_GSE156625_2.maf",
  "/data2/Final_GSE156625_SingleCell.maf",
  "/data2/Final_GSE282701_SingleCell.maf",
  "/data2/Final_SRP_SingleCell.maf",
  "/data2/Final_MACS_SingleCell.maf"
)

META_FILE      <- "/home/download/Monopogen/monopogen_snv/cell_annotation.csv"
DNDS_SIG_FILE  <- "/data2/HCC_SNV_Landscape/dNdS_significant_genes.csv"
DNDS_GLOB_FILE <- "/data2/HCC_SNV_Landscape/dNdS_global_results.csv"

OUT_DIR <- "/home/download/csc_article/fig3/dnds_lollipop"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "lollipop_all"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "lollipop_by_group"), showWarnings = FALSE)

maf_list <- lapply(MAF_FILES, function(f) {
  if (!file.exists(f)) return(NULL)
  fread(f, fill = TRUE)
})
maf_raw <- rbindlist(maf_list[!sapply(maf_list, is.null)], fill = TRUE)

meta <- fread(META_FILE)
if (!"Tumor_Sample_Barcode" %in% names(meta)) setnames(meta, 1, "Tumor_Sample_Barcode")
meta[cell_type == "CSC", cell_type := "CSC(classic)"]

maf_m <- merge(maf_raw, meta[, .(Tumor_Sample_Barcode, cell_type, patient_id)],
               by = "Tumor_Sample_Barcode", all.x = FALSE)
message(sprintf("  %d mutations | %d cells", nrow(maf_m), uniqueN(maf_m$Tumor_Sample_Barcode)))

# dN/dS
sig_dt <- fread(DNDS_SIG_FILE)
sig_genes <- unique(sig_dt$gene_name)

# ======================== Patient-level MAF for maftools ========================
maf_patient <- copy(maf_m)
maf_patient[, Tumor_Sample_Barcode := patient_id]
maf_patient <- unique(maf_patient, by = c("Tumor_Sample_Barcode", "Hugo_Symbol",
                                           "Chromosome", "Start_Position"))
tmp_maf <- file.path(OUT_DIR, ".tmp_all.maf")
fwrite(maf_patient, tmp_maf, sep = "\t")
maf_obj <- read.maf(tmp_maf)

for (g in sig_genes) {
  if (!g %in% maf_obj@data$Hugo_Symbol) {
    next
  }
  tryCatch({
    pdf(file.path(OUT_DIR, "lollipop_all", paste0(g, ".pdf")), width = 10, height = 4)
    lollipopPlot(maf = maf_obj, gene = g, showMutationRate = FALSE,
                 labelPos = "all", labPosSize = 0.8,
                 printCount = TRUE, cBioPortal = FALSE)
    dev.off()
    message(sprintf("     %s", g))
  }, error = function(e) {
    dev.off()
    message(sprintf("     %s: %s", g, e$message))
  })
}

group_def <- list(
  Normal = c("Bipotent Progenitor", "Hepatocytes", "Cholangiocytes"),
  Tumor  = c("Hepatocytes (cycling)", "Hepatocytes (Metabolic-H)",
             "Hepatocytes (Metabolic-L)", "Hepatocytes (Stem)",
             "Hepatocytes (Quiescent)", "Biliary-like"),
  CSC    = c("CSC(classic)", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)")
)

maf_objs <- list()
for (grp in names(group_def)) {
  sub <- maf_m[cell_type %in% group_def[[grp]]]
  if (nrow(sub) < 50) next
  sub[, Tumor_Sample_Barcode := patient_id]
  sub <- unique(sub, by = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position"))
  tmp_f <- file.path(OUT_DIR, paste0(".tmp_", grp, ".maf"))
  fwrite(sub, tmp_f, sep = "\t")
  maf_objs[[grp]] <- tryCatch(read.maf(tmp_f), error = function(e) NULL)
}

for (g in sig_genes) {
  for (grp in names(maf_objs)) {
    obj <- maf_objs[[grp]]
    if (is.null(obj) || !g %in% obj@data$Hugo_Symbol) next
    tryCatch({
      pdf(file.path(OUT_DIR, "lollipop_by_group", paste0(g, "_", grp, ".pdf")),
          width = 10, height = 4)
      lollipopPlot(maf = obj, gene = g, showMutationRate = FALSE,
                   labelPos = "all", labPosSize = 0.8)
      title(sub = paste0("[", grp, "] ", paste(group_def[[grp]], collapse = " + ")),
            cex.sub = 0.65, col.sub = "grey40")
      dev.off()
    }, error = function(e) { dev.off() })
  }
  
  if (!is.null(maf_objs[["Normal"]]) && !is.null(maf_objs[["CSC"]])) {
    if (g %in% maf_objs[["Normal"]]@data$Hugo_Symbol &&
        g %in% maf_objs[["CSC"]]@data$Hugo_Symbol) {
      tryCatch({
        pdf(file.path(OUT_DIR, "lollipop_by_group", paste0(g, "_Normal_vs_CSC.pdf")),
            width = 10, height = 5)
        lollipopPlot2(m1 = maf_objs[["Normal"]], m2 = maf_objs[["CSC"]],
                      gene = g, m1_name = "Normal", m2_name = "CSC",
                      showDomainLabel = TRUE)
        dev.off()
      }, error = function(e) { dev.off() })
    }
  }
  
  if (!is.null(maf_objs[["Tumor"]]) && !is.null(maf_objs[["CSC"]])) {
    if (g %in% maf_objs[["Tumor"]]@data$Hugo_Symbol &&
        g %in% maf_objs[["CSC"]]@data$Hugo_Symbol) {
      tryCatch({
        pdf(file.path(OUT_DIR, "lollipop_by_group", paste0(g, "_Tumor_vs_CSC.pdf")),
            width = 10, height = 5)
        lollipopPlot2(m1 = maf_objs[["Tumor"]], m2 = maf_objs[["CSC"]],
                      gene = g, m1_name = "Tumor(non-CSC)", m2_name = "CSC",
                      showDomainLabel = TRUE)
        dev.off()
      }, error = function(e) { dev.off() })
    }
  }
}

sig_dt[, sig_type := fcase(
  qmis_cv < 0.1 & qtrunc_cv < 0.1, "Both",
  qmis_cv < 0.1, "Missense",
  qtrunc_cv < 0.1, "Truncating",
  qallsubs_cv < 0.1, "All subs",
  default = "Other"
)]

# -log10(q)
sig_dt[, neg_log_q := -log10(qallsubs_cv + 1e-20)]

ct_order <- sig_dt[, .(n = .N), by = cell_type][order(-n)]$cell_type
gene_order <- sig_dt[, .(min_q = min(qallsubs_cv)), by = gene_name][order(min_q)]$gene_name

sig_dt[, cell_type := factor(cell_type, levels = ct_order)]
sig_dt[, gene_name := factor(gene_name, levels = rev(gene_order))]

pSummary <- ggplot(sig_dt, aes(x = cell_type, y = gene_name)) +
  geom_point(aes(size = neg_log_q, color = sig_type)) +
  scale_size_continuous(range = c(2, 8), name = expression(-log[10](q))) +
  scale_color_manual(values = c("Missense" = "#0369A1", "Truncating" = "#DC2626",
                                 "Both" = "#7B2D8E", "All subs" = "#F59E0B",
                                 "Other" = "#9CA3AF"),
                     name = "Selection type") +
  labs(title = "dN/dS significant genes across cell types (q < 0.1)",
       subtitle = "Bubble size = significance strength",
       x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
        axis.text.y = element_text(face = "italic", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.major = element_line(color = "#F1F5F9"),
        legend.position = "right")

ggsave(file.path(OUT_DIR, "dNdS_significant_genes_bubble.pdf"), pSummary, width = 10, height = 6)
ggsave(file.path(OUT_DIR, "dNdS_significant_genes_bubble.png"), pSummary, width = 10, height = 6, dpi = 300)

maf_m[, region := fcase(
  Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Silent",
                                 "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del",
                                 "Frame_Shift_Ins", "Nonstop_Mutation",
                                 "Translation_Start_Site"), "Exonic",
  Variant_Classification %in% c("Splice_Site", "Splice_Region"), "Splice",
  Variant_Classification %in% c("3'UTR", "5'UTR", "3'Flank", "5'Flank"), "UTR/Flank",
  Variant_Classification == "Intron", "Intronic",
  Variant_Classification == "IGR", "Intergenic",
  default = "Other"
)]

region_stats <- maf_m[, .N, by = region][order(-N)]
region_stats[, pct := round(N / sum(N) * 100, 1)]
print(region_stats)
fwrite(region_stats, file.path(OUT_DIR, "genomic_region_counts.csv"))

file.remove(list.files(OUT_DIR, pattern = "^\\.tmp", full.names = TRUE))

# 2d,e,f: SNV mutation burden (shared-site O/E, LMM + paired Wilcoxon)

# !/usr/bin/env Rscript
# 02f_shared_site_TMB_with_figures.R  (v7: LMM + Wilcoxon supplement)
#
#   1. CSC → CSC(classic)
# 3. : lmer(log(OE) ~ cell_type + (1|patient_id))  —
# Wilcoxon signed-rank (patient-level median) —

library(data.table)
library(MASS)
library(lmerTest)   # lmer + Satterthwaite P-value
library(ggplot2)
library(gridExtra)
library(scales)

base_dir       <- "/data2/nb_correction_results/shared_site_v3"
evaluable_file <- file.path(base_dir, "per_sample_evaluable_sites.rds")
mapping_file   <- file.path(base_dir, "cell_to_mono_sample_mapping_unified.csv")
out_dir        <- file.path(base_dir, "noSplice")
fig_dir        <- "/home/download/csc_article/fig3/LMM"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(fig_dir, "paired_shared_all"), showWarnings = FALSE)
dir.create(file.path(fig_dir, "paired_shared_merged"), showWarnings = FALSE)

maf_files <- c(
  "/data2/Final_24_SingleCell.maf",
  "/data2/Final_GSE156625_2.maf",
  "/data2/Final_GSE156625_SingleCell.maf",
  "/data2/Final_GSE282701_SingleCell.maf",
  "/data2/Final_SRP_SingleCell.maf",
  "/data2/Final_MACS_SingleCell.maf"
)

ct_colors <- c(
  "CSC(classic)"               = "#E41A1C",
  "CSC(progenitor_cycling)"    = "#FF6B6B",
  "CSC(progenitor_quiet)"      = "#FF9999",
  "CSC_merged"                 = "#C62828",
  "Bipotent Progenitor"        = "#FF7F00",
  "Hepatocytes (cycling)"      = "#984EA3",
  "Hepatocytes (Metabolic-L)"  = "#4DAF4A",
  "Hepatocytes (Metabolic-H)"  = "#377EB8",
  "Hepatocytes (Stem)"         = "#A65628",
  "Hepatocytes"                = "#999999",
  "Hepatocytes (Quiescent)"    = "#66C2A5",
  "Cholangiocytes"             = "#E7298A",
  "Biliary-like"               = "#FFED6F",
  "Progenitor-like"            = "#B3DE69"
)

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12, face = "bold"),
    plot.title  = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 15)
  )

message(paste(rep("=", 70), collapse = ""))
message(paste(rep("=", 70), collapse = ""))

cell_map <- fread(mapping_file)
cell_map[cell_type == "CSC", cell_type := "CSC(classic)"]
                nrow(cell_map), uniqueN(cell_map$cell_type)))
print(cell_map[, .N, by = cell_type][order(-N)])

# --- 2.  per-sample evaluable sites ---
if (!file.exists(evaluable_file)) stop("Run 02d_step1_sample_site_coverage.R first")
sample_evaluable <- readRDS(evaluable_file)

# --- 3.  MAF ---
maf_list <- lapply(maf_files, function(f) {
  if (!file.exists(f)) { message(sprintf("   %s not found", basename(f))); return(NULL) }
  fread(f, skip = "Hugo_Symbol",
        select = c("Chromosome", "Start_Position", "Tumor_Sample_Barcode", "Variant_Classification"))
})
maf_all <- rbindlist(maf_list[!sapply(maf_list, is.null)], fill = TRUE)

# Splice_Site
n_before <- nrow(maf_all)
maf_all <- maf_all[!Variant_Classification %in% c("Splice_Site", "Splice_Region")]
                n_before, nrow(maf_all), n_before - nrow(maf_all),
                (n_before - nrow(maf_all)) / n_before * 100))
maf_all[, Variant_Classification := NULL]
maf_all[, Chromosome := as.character(Chromosome)]
maf_all[, Chromosome := ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome))]
maf_all[, site_id := paste0(Chromosome, ":", Start_Position)]
                nrow(maf_all), uniqueN(maf_all$site_id),
                uniqueN(maf_all$Tumor_Sample_Barcode)))

run_pairwise <- function(anno, maf_all, sample_evaluable, label = "") {

  cell_types <- sort(unique(anno$cell_type))

  ct_evaluable <- list()
  for (ct in cell_types) {
    samples_ct <- unique(anno[cell_type == ct, mono_sample])
    samples_ct <- intersect(samples_ct, names(sample_evaluable))
    if (length(samples_ct) == 0) {
      ct_evaluable[[ct]] <- character(0); next
    }
    ct_evaluable[[ct]] <- unique(unlist(sample_evaluable[samples_ct]))
    message(sprintf("    %s: %d samples, %d evaluable sites",
                    ct, length(samples_ct), length(ct_evaluable[[ct]])))
  }

  cell_site_mut <- unique(maf_all[Tumor_Sample_Barcode %in% anno$Tumor_Sample_Barcode,
                                   .(Tumor_Sample_Barcode, site_id)])

  pairs <- combn(cell_types, 2, simplify = FALSE)
  results <- list()
  per_cell_oe <- list()

  for (i in seq_along(pairs)) {
    t1 <- pairs[[i]][1]; t2 <- pairs[[i]][2]

    shared_sites <- intersect(ct_evaluable[[t1]], ct_evaluable[[t2]])
    n_shared <- length(shared_sites)

    empty_row <- data.table(
      CellType1 = t1, CellType2 = t2, N_shared_sites = n_shared,
      N_patients = 0L, N_cells_T1 = 0L, N_cells_T2 = 0L,
      Med_OE_1 = NA_real_, Med_OE_2 = NA_real_,
      Med_Eval_1 = NA_real_, Med_Eval_2 = NA_real_,
      N_hotspot = NA_integer_, Med_BMR = NA_real_,
      # LMM
      LMM_beta = NA_real_, LMM_CI_lo = NA_real_, LMM_CI_hi = NA_real_,
      LMM_P = NA_real_, LMM_converged = NA,
      # Wilcoxon (supplementary)
      W_stat = NA_real_, Wilcox_P = NA_real_,
      Direction = NA_character_)

    if (n_shared < 10) { results[[i]] <- empty_row; next }

    message(sprintf("  [%d] %s vs %s: %d shared evaluable sites...",
                    i, t1, t2, n_shared))

    cells_pair <- anno[cell_type %in% c(t1, t2)]
    unique_samples <- unique(cells_pair$mono_sample)

    sample_shared_sites <- list()
    sample_n_covered <- integer(0)
    for (s in unique_samples) {
      ev <- sample_evaluable[[s]]
      if (is.null(ev)) {
        sample_shared_sites[[s]] <- character(0)
        sample_n_covered[s] <- 0L
      } else {
        cov_s <- intersect(ev, shared_sites)
        sample_shared_sites[[s]] <- cov_s
        sample_n_covered[s] <- length(cov_s)
      }
    }
    cells_pair[, N_covered := sample_n_covered[mono_sample]]
    cells_pair <- cells_pair[N_covered > 0]

    n1 <- cells_pair[cell_type == t1, .N]
    n2 <- cells_pair[cell_type == t2, .N]
    if (n1 < 5 || n2 < 5) {
      empty_row$N_cells_T1 <- n1; empty_row$N_cells_T2 <- n2
      results[[i]] <- empty_row; next
    }

    # ---  pooled BMR ---
    sample_ncells <- cells_pair[, .N, by = mono_sample]
    setkey(sample_ncells, mono_sample)

    site_cov_parts <- list()
    for (s in unique(cells_pair$mono_sample)) {
      ss <- sample_shared_sites[[s]]
      if (length(ss) == 0) next
      nc <- sample_ncells[s, N]
      site_cov_parts[[s]] <- data.table(site_id = ss, n_cells = nc)
    }
    site_cov <- rbindlist(site_cov_parts)
    site_cov <- site_cov[, .(Total_Covered = sum(n_cells)), by = site_id]

    site_mut <- cell_site_mut[site_id %in% shared_sites &
                                Tumor_Sample_Barcode %in% cells_pair$Tumor_Sample_Barcode]
    site_mut_count <- site_mut[, .(N_Mutated = .N), by = site_id]

    site_bmr <- merge(site_cov, site_mut_count, by = "site_id", all.x = TRUE)
    site_bmr[is.na(N_Mutated), N_Mutated := 0L]
    site_bmr[, mu_j := N_Mutated / Total_Covered]
    setkey(site_bmr, site_id)

    n_hotspot <- site_bmr[mu_j > 0.01, .N]
    n_with_mut <- site_bmr[N_Mutated > 0, .N]
    med_bmr <- fifelse(n_with_mut > 0, median(site_bmr[N_Mutated > 0, mu_j]), 0)

    message(sprintf("    BMR: %d/%d sites mutated, %d hotspots (>1%%), medBMR=%.5f",
                    n_with_mut, n_shared, n_hotspot, med_bmr))

    # --- Expected = Σ μ_j ---
    sample_expected <- numeric(length(unique(cells_pair$mono_sample)))
    names(sample_expected) <- unique(cells_pair$mono_sample)
    for (s in names(sample_expected)) {
      ss <- sample_shared_sites[[s]]
      if (length(ss) == 0) { sample_expected[s] <- 0; next }
      mu_vals <- site_bmr[ss, mu_j, on = "site_id"]
      mu_vals <- mu_vals[!is.na(mu_vals)]
      sample_expected[s] <- sum(mu_vals)
    }
    cells_pair[, Expected := sample_expected[mono_sample]]

    # --- Observed ---
    obs_per_cell <- cell_site_mut[site_id %in% shared_sites &
                                    Tumor_Sample_Barcode %in% cells_pair$Tumor_Sample_Barcode,
                                   .(Obs = .N), by = Tumor_Sample_Barcode]
    cells_pair <- merge(cells_pair, obs_per_cell, by = "Tumor_Sample_Barcode", all.x = TRUE)
    cells_pair[is.na(Obs), Obs := 0L]

    # --- O/E ---
    cells_pair[, OE := fifelse(Expected > 0, Obs / Expected, NA_real_)]
    cells_pair <- cells_pair[!is.na(OE)]

    n1 <- cells_pair[cell_type == t1, .N]
    n2 <- cells_pair[cell_type == t2, .N]
    if (n1 < 5 || n2 < 5) {
      empty_row$N_cells_T1 <- n1; empty_row$N_cells_T2 <- n2
      empty_row$N_hotspot <- n_hotspot; empty_row$Med_BMR <- med_bmr
      results[[i]] <- empty_row; next
    }

    patients_both <- intersect(
      cells_pair[cell_type == t1, unique(patient_id)],
      cells_pair[cell_type == t2, unique(patient_id)])
    n_patients <- length(patients_both)

    if (n_patients < 3) {
      empty_row$N_patients <- n_patients
      empty_row$N_cells_T1 <- n1; empty_row$N_cells_T2 <- n2
      empty_row$N_hotspot <- n_hotspot; empty_row$Med_BMR <- med_bmr
      results[[i]] <- empty_row; next
    }

    cells_lmm <- cells_pair[patient_id %in% patients_both]
    n1_lmm <- cells_lmm[cell_type == t1, .N]
    n2_lmm <- cells_lmm[cell_type == t2, .N]

    per_cell_oe[[i]] <- cells_lmm[, .(Tumor_Sample_Barcode, cell_type, patient_id,
                                       OE, N_covered, Obs, Expected,
                                       comparison = paste0(t1, " vs ", t2))]

    # === log(O/E)  ===
    # O/E=0 :  pseudocount = min(nonzero OE) / 2
    nonzero_oe <- cells_lmm[OE > 0, OE]
    if (length(nonzero_oe) > 0) {
      pseudo <- min(nonzero_oe) / 2
    } else {
      pseudo <- 1e-4
    }
    cells_lmm[, log_OE := log(OE + pseudo)]

    cells_lmm[, ct_factor := factor(cell_type, levels = c(t1, t2))]

    # ===================== LMM  =====================
    lmm_beta <- NA_real_; lmm_ci_lo <- NA_real_; lmm_ci_hi <- NA_real_
    lmm_p <- NA_real_; lmm_ok <- FALSE

    tryCatch({
      fit <- lmer(log_OE ~ ct_factor + (1 | patient_id), data = cells_lmm,
                  REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 20000)))

      sfit <- summary(fit)
      coef_row <- sfit$coefficients
      # ct_factor  (t2 vs t1)
      if (nrow(coef_row) >= 2) {
        lmm_beta <- coef_row[2, "Estimate"]
        lmm_se   <- coef_row[2, "Std. Error"]
        lmm_p    <- coef_row[2, "Pr(>|t|)"]
        lmm_ci_lo <- lmm_beta - 1.96 * lmm_se
        lmm_ci_hi <- lmm_beta + 1.96 * lmm_se
        lmm_ok <- TRUE
      }
    }, error = function(e) {
      message(sprintf("     LMM failed: %s", conditionMessage(e)))
    }, warning = function(w) {
      message(sprintf("     LMM warning: %s", conditionMessage(w)))
      tryCatch({
        fit <- suppressWarnings(
          lmer(log_OE ~ ct_factor + (1 | patient_id), data = cells_lmm,
               REML = TRUE,
               control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 20000))))
        sfit <- summary(fit)
        coef_row <- sfit$coefficients
        if (nrow(coef_row) >= 2) {
          lmm_beta  <<- coef_row[2, "Estimate"]
          lmm_se    <- coef_row[2, "Std. Error"]
          lmm_p     <<- coef_row[2, "Pr(>|t|)"]
          lmm_ci_lo <<- lmm_beta - 1.96 * lmm_se
          lmm_ci_hi <<- lmm_beta + 1.96 * lmm_se
          lmm_ok    <<- TRUE
        }
      }, error = function(e2) {
        message(sprintf("     LMM retry failed: %s", conditionMessage(e2)))
      })
    })

    # ===================== Wilcoxon  =====================
    pat <- cells_lmm[, .(Med_OE = median(OE, na.rm = TRUE)),
                      by = .(patient_id, cell_type)]
    d1 <- pat[cell_type == t1, .(patient_id, OE_1 = Med_OE)]
    d2 <- pat[cell_type == t2, .(patient_id, OE_2 = Med_OE)]
    paired <- merge(d1, d2, by = "patient_id")

    wt_2s <- wilcox.test(paired$OE_1, paired$OE_2, paired = TRUE, exact = FALSE)

    # patient-level median
    direction <- ifelse(median(paired$OE_1) > median(paired$OE_2),
                        paste0(t1, " > ", t2), paste0(t2, " > ", t1))
    # LMM : beta < 0 → t1 > t2; beta > 0 → t2 > t1
    if (!is.na(lmm_beta)) {
      direction_lmm <- ifelse(lmm_beta < 0,
                               paste0(t1, " > ", t2), paste0(t2, " > ", t1))
      if (direction != direction_lmm)
        message(sprintf("     Direction mismatch: Wilcox=%s, LMM=%s", direction, direction_lmm))
    }

    med_eval_1 <- median(cells_lmm[cell_type == t1, N_covered])
    med_eval_2 <- median(cells_lmm[cell_type == t2, N_covered])

    message(sprintf("    -> n=%d patients (%d+%d cells), LMM_P=%s, Wilcox_P=%.4f, %s",
                    n_patients, n1_lmm, n2_lmm,
                    ifelse(is.na(lmm_p), "NA", sprintf("%.4f", lmm_p)),
                    wt_2s$p.value, direction))

    results[[i]] <- data.table(
      CellType1 = t1, CellType2 = t2, N_shared_sites = n_shared,
      N_patients = n_patients, N_cells_T1 = n1_lmm, N_cells_T2 = n2_lmm,
      Med_OE_1 = round(median(paired$OE_1), 4),
      Med_OE_2 = round(median(paired$OE_2), 4),
      Med_Eval_1 = med_eval_1, Med_Eval_2 = med_eval_2,
      N_hotspot = n_hotspot, Med_BMR = med_bmr,
      LMM_beta = round(lmm_beta, 4), LMM_CI_lo = round(lmm_ci_lo, 4),
      LMM_CI_hi = round(lmm_ci_hi, 4), LMM_P = lmm_p, LMM_converged = lmm_ok,
      W_stat = wt_2s$statistic, Wilcox_P = wt_2s$p.value,
      Direction = direction)
  }

  pw_dt <- rbindlist(results, fill = TRUE)
  pw_dt[!is.na(LMM_P), LMM_P_adj := p.adjust(LMM_P, method = "BH")]
  pw_dt[!is.na(Wilcox_P), Wilcox_P_adj := p.adjust(Wilcox_P, method = "BH")]

  list(pw = pw_dt, per_cell_oe = rbindlist(per_cell_oe, fill = TRUE))
}

anno_orig <- unique(cell_map[, .(Tumor_Sample_Barcode, cell_type, patient_id, mono_sample)])
res_orig <- run_pairwise(anno_orig, maf_all, sample_evaluable, "original")

print(res_orig$pw[LMM_P < 0.1][order(LMM_P)])
fwrite(res_orig$pw, file.path(out_dir, "shared_site_pairwise_all_subtypes.csv"))

# 5. CSC_merged
anno_m <- copy(anno_orig)
anno_m[cell_type %in% c("CSC(classic)", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)"),
       cell_type := "CSC_merged"]
res_merged <- run_pairwise(anno_m, maf_all, sample_evaluable, "CSC_merged")

print(res_merged$pw[LMM_P < 0.1][order(LMM_P)])
fwrite(res_merged$pw, file.path(out_dir, "shared_site_pairwise_CSCmerged.csv"))

# ----  (patient-median + LMM P-value) ----
make_paired_plot <- function(pw_dt, anno_dt, maf_all, sample_evaluable, t1, t2) {

  row <- pw_dt[(CellType1 == t1 & CellType2 == t2) |
                 (CellType1 == t2 & CellType2 == t1)]
  if (nrow(row) == 0 || (is.na(row$LMM_P[1]) && is.na(row$Wilcox_P[1]))) return(NULL)

  samples_t1 <- unique(anno_dt[cell_type == t1, mono_sample])
  samples_t2 <- unique(anno_dt[cell_type == t2, mono_sample])
  eval_t1 <- unique(unlist(sample_evaluable[intersect(samples_t1, names(sample_evaluable))]))
  eval_t2 <- unique(unlist(sample_evaluable[intersect(samples_t2, names(sample_evaluable))]))
  shared <- intersect(eval_t1, eval_t2)
  if (length(shared) < 10) return(NULL)

  cells_pair <- anno_dt[cell_type %in% c(t1, t2)]
  unique_samples <- unique(cells_pair$mono_sample)

  sample_shared_sites <- list()
  for (s in unique_samples) {
    ev <- sample_evaluable[[s]]
    if (is.null(ev)) { sample_shared_sites[[s]] <- character(0); next }
    sample_shared_sites[[s]] <- intersect(ev, shared)
  }
  cells_pair[, N_covered := sapply(mono_sample, function(s) length(sample_shared_sites[[s]]))]
  cells_pair <- cells_pair[N_covered > 0]

  # pooled BMR
  cell_site_mut_local <- unique(maf_all[Tumor_Sample_Barcode %in% cells_pair$Tumor_Sample_Barcode &
                                          site_id %in% shared, .(Tumor_Sample_Barcode, site_id)])

  sample_ncells <- cells_pair[, .N, by = mono_sample]
  setkey(sample_ncells, mono_sample)

  site_cov_parts <- list()
  for (s in unique(cells_pair$mono_sample)) {
    ss <- sample_shared_sites[[s]]
    if (length(ss) == 0) next
    site_cov_parts[[s]] <- data.table(site_id = ss, n_cells = sample_ncells[s, N])
  }
  site_cov <- rbindlist(site_cov_parts)[, .(Total_Covered = sum(n_cells)), by = site_id]
  site_mut_count <- cell_site_mut_local[, .(N_Mutated = .N), by = site_id]
  site_bmr <- merge(site_cov, site_mut_count, by = "site_id", all.x = TRUE)
  site_bmr[is.na(N_Mutated), N_Mutated := 0L]
  site_bmr[, mu_j := N_Mutated / Total_Covered]
  setkey(site_bmr, site_id)

  sample_expected <- numeric(length(unique(cells_pair$mono_sample)))
  names(sample_expected) <- unique(cells_pair$mono_sample)
  for (s in names(sample_expected)) {
    ss <- sample_shared_sites[[s]]
    if (length(ss) == 0) { sample_expected[s] <- 0; next }
    mu_vals <- site_bmr[ss, mu_j, on = "site_id"]
    sample_expected[s] <- sum(mu_vals[!is.na(mu_vals)])
  }
  cells_pair[, Expected := sample_expected[mono_sample]]

  obs_pc <- cell_site_mut_local[, .(Obs = .N), by = Tumor_Sample_Barcode]
  cells_pair <- merge(cells_pair, obs_pc, by = "Tumor_Sample_Barcode", all.x = TRUE)
  cells_pair[is.na(Obs), Obs := 0L]
  cells_pair[, OE := fifelse(Expected > 0, Obs / Expected, NA_real_)]
  cells_pair <- cells_pair[!is.na(OE)]

  patients_both <- intersect(
    cells_pair[cell_type == t1, unique(patient_id)],
    cells_pair[cell_type == t2, unique(patient_id)])
  cells_pair <- cells_pair[patient_id %in% patients_both]

  pat <- cells_pair[, .(Med_OE = median(OE, na.rm = TRUE)), by = .(patient_id, cell_type)]
  d1 <- pat[cell_type == t1, .(patient_id, OE_1 = Med_OE)]
  d2 <- pat[cell_type == t2, .(patient_id, OE_2 = Med_OE)]
  paired <- merge(d1, d2, by = "patient_id")
  if (nrow(paired) < 3) return(NULL)

  lmm_p <- row$LMM_P[1]
  wlx_p <- row$Wilcox_P[1]
  n_t1 <- row$N_cells_T1[1]
  n_t2 <- row$N_cells_T2[1]

  # P  LMM,  LMM  Wilcoxon
  if (!is.na(lmm_p)) {
    p_show <- lmm_p
    method_lab <- "LMM"
  } else {
    p_show <- wlx_p
    method_lab <- "Wilcoxon"
  }

  plot_dt <- rbind(
    paired[, .(patient_id, cell_type = t1, Median_OE = OE_1)],
    paired[, .(patient_id, cell_type = t2, Median_OE = OE_2)])
  plot_dt[, cell_type := factor(cell_type, levels = c(t1, t2))]

  sig <- ifelse(p_show < 0.001, "***", ifelse(p_show < 0.01, "**",
           ifelse(p_show < 0.05, "*", "ns")))
  plab <- ifelse(p_show < 0.001, sprintf("P = %.1e (%s)", p_show, method_lab),
                 sprintf("P = %.3f (%s)", p_show, method_lab))

  ggplot(plot_dt, aes(x = cell_type, y = Median_OE)) +
    geom_boxplot(aes(fill = cell_type), alpha = 0.3, width = 0.45,
                 outlier.shape = NA, linewidth = 0.3) +
    geom_point(aes(color = cell_type), size = 2) +
    geom_line(aes(group = patient_id), color = "grey60", alpha = 0.5, linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 0.3) +
    scale_fill_manual(values = ct_colors) + scale_color_manual(values = ct_colors) +
    annotate("text", x = 1.5, y = max(plot_dt$Median_OE, na.rm = TRUE) * 1.12,
             label = paste0(sig, "\n", plab), size = 3.2) +
    labs(x = "", y = "Median O/E per patient",
         title = paste0(t1, " vs ", t2),
         subtitle = sprintf("n=%d patients (%d+%d cells), %d shared sites",
                            nrow(paired), n_t1, n_t2, length(shared))) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
          plot.title = element_text(size = 9, face = "bold"))
}

pw_orig <- res_orig$pw[(!is.na(LMM_P) | !is.na(Wilcox_P)) & N_patients >= 3]
plots_orig <- list()

for (k in seq_len(nrow(pw_orig))) {
  t1 <- pw_orig$CellType1[k]; t2 <- pw_orig$CellType2[k]
  p <- make_paired_plot(pw_orig, anno_orig, maf_all, sample_evaluable, t1, t2)
  if (!is.null(p)) {
    nm <- gsub("[^A-Za-z0-9_]", "_", paste0(t1, "_vs_", t2))
    plots_orig[[nm]] <- p
    ggsave(file.path(fig_dir, "paired_shared_all", paste0(nm, ".pdf")), p,
           width = 4.5, height = 4.5)
  }
}

if (length(plots_orig) > 0) {
  pdf(file.path(fig_dir, "SharedSite_all_paired.pdf"), width = 5, height = 5)
  for (pp in plots_orig) print(pp)
  dev.off()
}

# ---- 6b. CSC_merged  ----
pw_m <- res_merged$pw[(!is.na(LMM_P) | !is.na(Wilcox_P)) & N_patients >= 3]
plots_m <- list()

for (k in seq_len(nrow(pw_m))) {
  t1 <- pw_m$CellType1[k]; t2 <- pw_m$CellType2[k]
  p <- make_paired_plot(pw_m, anno_m, maf_all, sample_evaluable, t1, t2)
  if (!is.null(p)) {
    nm <- gsub("[^A-Za-z0-9_]", "_", paste0(t1, "_vs_", t2))
    plots_m[[nm]] <- p
    ggsave(file.path(fig_dir, "paired_shared_merged", paste0(nm, ".pdf")), p,
           width = 4.5, height = 4.5)
  }
}

if (length(plots_m) > 0) {
  pdf(file.path(fig_dir, "SharedSite_CSCmerged_paired.pdf"), width = 5, height = 5)
  for (pp in plots_m) print(pp)
  dev.off()
}

# ---- 6c. CSC  ----
csc_names <- grep("CSC_merged", names(plots_m), value = TRUE)
if (length(csc_names) >= 1) {
  nc <- min(3, length(csc_names))
  nr <- ceiling(length(csc_names) / nc)
  figCSC <- arrangeGrob(
    grobs = plots_m[csc_names], ncol = nc, nrow = nr,
    top = grid::textGrob("Per-site BMR + LMM: CSC_merged vs all",
                          gp = grid::gpar(fontsize = 13, fontface = "bold"),
                          x = 0.02, hjust = 0))
  ggsave(file.path(fig_dir, "SharedSite_CSC_focus.pdf"), figCSC,
         width = nc * 4.5, height = nr * 4.5)
}

orig_nb_file <- "/data2/nb_correction_results/paired_wilcoxon_results.csv"
if (file.exists(orig_nb_file)) {
  orig_nb <- fread(orig_nb_file)

  orig_nb[, CellType1 := gsub("^CSC$", "CSC(classic)", CellType1)]
  orig_nb[, CellType2 := gsub("^CSC$", "CSC(classic)", CellType2)]
  orig_nb[, Direction := gsub("\\bCSC\\b", "CSC(classic)", Direction)]

  compare <- merge(
    res_merged$pw[, .(CellType1, CellType2,
                       Dir_Shared = Direction, P_Shared = LMM_P)],
    orig_nb[, .(CellType1, CellType2,
                 Dir_NB = Direction, P_NB = P_value)],
    by = c("CellType1", "CellType2"), all = TRUE
  )
  compare <- compare[!is.na(Dir_Shared) & !is.na(Dir_NB) & Dir_NB != "" & Dir_Shared != ""]
  compare[, Consistent := (Dir_Shared == Dir_NB)]

  n_con <- sum(compare$Consistent, na.rm = TRUE)
  n_tot <- nrow(compare)
                  n_con, n_tot, 100 * n_con / max(n_tot, 1)))

  fwrite(compare, file.path(out_dir, "NB_vs_shared_site_direction.csv"))

  compare_sig <- compare[pmin(P_NB, P_Shared, na.rm = TRUE) < 0.15]

  if (nrow(compare_sig) > 0) {
    get_sign <- function(dir, t1) {
      ifelse(is.na(dir) | dir == "", 0,
             ifelse(startsWith(dir, t1), 1, -1))
    }

    compare_sig[, pair_label := paste0(CellType1, " vs\n", CellType2)]
    compare_sig[, Sign_NB := mapply(get_sign, Dir_NB, CellType1)]
    compare_sig[, Sign_Shared := mapply(get_sign, Dir_Shared, CellType1)]
    compare_sig[, Val_NB := Sign_NB * pmin(-log10(pmax(P_NB, 1e-10)), 5)]
    compare_sig[, Val_Shared := Sign_Shared * pmin(-log10(pmax(P_Shared, 1e-10)), 5)]

    heat_long <- rbind(
      compare_sig[, .(pair_label, Method = "NB regression", Value = Val_NB, P = P_NB)],
      compare_sig[, .(pair_label, Method = "Shared-site LMM", Value = Val_Shared, P = P_Shared)]
    )
    pair_order <- compare_sig[order(P_NB), pair_label]
    heat_long[, pair_label := factor(pair_label, levels = rev(pair_order))]

    figHeat <- ggplot(heat_long, aes(x = Method, y = pair_label, fill = Value)) +
      geom_tile(color = "white", linewidth = 0.8) +
      geom_text(aes(label = ifelse(abs(P) < 0.001, sprintf("%.0e", P),
                                    sprintf("%.3f", P))),
                size = 2.8, color = "black") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-5, 5),
                           name = expression(-log[10](P) %*% sign)) +
      labs(x = "", y = "", title = "Direction consistency",
           subtitle = "NB regression vs shared-site LMM") +
      theme_minimal(base_size = 11) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 10, face = "bold"),
            panel.grid = element_blank(),
            plot.title = element_text(size = 13, face = "bold"))

    ggsave(file.path(fig_dir, "SharedSite_direction_heatmap.pdf"), figHeat,
           width = 7, height = max(4, nrow(compare_sig) * 0.35 + 2))
  }
}

                nrow(res_orig$pw),
                sum(res_orig$pw$N_patients >= 3, na.rm = TRUE)))
                sum(res_orig$pw$LMM_P < 0.05, na.rm = TRUE),
                sum(res_orig$pw$Wilcox_P < 0.05, na.rm = TRUE)))

message("\n=== CSC_merged ===")
                nrow(res_merged$pw),
                sum(res_merged$pw$N_patients >= 3, na.rm = TRUE)))
                sum(res_merged$pw$LMM_P < 0.05, na.rm = TRUE),
                sum(res_merged$pw$Wilcox_P < 0.05, na.rm = TRUE)))

csc_rows <- res_merged$pw[grepl("CSC_merged", CellType1) | grepl("CSC_merged", CellType2)]
print(csc_rows[order(LMM_P)][,
  .(CellType1, CellType2, N_shared_sites, N_patients, N_cells_T1, N_cells_T2,
    Med_OE_1, Med_OE_2, LMM_beta, LMM_P = signif(LMM_P, 3),
    Wilcox_P = signif(Wilcox_P, 3), Direction)])

bp_rows <- res_merged$pw[grepl("Bipotent", CellType1) | grepl("Bipotent", CellType2)]
bp_rows <- bp_rows[!is.na(LMM_P) | !is.na(Wilcox_P)]
print(bp_rows[order(LMM_P)][,
  .(CellType1, CellType2, N_patients, N_cells_T1, N_cells_T2,
    Med_OE_1, Med_OE_2, LMM_beta, LMM_P = signif(LMM_P, 3),
    Wilcox_P = signif(Wilcox_P, 3), Direction)])

message("\n", paste(rep("=", 70), collapse = ""))
message("DONE!")
message(paste(rep("=", 70), collapse = ""))

# --- 2c: Mutational signature (de novo NMF extraction) ---

# !/usr/bin/env Rscript
# scRNA-seq  —— NMF De Novo
#
#
# STEP 2:  RNA  (REDIportal)
# STEP 4: NMF rank survey ()
# STEP 5: NMF  de novo
# STEP 6:  COSMIC v3.2  ()

suppressPackageStartupMessages({
  library(data.table)
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(NMF)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(RColorBrewer)
})

MAF_FILES <- c(
  "/data2/Final_24_SingleCell.maf",
  "/data2/Final_GSE156625_2.maf",
  "/data2/Final_GSE156625_SingleCell.maf",
  "/data2/Final_GSE282701_SingleCell.maf",
  "/data2/Final_SRP_SingleCell.maf"
)
META_FILE     <- "/home/download/Monopogen/monopogen_snv/cell_annotation.csv"
REDIPORTAL_FILE <- "/home/download/ref/REDIportal_hg38.txt"
OUT_DIR       <- "/home/download/csc_article/NG_refine/denovo/"

RANK_RANGE  <- 2:10
NRUN        <- 50
BEST_RANK   <- 3
SEED        <- 42

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

maf_list <- lapply(MAF_FILES, function(f) {
  if (file.exists(f)) fread(f, fill = TRUE) else NULL
})
maf_raw <- rbindlist(maf_list, fill = TRUE)
meta <- fread(META_FILE)
if (!"Tumor_Sample_Barcode" %in% names(meta)) setnames(meta, 1, "Tumor_Sample_Barcode")

maf_m <- merge(maf_raw, meta[, .(Tumor_Sample_Barcode, cell_type)],
               by = "Tumor_Sample_Barcode", all.x = FALSE)

maf_snp <- maf_m[Variant_Type == "SNP" |
                  (nchar(Reference_Allele) == 1 & nchar(Tumor_Seq_Allele2) == 1)]
maf_snp[, mut_type := paste0(Reference_Allele, ">", Tumor_Seq_Allele2)]

type_counts <- maf_snp[, .N, by = mut_type][order(-N)]
for (i in 1:nrow(type_counts)) {
  pct <- type_counts$N[i] / nrow(maf_snp) * 100
  flag <- ""
  if (type_counts$mut_type[i] %in% c("A>G", "T>C")) flag <- "  ← putative ADAR"
  message(sprintf("      %s: %d (%.1f%%)%s",
                  type_counts$mut_type[i], type_counts$N[i], pct, flag))
}

if (!is.null(REDIPORTAL_FILE) && file.exists(REDIPORTAL_FILE)) {
  redi <- fread(REDIPORTAL_FILE, select = c(2, 3), header = TRUE)
  setnames(redi, c("chrom", "pos"))

  chrom_snp <- as.character(maf_snp$Chromosome)
  if (!any(grepl("^chr", chrom_snp))) chrom_snp <- paste0("chr", chrom_snp)
  maf_snp[, key := paste0(chrom_snp, ":", Start_Position)]
  redi[, key := paste0(chrom, ":", pos)]

  n_before <- nrow(maf_snp)
  maf_snp <- maf_snp[!key %in% redi$key]
  n_removed <- n_before - nrow(maf_snp)
                  n_removed, n_removed / n_before * 100))
  maf_snp[, key := NULL]
}

type_counts2 <- maf_snp[, .N, by = mut_type][order(-N)]
for (i in 1:nrow(type_counts2)) {
  pct <- type_counts2$N[i] / nrow(maf_snp) * 100
  message(sprintf("      %s: %d (%.1f%%)", type_counts2$mut_type[i], type_counts2$N[i], pct))
}

chrom <- as.character(maf_snp$Chromosome)
if (!any(grepl("^chr", chrom))) chrom <- paste0("chr", chrom)
std_chr <- paste0("chr", c(1:22, "X", "Y"))

maf_gr <- GRanges(
  seqnames = chrom,
  ranges   = IRanges(start = maf_snp$Start_Position, end = maf_snp$End_Position),
  ref      = maf_snp$Reference_Allele,
  alt      = maf_snp$Tumor_Seq_Allele2,
  sample   = maf_snp$cell_type
)
genome(maf_gr) <- "hg38"
maf_gr <- maf_gr[seqnames(maf_gr) %in% std_chr]
seqlevels(maf_gr) <- seqlevelsInUse(maf_gr)

mut_mat <- mut_matrix(vcf_list = split(maf_gr, maf_gr$sample),
                      ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

fwrite(as.data.table(mut_mat, keep.rownames = "channel"),
       file.path(OUT_DIR, "mut_matrix_96channel.csv"))

# STEP 4: NMF Rank Survey ()
                min(RANK_RANGE), max(RANK_RANGE), NRUN))

mut_mat_nmf <- mut_mat + 0.0001

estimate <- nmf(mut_mat_nmf, rank = RANK_RANGE, method = "brunet",
                nrun = NRUN, seed = SEED, .opt = "vp")

# rank survey
pdf(file.path(OUT_DIR, "nmf_rank_survey.pdf"), width = 10, height = 7)
plot(estimate)
dev.off()

if (is.null(BEST_RANK)) {
  measures <- summary(estimate)
  coph <- measures[measures$measure == "cophenetic", ]

  coph_vals <- coph$value
  coph_ranks <- coph$rank
  if (length(coph_vals) > 2) {
    diffs <- diff(coph_vals)
    best_idx <- which.min(diffs)
    BEST_RANK <- coph_ranks[best_idx]
  } else {
    BEST_RANK <- min(RANK_RANGE) + 1
  }
} else {
}

# STEP 5:  De Novo

nmf_res <- extract_signatures(mut_mat_nmf, rank = BEST_RANK,
                               nrun = NRUN, single_core = FALSE)

# (96 x BEST_RANK)  (BEST_RANK x 12)
signatures <- nmf_res$signatures
colnames(signatures) <- paste0("DeNovo_", LETTERS[1:ncol(signatures)])

contribution <- nmf_res$contribution
rownames(contribution) <- paste0("DeNovo_", LETTERS[1:nrow(contribution)])

fwrite(as.data.table(signatures, keep.rownames = "channel"),
       file.path(OUT_DIR, "denovo_signatures.csv"))
fwrite(as.data.table(t(contribution), keep.rownames = "cell_type"),
       file.path(OUT_DIR, "denovo_contribution_matrix.csv"))

# STEP 6:  COSMIC v3.2

cancer_sigs <- get_known_signatures(source = "COSMIC_v3.2")

cos_sim_matrix <- cos_sim_matrix(signatures, cancer_sigs)

fwrite(as.data.table(cos_sim_matrix, keep.rownames = "DeNovo"),
       file.path(OUT_DIR, "cosine_similarity_denovo_vs_cosmic.csv"))

for (i in 1:ncol(signatures)) {
  sig_name <- colnames(signatures)[i]
  cos_vals <- cos_sim_matrix[i, ]
  top3_idx <- order(cos_vals, decreasing = TRUE)[1:3]
  top3_names <- colnames(cos_sim_matrix)[top3_idx]
  top3_vals  <- cos_vals[top3_idx]

  best_val <- top3_vals[1]
  interpretation <- ifelse(best_val > 0.85, "→ likely matches this COSMIC signature",
                   ifelse(best_val > 0.70, "→ partially similar, possible mixture/variant",
                                           "→  novel pattern, possible RNA artifact"))

  message(sprintf("    %s:", sig_name))
  message(sprintf("      #1 %s (cos=%.3f) %s", top3_names[1], top3_vals[1], interpretation))
  message(sprintf("      #2 %s (cos=%.3f)", top3_names[2], top3_vals[2]))
  message(sprintf("      #3 %s (cos=%.3f)", top3_names[3], top3_vals[3]))

  contrib_pct <- contribution[i, ] / colSums(contribution) * 100
  top_ct <- names(sort(contrib_pct, decreasing = TRUE))[1:3]
  top_ct_val <- sort(contrib_pct, decreasing = TRUE)[1:3]
                  top_ct[1], top_ct_val[1], top_ct[2], top_ct_val[2],
                  top_ct[3], top_ct_val[3]))
  message("")
}

reconstructed <- signatures %*% contribution
for (ct in colnames(mut_mat)) {
  cs <- crossprod(mut_mat[, ct], reconstructed[, ct]) /
    (sqrt(sum(mut_mat[, ct]^2)) * sqrt(sum(reconstructed[, ct]^2)))
  quality <- ifelse(cs > 0.95, "excellent",
             ifelse(cs > 0.90, "good",
             ifelse(cs > 0.85, "fair", "poor")))
  message(sprintf("      %s: %.4f (%s)", ct, cs, quality))
}

# STEP 7

# ---- 7a. De Novo  ----
p_sigs <- plot_96_profile(signatures, condensed = TRUE) +
  ggtitle("De Novo Mutational Signatures (NMF)")

ggsave(file.path(OUT_DIR, "denovo_signature_profiles.pdf"),
       plot = p_sigs, width = 16, height = 3 * BEST_RANK, dpi = 300)

rel_contrib <- sweep(contribution, 2, colSums(contribution), "/") * 100
rel_dt <- as.data.table(t(rel_contrib), keep.rownames = "cell_type")

plot_long <- melt(rel_dt, id.vars = "cell_type",
                  variable.name = "Signature", value.name = "Contribution")

ct_order <- c(
  "Biliary-like", "Bipotent Progenitor", "Cholangiocytes", "CSC",
  "CSC(progenitor_cycling)", "CSC(progenitor_quiet)",
  "Hepatocytes", "Hepatocytes (cycling)",
  "Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)",
  "Hepatocytes (Quiescent)", "Hepatocytes (Stem)"
)
ct_order <- ct_order[ct_order %in% unique(plot_long$cell_type)]
plot_long[, cell_type := factor(cell_type, levels = ct_order)]

denovo_colors <- brewer.pal(max(BEST_RANK, 3), "Set2")
if (BEST_RANK > length(denovo_colors)) {
  denovo_colors <- colorRampPalette(brewer.pal(8, "Set2"))(BEST_RANK)
}
names(denovo_colors) <- paste0("DeNovo_", LETTERS[1:BEST_RANK])

p_contrib <- ggplot(plot_long, aes(x = cell_type, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = denovo_colors, name = "De Novo\nSignature") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative contribution",
       title = "De Novo Signature Contributions per Cell Type") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_DIR, "denovo_contribution_barplot.pdf"),
       plot = p_contrib, width = 12, height = 6, dpi = 300)

# ---- 7c. COSMIC  ----
# de novo  cos_sim > 0.5  COSMIC
max_per_cosmic <- apply(cos_sim_matrix, 2, max)
relevant_cosmic <- names(max_per_cosmic[max_per_cosmic > 0.5])

if (length(relevant_cosmic) > 0) {
  heatmap_mat <- cos_sim_matrix[, relevant_cosmic, drop = FALSE]

  pdf(file.path(OUT_DIR, "denovo_vs_cosmic_heatmap.pdf"),
      width = max(8, length(relevant_cosmic) * 0.5),
      height = max(4, BEST_RANK * 0.8))
  pheatmap(
    heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("white", "#4DBBD5", "#E64B35"))(100),
    breaks = seq(0, 1, length.out = 101),
    fontsize = 10,
    fontsize_number = 8,
    main = "De Novo Signatures vs COSMIC v3.2 (Cosine Similarity)",
    angle_col = 45
  )
  dev.off()
}

# ---- 7d.  ----
fwrite(as.data.table(t(rel_contrib), keep.rownames = "cell_type"),
       file.path(OUT_DIR, "denovo_relative_contribution.csv"))

message("\n==================================================")
message(sprintf("  %s", OUT_DIR))
message("==================================================")

# --- 2c: Publication-quality signature plots ---

# !/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

DENOVO_DIR <- "/home/download/csc_article/NG_refine/denovo/"
OUT_DIR    <- "/home/download/csc_article/NG_refine/denovo/figures/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ct_order <- c(
  "Biliary-like", "Bipotent Progenitor", "Cholangiocytes", "CSC",
  "CSC(progenitor_cycling)", "CSC(progenitor_quiet)",
  "Hepatocytes", "Hepatocytes (cycling)",
  "Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)",
  "Hepatocytes (Quiescent)", "Hepatocytes (Stem)"
)

sig_colors <- c(
  "DeNovo_A" = "#C75B7A",
  "DeNovo_B" = "#5B8FB9",
  "DeNovo_C" = "#6DAB6D"
)

sig_labels <- c(
  "DeNovo_A" = "Sig-A (SBS7a-like, cos=0.73)",
  "DeNovo_B" = "Sig-B (SBS40-like, cos=0.55)",
  "DeNovo_C" = "Sig-C (SBS17b-like, cos=0.57)"
)

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.4, color = "grey30"),
    axis.ticks = element_line(linewidth = 0.3, color = "grey30"),
    axis.text = element_text(color = "grey20"),
    axis.title = element_text(color = "grey10"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 12, 8, 12)
  )

# denovo/  rank=5,
# 1.  run_signature_denovo_nmf.R  BEST_RANK <- 3
# 2.  Rscript run_signature_denovo_nmf.R

rel_dt <- fread(file.path(DENOVO_DIR, "denovo_relative_contribution.csv"))
sig_cols <- setdiff(names(rel_dt), "cell_type")

plot_long <- melt(rel_dt, id.vars = "cell_type",
                  variable.name = "Signature", value.name = "Contribution")
ct_order_use <- ct_order[ct_order %in% unique(plot_long$cell_type)]
plot_long[, cell_type := factor(cell_type, levels = ct_order_use)]
plot_long[, Signature := factor(Signature, levels = c("DeNovo_C", "DeNovo_B", "DeNovo_A"))]

p_contrib <- ggplot(plot_long, aes(x = cell_type, y = Contribution, fill = Signature)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = sig_colors, name = NULL,
    breaks = c("DeNovo_A", "DeNovo_B", "DeNovo_C"),
    labels = sig_labels
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.01)),
    breaks = seq(0, 100, 25)
  ) +
  labs(x = NULL, y = "Relative contribution (%)") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9.5, color = "grey15"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11.5, margin = margin(r = 8)),
    legend.position = c(0.80, 0.88),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.08, "cm")
  ) +
  guides(fill = guide_legend(reverse = FALSE, byrow = TRUE))

ggsave(file.path(OUT_DIR, "Fig_main_panelA_contribution.pdf"),
       plot = p_contrib, width = 9.5, height = 5.5, dpi = 300)

# Panel B: COSMIC

cos_dt <- fread(file.path(DENOVO_DIR, "cosine_similarity_denovo_vs_cosmic.csv"))
cos_mat <- as.matrix(cos_dt[, -1])
rownames(cos_mat) <- cos_dt$DeNovo

max_per_cosmic <- apply(cos_mat, 2, max)
relevant <- names(max_per_cosmic[max_per_cosmic > 0.35])
heatmap_mat <- cos_mat[, relevant, drop = FALSE]
col_order <- names(sort(apply(heatmap_mat, 2, max), decreasing = TRUE))
heatmap_mat <- heatmap_mat[, col_order]

rn <- rownames(heatmap_mat)
rownames(heatmap_mat) <- gsub("DeNovo_", "Sig-", rn)

heatmap_colors <- colorRampPalette(c(
  "#FFFFFF", "#F0F4F8", "#CADCE8", "#8CB4D2",
  "#4A90B8", "#2C6E99", "#D94F4F"
))(120)

pdf(file.path(OUT_DIR, "Fig_main_panelB_heatmap.pdf"), width = 10, height = 3.2)
pheatmap(
  heatmap_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f", number_color = "grey15",
  fontsize_number = 8,
  color = heatmap_colors,
  breaks = seq(0, 1, length.out = 121),
  fontsize = 10, fontsize_row = 11, fontsize_col = 9.5,
  border_color = "grey80",
  main = "",
  angle_col = 45,
  cellwidth = 32, cellheight = 26,
  legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0),
  legend_labels = c("0", "0.25", "0.50", "0.75", "1.0")
)
dev.off()

sig_dt <- fread(file.path(DENOVO_DIR, "denovo_signatures.csv"))
sig_mat <- as.matrix(sig_dt[, -1])
rownames(sig_mat) <- sig_dt$channel

p_raw <- plot_96_profile(sig_mat, condensed = TRUE)

p_profile <- p_raw +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11, face = "bold", color = "grey10",
                                margin = margin(b = 3, t = 6)),
    strip.text.y = element_text(size = 11, face = "bold", color = "grey10"),
    axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5,
                               color = "grey35", family = "mono"),
    axis.text.y = element_text(size = 8, color = "grey20"),
    axis.title.y = element_text(size = 10, color = "grey10", margin = margin(r = 6)),
    axis.title.x = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "grey40"),
    axis.ticks = element_line(linewidth = 0.2, color = "grey40"),
    axis.ticks.x = element_line(linewidth = 0.15),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(linewidth = 0.12, color = "grey90", linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave(file.path(OUT_DIR, "Fig_supp_96profile.pdf"),
       plot = p_profile, width = 15, height = 7.5, dpi = 300)

mut_dt <- fread(file.path(DENOVO_DIR, "mut_matrix_96channel.csv"))
mut_mat <- as.matrix(mut_dt[, -1])
rownames(mut_mat) <- mut_dt$channel

contrib_dt2 <- fread(file.path(DENOVO_DIR, "denovo_contribution_matrix.csv"))
contrib_mat <- as.matrix(contrib_dt2[, -1])
rownames(contrib_mat) <- contrib_dt2$cell_type
contrib_mat <- t(contrib_mat)

reconstructed <- sig_mat %*% contrib_mat

cos_sims <- sapply(colnames(mut_mat), function(ct) {
  crossprod(mut_mat[, ct], reconstructed[, ct]) /
    (sqrt(sum(mut_mat[, ct]^2)) * sqrt(sum(reconstructed[, ct]^2)))
})

recon_dt <- data.table(
  cell_type = names(cos_sims),
  cosine_similarity = as.numeric(cos_sims)
)
recon_dt[, cell_type := factor(cell_type, levels = ct_order_use)]

p_recon <- ggplot(recon_dt, aes(x = cell_type, y = cosine_similarity, fill = cosine_similarity)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#C75B7A",
             linewidth = 0.4, alpha = 0.7) +
  geom_text(aes(label = sprintf("%.3f", cosine_similarity)),
            vjust = -0.5, size = 2.8, color = "grey30") +
  scale_fill_gradient(low = "#8CB4D2", high = "#1D5B8A", guide = "none") +
  scale_y_continuous(
    limits = c(0.92, 1.008),
    breaks = seq(0.92, 1.0, 0.02),
    expand = expansion(mult = c(0, 0.05))
  ) +
  annotate("text", x = 11.8, y = 0.953, label = "0.95 threshold",
           size = 2.8, color = "#C75B7A", hjust = 1, fontface = "italic") +
  labs(x = NULL, y = "Cosine similarity\n(original vs reconstructed)") +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9, color = "grey15"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10.5, margin = margin(r = 6))
  )

ggsave(file.path(OUT_DIR, "Fig_supp_reconstruction.pdf"),
       plot = p_recon, width = 9, height = 5, dpi = 300)

# rank survey
rank_file <- file.path(DENOVO_DIR, "nmf_rank_survey.pdf")
if (file.exists(rank_file)) {
  file.copy(rank_file, file.path(OUT_DIR, "Fig_supp_rank_survey.pdf"), overwrite = TRUE)
}

message("\n==================================================")
message(sprintf("  %s", OUT_DIR))
message("         Fig_main_panelB_heatmap.pdf")
message("         Fig_supp_reconstruction.pdf")
message("         Fig_supp_rank_survey.pdf")
message("==================================================")

# --- 2b: Distribution of mutations across genomic regions ---

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

MAF_FILES_2b <- c(
  "/data2/Final_24_SingleCell.maf",
  "/data2/Final_GSE156625_2.maf",
  "/data2/Final_GSE156625_SingleCell.maf",
  "/data2/Final_GSE282701_SingleCell.maf",
  "/data2/Final_SRP_SingleCell.maf"
)
META_FILE_2b <- "/home/download/Monopogen/monopogen_snv/cell_annotation.csv"
OUT_DIR_2b   <- "/home/download/csc_article/fig3/main"
dir.create(OUT_DIR_2b, showWarnings = FALSE, recursive = TRUE)

maf_list <- lapply(MAF_FILES_2b, function(f) {
  if (file.exists(f)) fread(f, fill = TRUE) else NULL
})
maf_raw <- rbindlist(maf_list[!sapply(maf_list, is.null)], fill = TRUE)

meta <- fread(META_FILE_2b)
if (!"Tumor_Sample_Barcode" %in% names(meta)) setnames(meta, 1, "Tumor_Sample_Barcode")

maf_m <- merge(maf_raw, meta[, .(Tumor_Sample_Barcode, cell_type)],
               by = "Tumor_Sample_Barcode", all.x = FALSE)
maf_snp <- maf_m[Variant_Type == "SNP" |
                  (nchar(Reference_Allele) == 1 & nchar(Tumor_Seq_Allele2) == 1)]

# Simplify Variant_Classification
maf_snp[, Region := fcase(
  Variant_Classification %in% c("Missense_Mutation"), "Missense",
  Variant_Classification %in% c("Nonsense_Mutation"),  "Nonsense",
  Variant_Classification %in% c("Silent"),             "Synonymous",
  Variant_Classification %in% c("Splice_Site", "Splice_Region"), "Splice",
  Variant_Classification %in% c("3'UTR"),              "3'UTR",
  Variant_Classification %in% c("5'UTR"),              "5'UTR",
  Variant_Classification %in% c("3'Flank", "5'Flank"), "Flank",
  Variant_Classification %in% c("Intron"),             "Intron",
  default = "Other"
)]

# Count per cell type
region_ct <- maf_snp[, .N, by = .(cell_type, Region)]
region_ct[, Total := sum(N), by = cell_type]
region_ct[, Proportion := N / Total * 100]

# Order
ct_order_2b <- c(
  "Biliary-like", "Bipotent Progenitor", "Cholangiocytes",
  "CSC", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)",
  "Hepatocytes", "Hepatocytes (cycling)",
  "Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)",
  "Hepatocytes (Quiescent)", "Hepatocytes (Stem)"
)
ct_order_2b <- ct_order_2b[ct_order_2b %in% unique(region_ct$cell_type)]
region_ct[, cell_type := factor(cell_type, levels = ct_order_2b)]

region_order <- c("Missense", "Nonsense", "Splice", "Synonymous",
                  "3'UTR", "5'UTR", "Flank", "Intron", "Other")
region_ct[, Region := factor(Region, levels = rev(region_order))]

region_colors <- c(
  "Missense"    = "#E64B35",
  "Nonsense"    = "#4DBBD5",
  "Splice"      = "#00A087",
  "Synonymous"  = "#3C5488",
  "3'UTR"       = "#F39B7F",
  "5'UTR"       = "#8491B4",
  "Flank"       = "#91D1C2",
  "Intron"      = "#B09C85",
  "Other"       = "#D3D3D3"
)

p_region <- ggplot(region_ct, aes(x = cell_type, y = Proportion, fill = Region)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75,
           color = "white", linewidth = 0.2) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Proportion", fill = "Genomic Region") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    axis.title  = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUT_DIR_2b, "Fig2b_variant_classification.pdf"),
       plot = p_region, width = 10, height = 5.5, dpi = 300)

message("Fig 2b saved.")
