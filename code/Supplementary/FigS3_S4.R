#!/usr/bin/env Rscript
# Figure S3 — Single-cell SNV landscape, mutational signatures
# Figure S4 — CNV phylogenetic trees

# --- S3b,c,d: Mutational signature NMF + COSMIC comparison ---
# (signature extraction is in Fig2.R; plots below are the same pipeline)

# --- S4a,b,c: CNV phylogenetic trees (NJ, ggtree) ---

#!/usr/bin/env Rscript
library(phangorn); library(ggtree); library(data.table); library(dplyr); library(ggplot2); library(ape)

input_dir  <- "/home/download/tree_cnv/copykat_tree_input"
output_dir <- "/home/download/tree_cnv/tree_plots_v5"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

patients <- c("GSE242889_P5", "GSE242889_P1", "GSE242889_P2", "GSE242889_P4",
              "GSE282701_P10", "GSE242889_P3", "GSE156625_P9", "GSE282701_P12",
              "GSE149614_HCC09", "GSE282701_P14")

color_rules <- list(
  list(pattern = "CSC.*cycling",            color = "#D35E5B"),
  list(pattern = "CSC.*quiet",              color = "#C9AACB"),
  list(pattern = "CSC",                     color = "#DE8A3A"),
  list(pattern = "Bipotent",                color = "#DFB741"),
  list(pattern = "Biliary",                 color = "#F29E9A"),
  list(pattern = "Cholangiocyte",           color = "#A1779C"),
  list(pattern = "Metabolic.H|Metabolic-H", color = "#5B9B53"),
  list(pattern = "Metabolic.L|Metabolic-L", color = "#97CB7E"),
  list(pattern = "Quiescent",               color = "#A8CAE3"),
  list(pattern = "Stem",                    color = "#579B9A"),
  list(pattern = "cycling",                 color = "#94BABB"),
  list(pattern = "Hepatocyte",              color = "#5F80A5")
)

assign_color <- function(cell_type) {
  for (rule in color_rules) {
    if (grepl(rule$pattern, cell_type, ignore.case = TRUE)) return(rule$color)
  }
  return("#AAAAAA")
}

build_tree <- function(sample_id) {
  cna_file  <- file.path(input_dir, paste0(sample_id, "_CNA_mat.csv"))
  meta_file <- file.path(input_dir, paste0(sample_id, "_metadata.csv"))
  ft_file   <- file.path(input_dir, paste0(sample_id, "_final_type.csv"))
  if (!file.exists(cna_file) | !file.exists(meta_file)) { message("  [SKIP] ", sample_id); return(NULL) }

  cna_mat <- fread(cna_file, data.table = FALSE)
  rownames(cna_mat) <- cna_mat[,1]; cna_mat <- as.matrix(cna_mat[,-1])
  meta <- read.csv(meta_file, row.names = 1)
  meta$Status <- ifelse(meta$malignant == "tumor", "Tumor", "Normal")
  meta$Status[is.na(meta$Status)] <- "Tumor"

  if (file.exists(ft_file)) {
    ft <- read.csv(ft_file, row.names = 1)
    common <- intersect(rownames(meta), rownames(ft))
    if (length(common) > 0) { meta$final_type <- NA; meta[common, "final_type"] <- ft[common, "final_type"] }
  } else { meta$final_type <- meta$merged_cell_type }

  set.seed(123)
  if (nrow(meta) > 3000) {
    keep <- sample(rownames(meta), 3000)
    meta <- meta[keep, , drop = FALSE]; cna_mat <- cna_mat[, keep, drop = FALSE]
  }

  cna_char <- matrix("N", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_char[cna_mat > 0.1] <- "A"; cna_char[cna_mat < -0.1] <- "D"
  cna_char_t <- t(cna_char); rownames(cna_char_t) <- colnames(cna_mat)
  phy <- phyDat(cna_char_t, type = "USER", levels = c("N", "A", "D"))
  tree <- NJ(dist.hamming(phy))
  nr <- intersect(rownames(meta)[meta$Status == "Normal"], tree$tip.label)
  if (length(nr) > 0) { tree <- root(tree, outgroup = nr[1], resolve.root = TRUE) } else { tree <- midpoint(tree) }

  message("  Cell types in data: ", paste(sort(unique(meta$final_type)), collapse = " | "))

  return(list(tree = tree, meta = meta))
}

plot_tree <- function(tree, meta, sample_id, highlight_patterns, title_suffix = "") {
  type_col <- if ("final_type" %in% colnames(meta)) "final_type" else "merged_cell_type"
  plot_df <- data.frame(label = rownames(meta), CellType = meta[[type_col]], Status = meta$Status, stringsAsFactors = FALSE)

  plot_df$is_target <- "No"
  for (pat in highlight_patterns) plot_df$is_target[grepl(pat, plot_df$CellType, ignore.case = TRUE)] <- "Yes"
  target_types <- sort(unique(plot_df$CellType[plot_df$is_target == "Yes"]))
  if (length(target_types) == 0) { message("  [WARN] No match for: ", paste(highlight_patterns, collapse=", ")); return(NULL) }

  color_map <- sapply(target_types, assign_color)
  names(color_map) <- target_types

  message("  Highlight types & colors:")
  for (i in seq_along(color_map)) message("    ", names(color_map)[i], " -> ", color_map[i])

  p <- ggtree(tree, layout = "rectangular", size = 0) +
    geom_tree(color = "black", linewidth = 0.01, lineend = "round")

  p <- p %<+% plot_df

  p <- p + geom_tippoint(
    data = td_filter(is_target == "No"),
    color = "#EEEEEE", shape = 16, size = 0.15, alpha = 0.15
  )

  p <- p + geom_tippoint(
    data = td_filter(is_target == "Yes"),
    aes(color = CellType, shape = Status),
    size = 1.0, alpha = 0.85, stroke = 0.1
  ) +
    scale_shape_manual(
      values = c("Normal" = 16, "Tumor" = 17),
      guide = guide_legend(order = 2, override.aes = list(size = 2))
    ) +
    scale_color_manual(
      values = color_map,
      guide = guide_legend(order = 1, override.aes = list(size = 2, alpha = 1, shape = 16))
    ) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(face = "bold", size = 9, color = "#222222", margin = margin(t = 3, b = 0)),
      plot.subtitle    = element_text(face = "italic", size = 6.5, color = "#999999", margin = margin(b = 2)),
      legend.position  = "right",
      legend.key.size  = unit(0.25, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.title     = element_text(face = "bold", size = 6.5),
      legend.text      = element_text(size = 5.5),
      plot.margin      = margin(4, 2, 2, 2)
    ) +
    labs(title = sample_id, subtitle = title_suffix, color = "Cell Type", shape = "Status")
  return(p)
}

fig_w <- 4.5; fig_h <- 12

for (pid in patients) {
  message("\n=== ", pid, " ===")
  res <- tryCatch(build_tree(pid), error = function(e) { message("  [ERROR] ", e$message); NULL })
  if (is.null(res)) next

  p1 <- plot_tree(res$tree, res$meta, pid, c("CSC"), "CSC subtypes")
  if (!is.null(p1)) { ggsave(file.path(output_dir, paste0(pid, "_1_CSC.pdf")), p1, width = fig_w, height = fig_h, device = cairo_pdf); message("  1/3") }

  p2 <- plot_tree(res$tree, res$meta, pid, c("Cholangiocyte", "Hepatocyte", "Bipotent"), "Cholangiocytes / Hepatocytes / Bipotent Progenitor")
  if (!is.null(p2)) { ggsave(file.path(output_dir, paste0(pid, "_2_Chol_Hepa_Bipo.pdf")), p2, width = fig_w, height = fig_h, device = cairo_pdf); message("  2/3") }

  p3 <- plot_tree(res$tree, res$meta, pid, c("CSC", "Bipotent", "Metabolic.H|Metabolic-H|Metabolic H", "Metabolic.L|Metabolic-L|Metabolic L"), "CSC / Bipotent / Hepatocytes (Metabolic-H & L)")
  if (!is.null(p3)) { ggsave(file.path(output_dir, paste0(pid, "_3_CSC_Bipo_MetaHL.pdf")), p3, width = fig_w, height = fig_h, device = cairo_pdf); message("  3/3") }
}

message("\nDone -> ", output_dir)

# --- S3a: SNV distribution on genes (lollipop plots, e.g. NDUFA5) ---
# Note: this code is shared with Fig 2a (dN/dS analysis)

#!/usr/bin/env Rscript
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

# --- S3e: Mutation burden comparison among tumor subtypes ---
# Note: these comparisons are generated by the shared-site TMB pipeline (Fig2.R)
# and the CNV score pipeline (Fig3_cnv.py). See those scripts for details.
