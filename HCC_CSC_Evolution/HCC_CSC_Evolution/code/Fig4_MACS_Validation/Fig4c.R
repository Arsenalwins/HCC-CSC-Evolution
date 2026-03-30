#!/usr/bin/env Rscript
# Figure 4c — SNP-transcriptome integrated pseudotime (Monocle3)

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

work_dir <- "/home/download/tbsp/MACS_MP_output"
snp_base <- "/home/download/tbsp"
out_dir  <- file.path(work_dir, "monocle3_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patients <- c("MACS_Patient1", "MACS_Patient2")
OVERSAMPLE <- 50

celltype_colors <- c(
  "Biliary-like"                  = "#e8a0b4",
  "Bipotent"                      = "#f0c75e",
  "Bipotent Progenitor"           = "#f0c75e",
  "CSC"                           = "#e8943a",
  "CSC(progenitor_prolif)"        = "#d94735",
  "CSC(progenitor_quiet)"         = "#c44e52",
  "Cholangiocyte"                 = "#9b7bb8",
  "Hepatocytes"                   = "#4878a8",
  "Hepatocytes (Metabolic High)"  = "#55a868",
  "Hepatocytes (Metabolic Low)"   = "#a8bf52",
  "Hepatocytes (Cycling)"         = "#5bbfbf",
  "Progenitor-like"               = "#b0b0b0"
)

for (tag in patients) {
  cat(paste0("\n=========================================\n"))
  cat(paste0(" Monocle3: ", tag, "\n"))
  cat(paste0("=========================================\n"))

  expr_file <- file.path(work_dir, paste0("MP_Matrix_", tag, ".csv"))
  meta_file <- file.path(work_dir, paste0("MP_Meta_", tag, ".csv"))
  umap_file <- file.path(work_dir, paste0("MP_UMAP_", tag, ".csv"))

  if (!file.exists(expr_file)) { cat("  [skip] no expression matrix:", expr_file, "\n"); next }
  if (!file.exists(meta_file)) { cat("  [skip] no meta file:", meta_file, "\n"); next }
  if (!file.exists(umap_file)) { cat("  [skip] no UMAP file:", umap_file, "\n"); next }

  expr_mat  <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
  cell_meta <- read.csv(meta_file, row.names = 1)
  umap_df   <- read.csv(umap_file, row.names = 1)
  cat(sprintf("  Meta: %d cells\n", nrow(cell_meta)))
  cat(sprintf("  UMAP: %d cells (MP-based)\n", nrow(umap_df)))
  print(table(cell_meta$CellType))

  snp_file <- NULL
  tbsp_dir <- paste0("tbsp_output_", tag)
  candidates <- c(
    file.path(snp_base, "tbsp_out", gsub("Patient1", "P1", gsub("Patient2", "P2", tbsp_dir)), "SNP_matrix.tsv"),
    file.path(snp_base, "tbsp_out", tbsp_dir, "SNP_matrix.tsv"),
    file.path(snp_base, tbsp_dir, "SNP_matrix.tsv"),
    file.path(snp_base, gsub("Patient1", "P1", gsub("Patient2", "P2", tbsp_dir)), "SNP_matrix.tsv"),
    file.path(snp_base, paste0(tbsp_dir, "_all"), "SNP_matrix.tsv")
  )
  for (f in candidates) {
    if (file.exists(f)) { snp_file <- f; break }
  }
  if (is.null(snp_file)) {
    for (f in candidates) cat("    ", f, "\n")
    next
  }

  snp_df <- read.delim(snp_file, row.names = 1, check.names = FALSE)
  cat(sprintf("  SNP: %d SNPs × %d cells from %s\n", nrow(snp_df), ncol(snp_df), snp_file))

  common_cells <- Reduce(intersect, list(
    colnames(expr_mat), rownames(cell_meta), rownames(umap_df), colnames(snp_df)
  ))
  if (length(common_cells) < 10) { cat("  [skip] insufficient cells\n"); next }

  expr_mat  <- expr_mat[, common_cells, drop = FALSE]
  cell_meta <- cell_meta[common_cells, , drop = FALSE]
  umap_df   <- umap_df[common_cells, , drop = FALSE]
  snp_df    <- snp_df[, common_cells, drop = FALSE]

  snp_mat <- as.matrix(snp_df)
  n_snps <- nrow(snp_mat)
              nrow(expr_mat), n_snps, OVERSAMPLE, nrow(expr_mat) + n_snps * OVERSAMPLE))

  snp_rep_list <- lapply(1:OVERSAMPLE, function(r) {
    m <- snp_mat
    rownames(m) <- paste0(rownames(snp_mat), "_rep", r)
    m
  })
  snp_oversampled <- do.call(rbind, snp_rep_list)
  fused_mat <- rbind(expr_mat, snp_oversampled)

  valid <- colSums(fused_mat) > 0
  if (sum(!valid) > 0) {
    fused_mat <- fused_mat[, valid, drop = FALSE]
    cell_meta <- cell_meta[colnames(fused_mat), , drop = FALSE]
    umap_df   <- umap_df[colnames(fused_mat), , drop = FALSE]
  }

  fused_sparse <- as(fused_mat, "sparseMatrix")
  gene_meta <- data.frame(gene_short_name = rownames(fused_mat), row.names = rownames(fused_mat))
  cds <- new_cell_data_set(fused_sparse, cell_metadata = cell_meta, gene_metadata = gene_meta)
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 10)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")

  # ──  step3  MP-based UMAP  ──
  aligned_umap <- as.matrix(umap_df[colnames(cds), c("UMAP_1", "UMAP_2")])
  reducedDims(cds)[["UMAP"]] <- aligned_umap

  cds <- cluster_cells(cds)
  cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

  hep_types <- c("Hepatocytes", "Hepatocytes (Metabolic Low)",
                 "Hepatocytes (Metabolic High)", "Hepatocytes (Cycling)")
  hep_cells <- rownames(cell_meta)[cell_meta$CellType %in% hep_types]

  if (length(hep_cells) < 5) {
    hep_cells <- rownames(cell_meta)[order(cell_meta$Total_Mutations)[1:min(50, nrow(cell_meta))]]
  }

  hep_meta <- cell_meta[hep_cells, , drop = FALSE]
  mut_threshold <- quantile(hep_meta$Total_Mutations, 0.5)
  root_candidates <- rownames(hep_meta)[hep_meta$Total_Mutations <= mut_threshold]
  if (length(root_candidates) < 5) root_candidates <- hep_cells

  root_umap <- aligned_umap[root_candidates, , drop = FALSE]
  centroid <- colMeans(root_umap)
              length(root_candidates), centroid[1], centroid[2]))

  node_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  node_names <- rownames(node_coords)
  dists <- apply(node_coords, 1, function(x) sqrt(sum((x - centroid)^2)))
  root_node <- node_names[which.min(dists)]

  cds <- order_cells(cds, root_pr_nodes = root_node)

  rx <- centroid[1]; ry <- centroid[2]
  root_df <- data.frame(data_dim_1 = rx, data_dim_2 = ry)
  y_range <- diff(range(aligned_umap[, 2]))
  label_offset <- y_range * 0.03

  types_present <- unique(cell_meta$CellType)
  colors_use <- celltype_colors[names(celltype_colors) %in% types_present]

  # ── A: CellType + Trajectory ──
  p1 <- suppressWarnings(
    plot_cells(cds, color_cells_by = "CellType",
               label_groups_by_cluster = FALSE, label_leaves = FALSE,
               label_branch_points = FALSE, cell_size = 1.2,
               trajectory_graph_segment_size = 1.0, trajectory_graph_color = "grey30") +
    scale_color_manual(values = colors_use, na.value = "#cccccc") +
    geom_point(data = root_df, aes(x = data_dim_1, y = data_dim_2),
               shape = 21, size = 5, fill = "red", color = "black", stroke = 1.2) +
    annotate("text", x = rx, y = ry + label_offset, label = "ROOT",
             color = "red", fontface = "bold", size = 4) +
    ggtitle(paste0("Cell Type — ", tag, " (MP-based UMAP)")) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )
  ggsave(file.path(out_dir, paste0("CellType_", tag, ".pdf")), p1, width = 10, height = 7)
  ggsave(file.path(out_dir, paste0("CellType_", tag, ".png")), p1, width = 10, height = 7, dpi = 300)

  # ── B: Pseudotime ──
  p2 <- suppressWarnings(
    plot_cells(cds, color_cells_by = "pseudotime",
               label_cell_groups = FALSE, label_leaves = FALSE,
               label_branch_points = FALSE, cell_size = 1.2,
               trajectory_graph_segment_size = 1.0, trajectory_graph_color = "grey30") +
    geom_point(data = root_df, aes(x = data_dim_1, y = data_dim_2),
               shape = 21, size = 5, fill = "red", color = "black", stroke = 1.2) +
    annotate("text", x = rx, y = ry + label_offset, label = "ROOT",
             color = "red", fontface = "bold", size = 4) +
    scale_color_viridis_c() +
    ggtitle(paste0("Pseudotime — ", tag, " (MP-based UMAP)")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )
  ggsave(file.path(out_dir, paste0("Pseudotime_", tag, ".pdf")), p2, width = 10, height = 7)
  ggsave(file.path(out_dir, paste0("Pseudotime_", tag, ".png")), p2, width = 10, height = 7, dpi = 300)

  # ── C: Total Mutations ──
  p3 <- suppressWarnings(
    plot_cells(cds, color_cells_by = "Total_Mutations",
               label_cell_groups = FALSE, label_leaves = FALSE,
               label_branch_points = FALSE, cell_size = 1.2,
               trajectory_graph_segment_size = 1.0, trajectory_graph_color = "grey30") +
    scale_color_viridis_c(option = "inferno") +
    ggtitle(paste0("Total Mutations — ", tag, " (MP-based UMAP)")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )
  ggsave(file.path(out_dir, paste0("Mutations_", tag, ".pdf")), p3, width = 10, height = 7)
  ggsave(file.path(out_dir, paste0("Mutations_", tag, ".png")), p3, width = 10, height = 7, dpi = 300)

  if ("sample" %in% colnames(cell_meta)) {
    p4 <- suppressWarnings(
      plot_cells(cds, color_cells_by = "sample",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE, cell_size = 1.2) +
      ggtitle(paste0("Sample — ", tag)))
    ggsave(file.path(out_dir, paste0("Sample_", tag, ".pdf")), p4, width = 10, height = 7)
    ggsave(file.path(out_dir, paste0("Sample_", tag, ".png")), p4, width = 10, height = 7, dpi = 300)
  }

  pt <- pseudotime(cds)
  pt_df <- data.frame(
    cell = names(pt),
    pseudotime = pt,
    CellType = cell_meta[names(pt), "CellType"],
    Total_Mutations = cell_meta[names(pt), "Total_Mutations"]
  )
  write.csv(pt_df, file.path(out_dir, paste0("pseudotime_", tag, ".csv")), row.names = FALSE)

  # ── F: Pseudotime boxplot by CellType ──
  p_box <- ggplot(pt_df, aes(x = reorder(CellType, pseudotime, FUN = median),
                              y = pseudotime, fill = CellType)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = colors_use, na.value = "#cccccc") +
    coord_flip() +
    ggtitle(paste0("Pseudotime by Cell Type — ", tag)) +
    xlab("") + ylab("Pseudotime") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ggsave(file.path(out_dir, paste0("Pseudotime_boxplot_", tag, ".pdf")), p_box, width = 9, height = 6)
  ggsave(file.path(out_dir, paste0("Pseudotime_boxplot_", tag, ".png")), p_box, width = 9, height = 6, dpi = 300)

  # ── G: Pseudotime vs Mutations scatter ──
  cor_test <- cor.test(pt_df$pseudotime, pt_df$Total_Mutations, method = "spearman")
  cor_label <- sprintf("Spearman ρ = %.3f\np = %.2e", cor_test$estimate, cor_test$p.value)

  p_scatter <- ggplot(pt_df, aes(x = pseudotime, y = Total_Mutations, color = CellType)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_manual(values = colors_use, na.value = "#cccccc") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    annotate("text", x = min(pt_df$pseudotime, na.rm = TRUE),
             y = max(pt_df$Total_Mutations, na.rm = TRUE),
             label = cor_label, hjust = 0, vjust = 1, size = 4, fontface = "italic") +
    ggtitle(paste0("Pseudotime vs Mutations — ", tag)) +
    xlab("Pseudotime") + ylab("Total Mutations") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ggsave(file.path(out_dir, paste0("Pseudotime_vs_Mutations_", tag, ".pdf")), p_scatter, width = 9, height = 6)
  ggsave(file.path(out_dir, paste0("Pseudotime_vs_Mutations_", tag, ".png")), p_scatter, width = 9, height = 6, dpi = 300)

  saveRDS(cds, file.path(out_dir, paste0("cds_", tag, ".rds")))
  cat(paste0("   ", tag, " done\n"))
  rm(cds, fused_mat, fused_sparse, expr_mat, snp_mat, cell_meta, gene_meta); gc()
}

cat("\n=========================================\n")
cat(paste0(" All done! Results in: ", out_dir, "\n"))
cat("=========================================\n")
