#!/usr/bin/env Rscript
# Figure S2d — Monocle3 pseudotime (atlas normal parenchymal cells)

# Monocle3  — Normal Cells

library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

base_path <- "/data2/adata/normal_to_monocle3/"

expression_matrix <- readMM(paste0(base_path, "matrix.mtx"))
expression_matrix <- as(expression_matrix, "dgCMatrix")
cell_metadata   <- read.csv(paste0(base_path, "cells.csv"), row.names = 1)
gene_annotation <- read.csv(paste0(base_path, "genes.csv"), row.names = 1)
umap_coords     <- read.csv(paste0(base_path, "umap_coords.csv"), row.names = 1)

dimnames(expression_matrix) <- list(rownames(gene_annotation), rownames(cell_metadata))

valid_cells <- Matrix::colSums(expression_matrix) > 0
expression_matrix <- expression_matrix[, valid_cells, drop = FALSE]
cell_metadata     <- cell_metadata[valid_cells, , drop = FALSE]
umap_coords       <- umap_coords[valid_cells, , drop = FALSE]

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)
cds <- preprocess_cds(cds, num_dim = 50)

common_cells <- intersect(colnames(cds), rownames(umap_coords))
cds <- cds[, common_cells]
reducedDims(cds)[["UMAP"]] <- as.matrix(umap_coords[common_cells, ])

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds,
                   use_partition = FALSE,
                   close_loop = FALSE,
                   learn_graph_control = list(
                     rann.k = 100,
                     minimal_branch_len = 15
                   ))

# 6. Bipotent Progenitor
pr_graph <- principal_graph(cds)[["UMAP"]]
pr_node_names <- igraph::V(pr_graph)$name
pr_node_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)

bp_cells <- colnames(cds)[colData(cds)$cell_type_final == "Bipotent Progenitor"]
bp_center <- colMeans(reducedDims(cds)[["UMAP"]][bp_cells, , drop = FALSE])

dists <- sqrt(rowSums((pr_node_coords - matrix(bp_center, nrow = nrow(pr_node_coords), ncol = 2, byrow = TRUE))^2))
names(dists) <- pr_node_names
root_node <- names(which.min(dists))

cds <- order_cells(cds, root_pr_nodes = root_node)

my_colors <- c(
  "Hepatocytes"              = "#4E79A7",
  "Hepatocytes (Metabolic-H)"= "#59A14F",
  "Hepatocytes (Metabolic-L)"= "#8CD17D",
  "Hepatocytes (Stem)"       = "#499894",
  "Hepatocytes (cycling)"    = "#86BCB6",
  "Hepatocytes (Quiescent)"  = "#A0CBE8",
  "CSC(progenitor_cycling)"  = "#E15759",
  "CSC(progenitor_quiet)"    = "#D4A6C8",
  "CSC"                      = "#F28E2B",
  "Cholangiocytes"           = "#B07AA1",
  "Biliary-like"             = "#FF9D9A",
  "Bipotent Progenitor"      = "#EDC948"
)

theme_pub <- theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 8)),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

p1 <- plot_cells(cds,
                 color_cells_by = "cell_type_final",
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 cell_size = 0.6, cell_stroke = 0,
                 trajectory_graph_segment_size = 0.7,
                 trajectory_graph_color = "grey30") +
  scale_color_manual(values = my_colors) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "Cell Type & Trajectory")

p2 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 cell_size = 0.6, cell_stroke = 0,
                 trajectory_graph_segment_size = 0.7,
                 trajectory_graph_color = "grey30") +
  scale_color_viridis_c(option = "viridis", name = "Pseudotime",
                        guide = guide_colorbar(barwidth = 0.8, barheight = 8,
                                               title.position = "top")) +
  theme_pub +
  labs(title = "Pseudotime")

p3 <- plot_cells(cds,
                 color_cells_by = "cell_type_final",
                 label_groups_by_cluster = FALSE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 cell_size = 0.4, cell_stroke = 0,
                 graph_label_size = 3,
                 trajectory_graph_segment_size = 1.0,
                 trajectory_graph_color = "grey20") +
  scale_color_manual(values = my_colors) +
  theme_pub +
  theme(legend.position = "none") +
  labs(title = "Branch Points & Leaves")

p4 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 cell_size = 0.5, cell_stroke = 0,
                 trajectory_graph_segment_size = 0.5,
                 trajectory_graph_color = "grey50") +
  facet_wrap(~cell_type_final, ncol = 3) +
  scale_color_viridis_c(option = "viridis", name = "Pseudotime") +
  theme_pub +
  theme(strip.text = element_text(size = 9, face = "bold", margin = margin(b = 4)),
        legend.position = "bottom", legend.direction = "horizontal") +
  guides(color = guide_colorbar(barwidth = 12, barheight = 0.6, title.position = "left")) +
  labs(title = "Pseudotime by Cell Type")

combined <- (p1 + p2 + plot_layout(ncol = 2, widths = c(1.2, 1))) /
  p3 /
  p4 +
  plot_layout(heights = c(1, 0.9, 1.2)) +
  plot_annotation(
    title = "Monocle3 Trajectory — Normal Cells",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 12)))
  )

ggsave("/home/download/csc_article/fig2/normal.monocle.pdf",
       combined, width = 14, height = 18, dpi = 300)

print(combined)
