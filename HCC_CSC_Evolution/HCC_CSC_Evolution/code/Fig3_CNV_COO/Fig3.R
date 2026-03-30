#!/usr/bin/env Rscript
# Figure 3 — CNV accumulation and cell-of-origin analyses

# --- 3d,e: MEDALT evolutionary distance ---

#!/usr/bin/env Rscript
# MEDALT MED Score: Full Visualization & Per-Patient Pairwise Tests
# Style: paired boxplot per-patient violin

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# ======================== CONFIGURATION ========================
INPUT_DIR  <- "/home/download/csc_article/fig3/MEADLT"
OUTPUT_DIR <- "/home/download/csc_article/fig3/MEADLT/plots"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# All 17 comparison pairs
COMPARISONS <- list(
  c("CSC", "Bipotent Progenitor"),
  c("Bipotent Progenitor", "Cholangiocytes"),
  c("Bipotent Progenitor", "Hepatocytes"),
  c("Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)"),
  c("CSC", "Hepatocytes (Metabolic-H)"),
  c("CSC", "Hepatocytes (Metabolic-L)"),
  c("Bipotent Progenitor", "Hepatocytes (Metabolic-H)"),
  c("Bipotent Progenitor", "Hepatocytes (Metabolic-L)"),
  c("CSC", "CSC(progenitor_cycling)"),
  c("CSC", "CSC(progenitor_quiet)"),
  c("CSC(progenitor_cycling)", "CSC(progenitor_quiet)"),
  c("CSC(progenitor_cycling)", "Bipotent Progenitor"),
  c("CSC(progenitor_quiet)", "Bipotent Progenitor"),
  c("CSC(progenitor_cycling)", "Hepatocytes (Metabolic-H)"),
  c("CSC(progenitor_cycling)", "Hepatocytes (Metabolic-L)"),
  c("CSC(progenitor_quiet)", "Hepatocytes (Metabolic-H)"),
  c("CSC(progenitor_quiet)", "Hepatocytes (Metabolic-L)")
)

# ----------  (Okabe-Ito inspired, ) ----------
# per-patient violin
TYPE_COLORS <- c(
  "CSC"                        = "#D62728",
  "CSC(progenitor_cycling)"    = "#FF9896",
  "CSC(progenitor_quiet)"      = "#FF7F0E",
  "Bipotent Progenitor"        = "#8C564B",
  "Cholangiocytes"             = "#9467BD",
  "Hepatocytes"                = "#1F77B4",
  "Hepatocytes (Metabolic-H)"  = "#AEC7E8",
  "Hepatocytes (Metabolic-L)"  = "#98DF8A",
  "Hepatocytes (Stem)"         = "#17BECF",
  "Hepatocytes (cycling)"      = "#2CA02C",
  "Hepatocytes (Quiescent)"    = "#BCBD22",
  "Biliary-like"               = "#C5B0D5"
)

# ======================== READ DATA ========================
all_data <- read.csv(file.path(INPUT_DIR, "all_patients_med_scores.csv"),
                     stringsAsFactors = FALSE)

cat("Loaded", nrow(all_data), "cells from",
    length(unique(all_data$patient)), "patients\n")
cat("Cell types:", paste(sort(unique(all_data$final_type)), collapse = ", "), "\n\n")

safe_name <- function(s) gsub("[() ]", "_", s)

# 1. PER-PATIENT VIOLIN PLOTS ()
cat("=== 1. Generating per-patient violin plots ===\n")

patients <- sort(unique(all_data$patient))

for (pid in patients) {
  pdf_data <- all_data %>% filter(patient == pid)
  type_counts <- pdf_data %>% count(final_type) %>% filter(n >= 3)
  pdf_data <- pdf_data %>% filter(final_type %in% type_counts$final_type)

  if (nrow(type_counts) < 2) next

  type_order <- pdf_data %>%
    group_by(final_type) %>%
    summarise(med = median(med_to_root), .groups = "drop") %>%
    arrange(med) %>%
    pull(final_type)

  pdf_data$final_type <- factor(pdf_data$final_type, levels = type_order)

  fill_cols <- TYPE_COLORS[type_order]
  fill_cols[is.na(fill_cols)] <- "#AAAAAA"

  p <- ggplot(pdf_data, aes(x = final_type, y = med_to_root, fill = final_type)) +
    geom_violin(alpha = 0.45, trim = FALSE, scale = "width", color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85,
                 fill = "white", color = "grey40", linewidth = 0.35) +
    geom_jitter(width = 0.12, size = 0.4, alpha = 0.25, color = "grey40") +
    scale_fill_manual(values = fill_cols) +
    theme_classic(base_size = 12) +
    labs(title = paste("MED Score —", pid),
         x = NULL, y = "MED to Diploid Root") +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold",
                                   size = 9, color = "grey20"),
      axis.line     = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks    = element_line(color = "grey40", linewidth = 0.3),
      legend.position = "none",
      plot.title    = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(OUTPUT_DIR, paste0(pid, "_all_types_violin.pdf")),
         p, width = max(6, length(type_order) * 0.9), height = 6)
}
cat("  Done.\n\n")

# 2. PER-COMPARISON: Split violin + Pooled boxplot (17 pairs)
cat("=== 2. Generating per-comparison plots (17 pairs) ===\n")

for (idx in seq_along(COMPARISONS)) {
  comp  <- COMPARISONS[[idx]]
  type1 <- comp[1]; type2 <- comp[2]

  comp_data <- all_data %>% filter(final_type %in% c(type1, type2))
  if (nrow(comp_data) < 10) {
    cat("  Skipping", type1, "vs", type2, "- insufficient data\n"); next
  }

  patient_both <- comp_data %>%
    group_by(patient, final_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = final_type, values_from = n, values_fill = 0)

  valid_patients <- patient_both %>%
    filter(.data[[type1]] >= 3 & .data[[type2]] >= 3) %>%
    pull(patient)

  fill_vals <- c(TYPE_COLORS[type1], TYPE_COLORS[type2])
  fill_vals[is.na(fill_vals)] <- c("#CC6677", "#6699CC")[is.na(fill_vals)]
  names(fill_vals) <- c(type1, type2)

  # --- Plot A: Split violin per patient ---
  if (length(valid_patients) >= 1) {
    plot_data <- comp_data %>%
      filter(patient %in% valid_patients) %>%
      mutate(patient    = factor(patient, levels = sort(valid_patients)),
             final_type = factor(final_type, levels = c(type1, type2)))

    p1 <- ggplot(plot_data, aes(x = patient, y = med_to_root, fill = final_type)) +
      geom_violin(alpha = 0.5, position = position_dodge(width = 0.8),
                  trim = FALSE, scale = "width") +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.8),
                   outlier.shape = NA, alpha = 0.9) +
      stat_compare_means(aes(group = final_type), method = "wilcox.test",
                         label = "p.signif", label.y.npc = 0.95, size = 3.5) +
      scale_fill_manual(values = fill_vals) +
      theme_classic(base_size = 11) +
      labs(title = paste0(idx, ". ", type1, " vs ", type2, " (per patient)"),
           x = "Patient", y = "MED to Root", fill = "Cell Type") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.line     = element_line(color = "grey40", linewidth = 0.4),
            legend.position = "top")

    fname1 <- paste0(sprintf("%02d", idx), "_split_violin_",
                     safe_name(type1), "_vs_", safe_name(type2), ".pdf")
    ggsave(file.path(OUTPUT_DIR, fname1), p1,
           width = max(7, length(valid_patients) * 1.0), height = 6)
  }

  # --- Plot B: Pooled comparison ---
  pool_data <- comp_data %>%
    mutate(final_type = factor(final_type, levels = c(type1, type2)))

  n1 <- sum(pool_data$final_type == type1)
  n2 <- sum(pool_data$final_type == type2)

  p2 <- ggplot(pool_data, aes(x = final_type, y = med_to_root, fill = final_type)) +
    geom_violin(alpha = 0.4, trim = FALSE, scale = "width", color = NA) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.85,
                 fill = "white", color = "grey40", linewidth = 0.35) +
    stat_compare_means(comparisons = list(c(type1, type2)),
                       method = "wilcox.test", label = "p.format",
                       tip.length = 0.01) +
    scale_fill_manual(values = fill_vals) +
    theme_classic(base_size = 13) +
    labs(title = paste0(idx, ". ", type1, " vs ", type2, " (pooled)"),
         subtitle = paste0("n=", n1, " vs n=", n2,
                           " from ", length(unique(pool_data$patient)), " patients"),
         x = NULL, y = "MED to Root") +
    theme(axis.text.x  = element_text(face = "bold", size = 11),
          axis.line     = element_line(color = "grey40", linewidth = 0.4),
          legend.position = "none")

  fname2 <- paste0(sprintf("%02d", idx), "_pooled_",
                   safe_name(type1), "_vs_", safe_name(type2), ".pdf")
  ggsave(file.path(OUTPUT_DIR, fname2), p2, width = 5, height = 6)

  cat("  ", idx, ".", type1, "vs", type2, "- done\n")
}
cat("\n")

# 3. PER-PATIENT PAIRWISE WILCOXON TESTS (detailed table)
cat("=== 3. Per-patient pairwise Wilcoxon tests ===\n")

all_tests <- data.frame()

for (idx in seq_along(COMPARISONS)) {
  comp  <- COMPARISONS[[idx]]
  type1 <- comp[1]; type2 <- comp[2]

  for (pid in patients) {
    pdf_data <- all_data %>% filter(patient == pid)
    vals1 <- pdf_data %>% filter(final_type == type1) %>% pull(med_to_root)
    vals2 <- pdf_data %>% filter(final_type == type2) %>% pull(med_to_root)

    row <- data.frame(
      comparison_id = idx, group1 = type1, group2 = type2, patient = pid,
      n1 = length(vals1), n2 = length(vals2),
      median1 = ifelse(length(vals1) > 0, median(vals1), NA),
      median2 = ifelse(length(vals2) > 0, median(vals2), NA),
      mean1   = ifelse(length(vals1) > 0, mean(vals1), NA),
      mean2   = ifelse(length(vals2) > 0, mean(vals2), NA),
      stringsAsFactors = FALSE
    )

    if (length(vals1) >= 3 & length(vals2) >= 3) {
      wt <- wilcox.test(vals1, vals2)
      row$statistic    <- wt$statistic
      row$p_value      <- wt$p.value
      row$significance <- ifelse(wt$p.value < 0.001, "***",
                           ifelse(wt$p.value < 0.01, "**",
                           ifelse(wt$p.value < 0.05, "*", "ns")))
      row$direction    <- ifelse(median(vals1) > median(vals2),
                                 paste0(type1, " > ", type2),
                           ifelse(median(vals1) < median(vals2),
                                  paste0(type1, " < ", type2), "="))
    } else {
      row$statistic    <- NA
      row$p_value      <- NA
      row$significance <- "insufficient"
      row$direction    <- NA
    }
    all_tests <- rbind(all_tests, row)
  }
}

write.csv(all_tests, file.path(OUTPUT_DIR, "all_pairwise_tests_per_patient.csv"),
          row.names = FALSE)

# Print summary
cat("\nSummary: significant patients per comparison\n")
summary_sig <- all_tests %>%
  filter(!is.na(p_value)) %>%
  group_by(comparison_id, group1, group2) %>%
  summarise(
    n_patients  = n(),
    n_sig_005   = sum(p_value < 0.05),
    n_sig_001   = sum(p_value < 0.01),
    n_sig_0001  = sum(p_value < 0.001),
    median_p    = median(p_value),
    .groups     = "drop"
  )

for (i in 1:nrow(summary_sig)) {
  r <- summary_sig[i, ]
  cat(sprintf("  %2d. %-30s vs %-30s | %2d patients | sig(.05):%2d | sig(.01):%2d | sig(.001):%2d\n",
              r$comparison_id, r$group1, r$group2,
              r$n_patients, r$n_sig_005, r$n_sig_001, r$n_sig_0001))
}

write.csv(summary_sig, file.path(OUTPUT_DIR, "summary_significance_per_comparison.csv"),
          row.names = FALSE)

# 4. P-VALUE HEATMAP (patients x comparisons)
cat("\n=== 4. Generating p-value heatmap ===\n")

heatmap_data <- all_tests %>%
  filter(!is.na(p_value)) %>%
  mutate(
    comparison = paste0(comparison_id, ". ", group1, "\nvs ", group2),
    log10p     = pmin(-log10(p_value), 10),
    sig_label  = significance
  )

if (nrow(heatmap_data) > 0) {
  comp_order <- heatmap_data %>%
    distinct(comparison_id, comparison) %>%
    arrange(comparison_id) %>%
    pull(comparison)
  heatmap_data$comparison <- factor(heatmap_data$comparison, levels = comp_order)

  p_hm <- ggplot(heatmap_data, aes(x = patient, y = comparison, fill = log10p)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig_label), size = 2.5) +
    scale_fill_gradient2(low = "white", mid = "#FFFFCC", high = "#D73027",
                         midpoint = 1.3, name = "-log10(p)",
                         limits = c(0, 10)) +
    theme_minimal(base_size = 10) +
    labs(title = "Per-patient significance of MED comparisons",
         x = "Patient", y = "") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_blank())

  ggsave(file.path(OUTPUT_DIR, "pvalue_heatmap_all_comparisons.pdf"), p_hm,
         width = max(14, length(patients) * 0.5),
         height = max(10, length(COMPARISONS) * 0.6))
  cat("  Saved.\n")
}

# 5. FOREST PLOT: pooled effect per comparison
cat("\n=== 5. Generating forest plot ===\n")

forest_data <- data.frame()
for (idx in seq_along(COMPARISONS)) {
  comp  <- COMPARISONS[[idx]]
  type1 <- comp[1]; type2 <- comp[2]
  vals1 <- all_data %>% filter(final_type == type1) %>% pull(med_to_root)
  vals2 <- all_data %>% filter(final_type == type2) %>% pull(med_to_root)

  if (length(vals1) >= 3 & length(vals2) >= 3) {
    wt <- wilcox.test(vals1, vals2)
    forest_data <- rbind(forest_data, data.frame(
      idx = idx,
      comparison  = paste0(idx, ". ", type1, "\nvs ", type2),
      median_diff = median(vals1) - median(vals2),
      mean_diff   = mean(vals1) - mean(vals2),
      n1 = length(vals1), n2 = length(vals2),
      p_value = wt$p.value,
      sig = ifelse(wt$p.value < 0.001, "***",
              ifelse(wt$p.value < 0.01, "**",
              ifelse(wt$p.value < 0.05, "*", "ns")))
    ))
  }
}

if (nrow(forest_data) > 0) {
  forest_data$comparison <- factor(forest_data$comparison,
                                    levels = rev(forest_data$comparison))

  p_forest <- ggplot(forest_data, aes(x = mean_diff, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = n1 + n2, color = -log10(p_value)), shape = 16) +
    geom_text(aes(label = sig), hjust = -0.5, size = 4) +
    scale_color_gradient(low = "grey70", high = "red", name = "-log10(p)") +
    scale_size_continuous(range = c(2, 8), name = "Total cells") +
    theme_minimal(base_size = 11) +
    labs(title = "Pooled MED difference (group1 − group2)",
         x = "Mean MED difference", y = "") +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(OUTPUT_DIR, "forest_plot_pooled.pdf"), p_forest,
         width = 11, height = max(7, nrow(forest_data) * 0.45))
  cat("  Saved.\n")
}

# 6. OVERVIEW: all key types in one plot (pooled)
cat("\n=== 6. Generating overview plot ===\n")

key_types <- c("CSC", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)",
               "Bipotent Progenitor", "Cholangiocytes", "Hepatocytes",
               "Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)")

overview_data <- all_data %>% filter(final_type %in% key_types)

if (nrow(overview_data) > 0) {
  type_order <- overview_data %>%
    group_by(final_type) %>%
    summarise(med = median(med_to_root), .groups = "drop") %>%
    arrange(med) %>%
    pull(final_type)

  overview_data$final_type <- factor(overview_data$final_type, levels = type_order)

  fill_cols <- TYPE_COLORS[type_order]
  fill_cols[is.na(fill_cols)] <- "#AAAAAA"

  p_ov <- ggplot(overview_data, aes(x = final_type, y = med_to_root,
                                     fill = final_type)) +
    geom_violin(alpha = 0.45, trim = FALSE, scale = "width", color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85,
                 fill = "white", color = "grey40", linewidth = 0.35) +
    scale_fill_manual(values = fill_cols) +
    theme_classic(base_size = 13) +
    labs(title = "MED Score Overview — All Patients Pooled",
         subtitle = paste0("N = ", nrow(overview_data), " cells, ",
                           length(unique(overview_data$patient)), " patients"),
         x = NULL, y = "MED to Diploid Root") +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, face = "bold",
                                   size = 10, color = "grey20"),
      axis.line     = element_line(color = "grey40", linewidth = 0.4),
      legend.position = "none",
      plot.title    = element_text(face = "bold", size = 14)
    )

  ggsave(file.path(OUTPUT_DIR, "overview_all_key_types.pdf"), p_ov,
         width = 10, height = 8)
  cat("  Saved.\n")
}

# 7. PATIENT-LEVEL PAIRED COMPARISONS —
cat("\n=== 7. Patient-level paired comparisons (reference-style) ===\n")

median_csv <- file.path(INPUT_DIR, "patient_type_median_med.csv")
if (file.exists(median_csv)) {
  patient_med <- read.csv(median_csv, stringsAsFactors = FALSE)
  cat("  Loaded pre-computed patient-level medians from CSV\n")
} else {
  patient_med <- all_data %>%
    group_by(patient, final_type) %>%
    summarise(
      median_med = median(med_to_root),
      mean_med   = mean(med_to_root),
      n_cells    = n(),
      .groups    = "drop"
    )
  cat("  Computed patient-level medians from raw data\n")
}

paired_results <- data.frame()

for (idx in seq_along(COMPARISONS)) {
  comp  <- COMPARISONS[[idx]]
  type1 <- comp[1]; type2 <- comp[2]

  df1 <- patient_med %>%
    filter(final_type == type1, n_cells >= 3) %>%
    select(patient, median1 = median_med, mean1 = mean_med, n1 = n_cells)
  df2 <- patient_med %>%
    filter(final_type == type2, n_cells >= 3) %>%
    select(patient, median2 = median_med, mean2 = mean_med, n2 = n_cells)

  paired   <- inner_join(df1, df2, by = "patient")
  n_pairs  <- nrow(paired)

  if (n_pairs < 3) {
    cat("  ", idx, ".", type1, "vs", type2, "- only", n_pairs,
        "paired patients, skipping\n")
    paired_results <- rbind(paired_results, data.frame(
      comparison_id = idx, group1 = type1, group2 = type2,
      n_patients = n_pairs,
      median_of_median1 = NA, median_of_median2 = NA,
      statistic = NA, p_value = NA, significance = "insufficient",
      stringsAsFactors = FALSE
    ))
    next
  }

  wt <- tryCatch(
    wilcox.test(paired$median1, paired$median2, paired = TRUE),
    error = function(e) list(statistic = NA, p.value = NA)
  )
  pval <- wt$p.value
  sig  <- ifelse(is.na(pval), "NA",
           ifelse(pval < 0.001, "***",
           ifelse(pval < 0.01, "**",
           ifelse(pval < 0.05, "*", "ns"))))

  paired_results <- rbind(paired_results, data.frame(
    comparison_id = idx, group1 = type1, group2 = type2,
    n_patients = n_pairs,
    median_of_median1 = median(paired$median1),
    median_of_median2 = median(paired$median2),
    statistic = ifelse(is.null(wt$statistic), NA, wt$statistic),
    p_value = pval, significance = sig,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %2d. %-30s vs %-30s | %2d pairs | p=%.2e %s\n",
              idx, type1, type2, n_pairs, pval, sig))

  # ·  shape=21
  plot_df <- paired %>%
    select(patient, median1, median2) %>%
    pivot_longer(cols = c(median1, median2),
                 names_to = "group", values_to = "median_med") %>%
    mutate(
      cell_type = factor(
        ifelse(group == "median1", type1, type2),
        levels = c(type1, type2)
      )
    )

  fill_vals <- c(TYPE_COLORS[type1], TYPE_COLORS[type2])
  fill_vals[is.na(fill_vals)] <- c("#CC6677", "#6699CC")[is.na(fill_vals)]
  names(fill_vals) <- c(type1, type2)

  if (is.na(pval)) {
    p_sig_label <- "ns"
    p_num_label <- "P = NA"
  } else if (pval < 0.001) {
    p_sig_label <- "***"
    p_num_label <- sprintf("italic(P) == %.4f", pval)
  } else if (pval < 0.01) {
    p_sig_label <- "**"
    p_num_label <- sprintf("italic(P) == %.4f", pval)
  } else if (pval < 0.05) {
    p_sig_label <- "*"
    p_num_label <- sprintf("italic(P) == %.4f", pval)
  } else {
    p_sig_label <- "ns"
    p_num_label <- sprintf("italic(P) == %.3f", pval)
  }

  y_max   <- max(plot_df$median_med, na.rm = TRUE)
  y_range <- diff(range(plot_df$median_med, na.rm = TRUE))
  bar_y   <- y_max + y_range * 0.08
  sig_y   <- y_max + y_range * 0.18
  pnum_y  <- y_max + y_range * 0.13

  p_paired <- ggplot(plot_df, aes(x = cell_type, y = median_med)) +
    geom_line(aes(group = patient), color = "grey75", alpha = 0.6,
              linewidth = 0.35) +
    geom_boxplot(aes(fill = cell_type), width = 0.45, alpha = 0.55,
                 outlier.shape = NA, color = "black", linewidth = 0.4,
                 fatten = 1.8) +
    geom_point(aes(fill = cell_type), shape = 21, size = 2.5,
               alpha = 0.75, color = "grey30", stroke = 0.3,
               position = position_jitter(width = 0.08, seed = 42)) +
    annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y,
             color = "black", linewidth = 0.45) +
    annotate("text", x = 1.5, y = sig_y, label = p_sig_label,
             size = 5, fontface = "bold") +
    annotate("text", x = 1.5, y = pnum_y,
             label = p_num_label, parse = TRUE,
             size = 3.2, color = "grey20") +
    scale_fill_manual(values = fill_vals) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    theme_classic(base_size = 12) +
    labs(
      title    = NULL,
      x        = NULL,
      y        = "MED to Diploid Root"
    ) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(face = "bold", size = 11,
                                        color = "grey10"),
      axis.title.y      = element_text(size = 11),
      axis.line         = element_line(color = "grey30", linewidth = 0.4),
      axis.ticks        = element_line(color = "grey30", linewidth = 0.3),
      plot.margin       = margin(t = 15, r = 10, b = 5, l = 10)
    )

  fname_paired <- paste0(sprintf("%02d", idx), "_paired_",
                         safe_name(type1), "_vs_", safe_name(type2), ".pdf")
  ggsave(file.path(OUTPUT_DIR, fname_paired), p_paired,
         width = 4.2, height = 5.5)

  diff_df <- paired %>%
    mutate(diff = median1 - median2) %>%
    arrange(diff) %>%
    mutate(patient = factor(patient, levels = patient))

  p_label_diff <- ifelse(is.na(pval), "P = NA",
                   ifelse(pval < 0.001, sprintf("P = %.1e ***", pval),
                   ifelse(pval < 0.01,  sprintf("P = %.4f **", pval),
                   ifelse(pval < 0.05,  sprintf("P = %.4f *", pval),
                                         sprintf("P = %.3f ns", pval)))))

  p_diff <- ggplot(diff_df, aes(x = patient, y = diff, fill = diff > 0)) +
    geom_col(alpha = 0.65, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(
      values = c("TRUE" = fill_vals[type1], "FALSE" = fill_vals[type2]),
      labels = c("TRUE"  = paste(type1, ">", type2),
                 "FALSE" = paste(type2, ">", type1)),
      name = "Direction"
    ) +
    annotate("text", x = nrow(diff_df) / 2,
             y = max(abs(diff_df$diff)) * 1.15,
             label = p_label_diff, size = 3.5, fontface = "bold") +
    theme_classic(base_size = 11) +
    labs(
      title    = paste0(idx, ". Patient-level MED difference"),
      subtitle = paste0(type1, " median − ", type2, " median"),
      x = "Patient", y = expression(Delta ~ "Median MED")
    ) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.line        = element_line(color = "grey40", linewidth = 0.4),
      legend.position  = "top"
    )

  fname_diff <- paste0(sprintf("%02d", idx), "_paired_diff_",
                        safe_name(type1), "_vs_", safe_name(type2), ".pdf")
  ggsave(file.path(OUTPUT_DIR, fname_diff), p_diff,
         width = max(5, n_pairs * 0.6), height = 5)
}

write.csv(paired_results,
          file.path(OUTPUT_DIR, "paired_wilcoxon_signed_rank_tests.csv"),
          row.names = FALSE)

cat("\n============================================\n")
cat("All outputs saved to:", OUTPUT_DIR, "\n")
cat("Key files:\n")
cat("  - all_pairwise_tests_per_patient.csv\n")
cat("  - summary_significance_per_comparison.csv\n")
cat("  - pvalue_heatmap_all_comparisons.pdf\n")
cat("  - forest_plot_pooled.pdf\n")
cat("  - overview_all_key_types.pdf\n")
cat("  - 01-17_split_violin_*.pdf\n")
cat("  - 01-17_pooled_*.pdf\n")
cat("  - <patient>_all_types_violin.pdf\n")
cat("  - paired_wilcoxon_signed_rank_tests.csv\n")
cat("  - XX_paired_*.pdf   (reference-style paired boxplot)\n")
cat("  - XX_paired_diff_*.pdf\n")

# --- 3f,g: Cell-of-origin (TCR-based negative binomial regression) ---

# COO  — Panja et al. 2025 Nature Genetics
# Barplot: Relative Risk (1),

OUTDIR <- "/home/download/csc_article/fig3/COO/"
INDIR  <- "/home/download/COO/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTDIR)

suppressPackageStartupMessages({
  library(MASS); library(data.table); library(dplyr); library(readxl)
  library(readr); library(ggplot2); library(GenomeInfoDb); library(GenomicRanges)
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(MutationalPatterns)
  library(org.Hs.eg.db); library(rtracklayer); library(ComplexHeatmap)
  library(circlize); library(parallel)
})

# ║  PART 1                                                                    ║
case_raw <- read_excel(paste0(INDIR, "Cases_20260117.xlsx"))
clin <- case_raw %>%
  mutate(
    HBV_Group = case_when(
      is.na(Hepatitis)~"Unknown",
      grepl("HBV",Hepatitis,ignore.case=TRUE)&!grepl("HBV-HCV-",Hepatitis,fixed=TRUE)~"HBV_Pos",
      TRUE~"HBV_Neg"),
    Virus_Status = case_when(
      grepl("HBV",Hepatitis)&grepl("HCV",Hepatitis)~"HBV+HCV",
      grepl("HBV",Hepatitis)&!grepl("HBV-",Hepatitis)~"HBV",
      grepl("HCV",Hepatitis)&!grepl("HCV-",Hepatitis)~"HCV",
      TRUE~"NBNC")
  ) %>% dplyr::select(CaseID, HBV_Group, Virus_Status, Gender, everything())
fwrite(as.data.table(clin), paste0(OUTDIR, "clinical_info.csv"))

mut_raw <- read_excel(paste0(INDIR, "Mutations_20260117.xlsx"))
mut_snp <- mut_raw %>% filter(Type=="SNP", nchar(Ref)==1, nchar(Allele)==1)
mut_gr_df <- as.data.frame(mut_snp) %>%
  dplyr::rename(chromosome=Chr,start=Start,end=End,REF=Ref,ALT=Allele) %>%
  mutate(chromosome=paste0("chr",sub("^chr","",chromosome)))
gr_all <- makeGRangesFromDataFrame(mut_gr_df,keep.extra.columns=TRUE,ignore.strand=TRUE,
  seqnames.field="chromosome",start.field="start",end.field="end")
genome(gr_all)<-"hg19"; seqlevelsStyle(gr_all)<-"UCSC"
gr_all <- keepStandardChromosomes(gr_all,pruning.mode="coarse")

# ║  PART 2: APOBEC                                                            ║
vcfs_by_patient <- split(gr_all, gr_all$CaseID)
mut_mat <- mut_matrix(vcfs_by_patient, "BSgenome.Hsapiens.UCSC.hg19")
cosmic_all <- get_known_signatures(source="COSMIC_v3.2")
fit_res <- fit_to_signatures(mut_mat, cosmic_all)
sig_contrib <- fit_res$contribution

# Signature counts
sig_counts <- as.data.frame(t(sig_contrib))
sig_counts$CaseID <- rownames(sig_counts)
fwrite(sig_counts, paste0(OUTDIR, "Signature_Counts.csv"))

# Motif-based APOBEC (HCC  posterior ,  motif)
tc_res <- type_context(gr_all, "BSgenome.Hsapiens.UCSC.hg19")
is_apobec <- (tc_res$types %in% c("C>T","C>G")) &
  (substr(tc_res$context,1,1)=="T") &
  (substr(tc_res$context,3,3) %in% c("A","T"))
message(sprintf("   APOBEC (TCW motif): %d (%.1f%%)", sum(is_apobec), 100*mean(is_apobec)))

gr_clean <- gr_all[!is_apobec]
mut_clean <- as.data.table(as.data.frame(gr_clean))
if ("seqnames" %in% colnames(mut_clean)) setnames(mut_clean,"seqnames","chromosome")

# ║  PART 3:                                                             ║
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; genes_gr <- genes(txdb)
gene_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, genes_gr)
gc <- letterFrequency(gene_seqs, letters="GC", as.prob=TRUE)[,"G|C"]
exons_reduced <- reduce(exonsBy(txdb, by="gene"))
exon_len <- sum(width(exons_reduced)); exon_len <- exon_len[names(genes_gr)]
exon_len[is.na(exon_len)] <- 0; gene_w <- width(genes_gr)

bw_file <- paste0(INDIR, "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig")
if (!file.exists(bw_file)) download.file(paste0(
  "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/",
  "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig"), destfile=bw_file, mode="wb")
rt_track <- import(bw_file)
ov <- findOverlaps(genes_gr, rt_track)
avg_rt <- tapply(score(rt_track)[subjectHits(ov)], queryHits(ov), mean)
rt_vec <- rep(0, length(genes_gr)); rt_vec[as.numeric(names(avg_rt))] <- avg_rt

covariates <- data.table(entrez_id=names(genes_gr), width=as.numeric(gene_w),
  exon_width=as.numeric(exon_len), gc=as.numeric(gc),
  exon_frac=as.numeric(exon_len/gene_w), rept=rt_vec)
symbols <- mapIds(org.Hs.eg.db, keys=covariates$entrez_id, column="SYMBOL",
                  keytype="ENTREZID", multiVals="first")
covariates$Gene <- symbols
covariates <- covariates[!is.na(Gene)][order(-width)][!duplicated(Gene)]
fwrite(covariates, paste0(OUTDIR, "gene_covariates.csv"))

# ║  PART 4:                                                           ║
run_coo_regression <- function(mut_dt, covariates, centroids, hvg_list, cell_types, label="") {
  mut_agg <- mut_dt[, .(mutation_count=.N), by=Gene]
  dat <- mut_agg %>% inner_join(covariates,by="Gene") %>%
    inner_join(centroids,by="Gene") %>% filter(Gene %in% hvg_list) %>% as.data.table()
  if (nrow(dat)<50) return(NULL)
  n_mut <- sum(dat$mutation_count)
  
  expr_mat <- as.matrix(dat[, ..cell_types]); avg_expr <- rowMeans(expr_mat)
  base <- data.frame(y=dat$mutation_count, log_avg=log(avg_expr+1),
    gc=dat$gc, rept=dat$rept, log_exon_w=log(dat$exon_width+1))
  
  results <- data.frame()
  for (ct in cell_types) {
    md <- base; md$log_expr <- log(dat[[ct]]+1)
    fit <- tryCatch(glm.nb(y~log_expr+log_avg+gc+rept+offset(log_exon_w),
      data=md, control=glm.control(maxit=100)), error=function(e) NULL)
    if (is.null(fit)) next
    coefs <- summary(fit)$coefficients
    if (!"log_expr" %in% rownames(coefs)) next
    est <- coefs["log_expr","Estimate"]; se <- coefs["log_expr","Std. Error"]
    results <- rbind(results, data.frame(CellType=ct, Estimate=est, StdError=se,
      PValue=coefs["log_expr","Pr(>|z|)"],
      RelativeRisk=exp(est), RR_Lower=exp(est-1.96*se), RR_Upper=exp(est+1.96*se),
      stringsAsFactors=FALSE))
  }
  if (nrow(results)==0) return(NULL)
  results$FDR <- p.adjust(results$PValue, method="BH")
  
  sig_neg <- results %>% filter(FDR<0.1 & Estimate<0)
  coo_set <- character(0)
  if (nrow(sig_neg)>0) {
    top <- sig_neg[which.min(sig_neg$Estimate),]
    coo_set <- top$CellType
    tw <- top$RR_Upper - top$RR_Lower
    for (i in seq_len(nrow(sig_neg))) {
      if (sig_neg$CellType[i]==top$CellType) next
      ov_w <- max(0, min(top$RR_Upper,sig_neg$RR_Upper[i])-max(top$RR_Lower,sig_neg$RR_Lower[i]))
      if (tw>0 && ov_w/tw>=0.6) coo_set <- c(coo_set, sig_neg$CellType[i])
    }
  }
  results$IsCOO <- results$CellType %in% coo_set
  results$Significance <- case_when(
    results$IsCOO~"Inferred COO", results$FDR<0.1&results$Estimate<0~"Sig. negative",
    results$FDR<0.1&results$Estimate>0~"Sig. positive", TRUE~"Not significant")
  
  message(sprintf("   [%s] COO: %s", label,
    ifelse(length(coo_set)>0, paste(coo_set,collapse=", "), "none")))
  attr(results,"label")<-label; attr(results,"n_mutations")<-n_mut; attr(results,"n_genes")<-nrow(dat)
  return(results)
}

# ║  PART 5: Aggregate                                                         ║
message("\n>> PART 5: Aggregate COO")
cent_atlas <- fread(paste0(OUTDIR,"liver_centroids_atlas.csv"))
hvg_atlas <- fread(paste0(OUTDIR,"hvg_atlas.csv"))$Gene
ct_atlas <- setdiff(colnames(cent_atlas),"Gene")

cent_hepcom <- fread(paste0(OUTDIR,"liver_centroids_hepcom.csv"))
hvg_hepcom <- fread(paste0(OUTDIR,"hvg_hepcom.csv"))$Gene
ct_hepcom <- setdiff(colnames(cent_hepcom),"Gene")

message("\n--- Atlas ---")
res_atlas <- run_coo_regression(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,"Atlas")
message("\n--- HepCom ---")
res_hepcom <- run_coo_regression(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,"HepCom")

for (res in list(res_atlas,res_hepcom)) {
  if (!is.null(res)) write.csv(res, paste0(OUTDIR,"COO_Result_",attr(res,"label"),".csv"), row.names=FALSE)
}

# ║  PART 6: Per-Patient                                                       ║
message("\n>> PART 6: Per-Patient")
run_per_patient <- function(mut_dt, covariates, centroids, hvg_list, cell_types, min_snv=80, label="") {
  shared <- Reduce(intersect, list(hvg_list, covariates$Gene, centroids$Gene))
  ref <- merge(covariates[Gene %in% shared], centroids[Gene %in% shared], by="Gene")
  per_gene <- mut_dt[Gene %in% shared, .(N=.N), by=.(CaseID,Gene)]
  valid_pid <- per_gene[,.(Tot=sum(N)),by=CaseID][Tot>=min_snv]$CaseID
  if (length(valid_pid)<5) return(NULL)
  run_one <- function(pid) {
    df <- copy(ref)
    df <- merge(df, per_gene[CaseID==pid,.(Gene,N)], by="Gene", all.x=TRUE); df[is.na(N),N:=0]
    avg_e <- rowMeans(as.matrix(df[,..cell_types]))
    df$log_avg <- log(avg_e+1); df$log_exon_w <- log(df$exon_width+1)
    beta <- setNames(rep(NA_real_,length(cell_types)),cell_types)
    for (ct in cell_types) {
      df$log_expr <- log(df[[ct]]+1)
      fit <- tryCatch(glm.nb(N~log_expr+log_avg+gc+rept+offset(log_exon_w),
        data=df,control=glm.control(maxit=50)),error=function(e) NULL)
      if (!is.null(fit)&&"log_expr" %in% names(coef(fit))) beta[ct] <- coef(fit)["log_expr"]
    }; return(beta)
  }
  nc <- min(detectCores()-2,8)
  blist <- mclapply(valid_pid, run_one, mc.cores=nc)
  bmat <- do.call(rbind,blist); rownames(bmat) <- valid_pid
  bmat <- bmat[complete.cases(bmat),,drop=FALSE]
}
beta_atlas <- run_per_patient(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,80,"Atlas")
beta_hepcom <- run_per_patient(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,80,"HepCom")
save(beta_atlas,beta_hepcom,sig_counts,file=paste0(OUTDIR,"PerPatient_Betas.RData"))

# ║  PART 7:  (Nature Genetics )                                       ║
clin <- fread(paste0(OUTDIR,"clinical_info.csv"))
sig_counts_dt <- fread(paste0(OUTDIR,"Signature_Counts.csv"))

# 7A. Aggregate Barplot — Relative Risk ( 1 ,  Fig2 )
plot_coo_bar <- function(res, fname) {
  if (is.null(res)) return(invisible(NULL))
  lbl <- attr(res,"label")
  
  cols <- c("Inferred COO"="#C0392B", "Sig. negative"="#E74C3C",
            "Sig. positive"="#2980B9", "Not significant"="#BDC3C7")
  
  res$LogRR    <- res$Estimate
  res$CI_lower <- res$Estimate - 1.96 * res$StdError
  res$CI_upper <- res$Estimate + 1.96 * res$StdError
  
  max_abs <- max(abs(c(res$CI_lower, res$CI_upper)), na.rm=TRUE) * 1.1
  step <- if (max_abs > 0.5) 0.2 else if (max_abs > 0.25) 0.1 else 0.05
  brks <- seq(-ceiling(max_abs/step)*step, ceiling(max_abs/step)*step, by=step)
  
  p <- ggplot(res, aes(x=reorder(CellType, LogRR), y=LogRR, fill=Significance)) +
    geom_bar(stat="identity", width=0.7, alpha=0.9) +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.25, linewidth=0.5, color="grey30") +
    scale_fill_manual(values=cols, drop=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.6) +
    coord_flip() +
    scale_y_continuous(breaks=brks, expand=expansion(mult=0.02)) +
    theme_bw(base_size=12) +
    labs(title=sprintf("Cell-of-Origin: %s", lbl),
         subtitle=sprintf("n = %s mutations | %d genes | APOBEC excluded",
                          format(attr(res,"n_mutations"),big.mark=","), attr(res,"n_genes")),
         x="", y="Log relative risk",
         caption="< 0 = TCR protection (potential origin) | > 0 = positive correlation") +
    theme(axis.text.y=element_text(size=10, face="bold", color="black"),
          axis.text.x=element_text(size=9),
          plot.title=element_text(hjust=0.5, face="bold", size=14),
          panel.grid.major.y=element_blank(),
          legend.position="bottom")
  
  ggsave(fname, p, width=9, height=max(6, length(unique(res$CellType))*0.4))
  message(sprintf("   %s", fname))
}

plot_coo_bar(res_atlas,  paste0(OUTDIR,"Fig_COO_Atlas.pdf"))
plot_coo_bar(res_hepcom, paste0(OUTDIR,"Fig_COO_HepCom.pdf"))

# 7B. Per-Patient Heatmap — Nature Fig 3a  ()
plot_heatmap <- function(bmat, lbl, fname) {
  if (is.null(bmat)||nrow(bmat)<5) { message(sprintf("   [%s] skipping heatmap",lbl)); return(invisible(NULL)) }
  
  # beta (= log RR),  0
  pmat <- t(bmat)
  
  adf <- data.frame(CaseID=colnames(pmat), stringsAsFactors=FALSE) %>%
    left_join(clin, by="CaseID") %>%
    left_join(sig_counts_dt, by="CaseID")
  
  tmb_per_patient <- mut_clean[, .N, by=CaseID]
  adf <- adf %>% left_join(as.data.frame(tmb_per_patient) %>% dplyr::rename(TMB=N), by="CaseID")
  adf$TMB[is.na(adf$TMB)] <- 0
  
  # Signature : Aging(SBS1+SBS5), APOBEC(SBS2+SBS13), AA(SBS22), Aflatoxin(SBS24)
  adf$Aging   <- ifelse("SBS1" %in% colnames(adf) & "SBS5" %in% colnames(adf),
                        adf$SBS1 + adf$SBS5, 0)
  adf$APOBEC  <- ifelse("SBS2" %in% colnames(adf) & "SBS13" %in% colnames(adf),
                        adf$SBS2 + adf$SBS13, 0)
  adf$AA      <- ifelse("SBS22" %in% colnames(adf), adf$SBS22, 0)
  adf$Aflatoxin <- ifelse("SBS24" %in% colnames(adf), adf$SBS24, 0)
  
  col_virus  <- c("HBV"="#d95f02","HCV"="#7570b3","HBV+HCV"="#e7298a","NBNC"="#66a61e")
  col_gender <- c("Male"="#1f78b4","Female"="#fb9a99")
  col_hbv    <- c("HBV_Pos"="#d95f02","HBV_Neg"="#66a61e","Unknown"="grey85")
  
  col_sig <- c("AA"="#984ea3", "APOBEC"="#e41a1c", "Aflatoxin"="#ffff33", "Aging"="#bababa")
  
  ha_top <- HeatmapAnnotation(
    "TMB" = anno_barplot(adf$TMB, gp=gpar(fill="black",col=NA),
                         height=unit(1.2,"cm"), border=FALSE,
                         axis_param=list(side="right",gp=gpar(fontsize=7))),
    "Signatures" = anno_barplot(
      as.matrix(adf[, c("AA","Aflatoxin","APOBEC","Aging")]),
      gp=gpar(fill=col_sig, col=NA), height=unit(1.5,"cm"), border=FALSE),
    "Virus"  = adf$Virus_Status,
    "HBV"    = adf$HBV_Group,
    "Gender" = adf$Gender,
    col = list("Virus"=col_virus, "HBV"=col_hbv, "Gender"=col_gender),
    na_col = "grey90",
    gap = unit(c(1,2,1,1), "mm"),
    annotation_name_side = "left",
    annotation_label = c("TMB","Signatures","Virus","HBV","Gender")
  )
  
  # log(RR)  0  (=teal=COO, =brown)
  lim <- max(abs(quantile(pmat, c(0.02,0.98), na.rm=TRUE)))
  col_fn <- colorRamp2(c(-lim, 0, lim), c("#01665e","#f7f7f7","#8c510a"))
  
  opt_k <- 4
  if (requireNamespace("factoextra",quietly=TRUE)) {
    library(factoextra)
    nb <- fviz_nbclust(as.data.frame(t(pmat)),kmeans,method="silhouette",
                        k.max=min(8,ncol(pmat)-1))
    opt_k <- which.max(nb$data$y)
  }
  
  # legend
  lg_breaks <- round(seq(-lim, lim, length.out=5), 2)
  
  ht <- Heatmap(pmat, name="Log relative\nrisk", col=col_fn,
    top_annotation=ha_top,
    cluster_columns=TRUE, clustering_distance_columns="pearson",
    clustering_method_columns="ward.D2",
    cluster_rows=TRUE, clustering_distance_rows="pearson",
    show_column_names=FALSE,
    row_names_gp=gpar(fontsize=9, fontface="plain"),
    row_names_side="left",
    column_title=sprintf("%s: Per-Patient COO (n=%d)", lbl, ncol(pmat)),
    column_split=opt_k, border=TRUE,
    heatmap_legend_param=list(title="Log relative\nrisk", at=lg_breaks)
  )
  
  pdf(fname, width=max(14, ncol(pmat)*0.15), height=max(8, nrow(pmat)*0.4))
  draw(ht, merge_legend=TRUE, padding=unit(c(2,2,2,10),"mm"))
  dev.off()
  message(sprintf("   %s", fname))
}

plot_heatmap(beta_atlas,  "Atlas",  paste0(OUTDIR,"Fig_Heatmap_Atlas.pdf"))
plot_heatmap(beta_hepcom, "HepCom", paste0(OUTDIR,"Fig_Heatmap_HepCom.pdf"))

# 7C. Sensitivity
message("\n>> Sensitivity")
run_sensitivity <- function(mut_dt, covariates, centroids, hvg_list, cell_types, label) {
  mut_agg <- mut_dt[,.(mutation_count=.N),by=Gene]
  dat <- mut_agg %>% inner_join(covariates,by="Gene") %>%
    inner_join(centroids,by="Gene") %>% filter(Gene %in% hvg_list) %>% as.data.table()
  expr_mat <- as.matrix(dat[,..cell_types]); avg_expr <- rowMeans(expr_mat)
  out <- list()
  for (otype in c("exon_width","gene_width")) {
    bd <- data.frame(y=dat$mutation_count, log_avg=log(avg_expr+1), gc=dat$gc, rept=dat$rept)
    if (otype=="exon_width") {
      bd$log_off <- log(dat$exon_width+1)
      fml <- y~log_expr+log_avg+gc+rept+offset(log_off)
    } else {
      bd$log_off <- log(dat$width); bd$log_ef <- log(dat$exon_frac+1e-6)
      fml <- y~log_expr+log_avg+gc+rept+log_ef+offset(log_off)
    }
    res <- data.frame()
    for (ct in cell_types) {
      md <- bd; md$log_expr <- log(dat[[ct]]+1)
      fit <- tryCatch(glm.nb(fml,data=md,control=glm.control(maxit=100)),error=function(e) NULL)
      if (!is.null(fit)&&"log_expr" %in% rownames(summary(fit)$coefficients))
        res <- rbind(res, data.frame(CellType=ct,
          LogRR=summary(fit)$coefficients["log_expr","Estimate"], Offset=otype))
    }; out[[otype]] <- res
  }
  combined <- do.call(rbind, out)
  r1 <- combined %>% filter(Offset=="exon_width") %>% arrange(LogRR) %>% mutate(rk1=row_number())
  r2 <- combined %>% filter(Offset=="gene_width") %>% arrange(LogRR) %>% mutate(rk2=row_number())
  rr <- inner_join(r1 %>% dplyr::select(CellType,rk1,LogRR1=LogRR),
                   r2 %>% dplyr::select(CellType,rk2,LogRR2=LogRR), by="CellType")
  rho <- cor(rr$rk1, rr$rk2, method="spearman")
  message(sprintf("   [%s] rho: %.3f", label, rho))
  
  p <- ggplot(combined, aes(x=reorder(CellType,LogRR), y=LogRR, fill=Offset)) +
    geom_bar(stat="identity", position=position_dodge(0.8), width=0.7, alpha=0.85) +
    scale_fill_manual(values=c("exon_width"="#C0392B","gene_width"="#2980B9"),
                      labels=c("exon_width"="Exon width (corrected)","gene_width"="Gene width (original)")) +
    geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme_bw() +
    labs(title=sprintf("Sensitivity: %s (rho=%.3f)", label, rho), x="", y="Log relative risk") +
    theme(axis.text.y=element_text(size=9,face="bold"), panel.grid.major.y=element_blank(),
          legend.position="bottom")
  ggsave(paste0(OUTDIR,"Fig_Sensitivity_",label,".pdf"), p, width=10, height=7)
  write.csv(rr, paste0(OUTDIR,"Sensitivity_",label,".csv"), row.names=FALSE)
}
run_sensitivity(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,"Atlas")
run_sensitivity(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,"HepCom")

# ║  PART 8:                                                               ║
message("\n",paste(rep("=",70),collapse=""))
for (res in list(res_atlas,res_hepcom)) {
  if (is.null(res)) next; lbl <- attr(res,"label"); coo <- res$CellType[res$IsCOO]
  if (length(coo)>0) {
    message(sprintf("  COO: %s", paste(coo,collapse=", ")))
    for (c in coo) { r <- res[res$CellType==c,]
      message(sprintf("    %s: RR=%.4f [%.4f-%.4f], FDR=%.2e", c, r$RelativeRisk, r$RR_Lower, r$RR_Upper, r$FDR)) }
  } else {
    t3 <- res %>% arrange(RelativeRisk) %>% head(3)
    message("  COO: none. Top 3:")
    for (i in 1:nrow(t3)) message(sprintf("    %s: RR=%.4f, FDR=%.2e", t3$CellType[i],t3$RelativeRisk[i],t3$FDR[i]))
  }
}

# --- 3f,g alternative: COO with gene-width offset ---

# COO  — Panja et al. 2025 Nature Genetics
# Barplot: Relative Risk (1),

OUTDIR <- "/home/download/csc_article/fig3/COO/"
INDIR  <- "/home/download/COO/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTDIR)

suppressPackageStartupMessages({
  library(MASS); library(data.table); library(dplyr); library(readxl)
  library(readr); library(ggplot2); library(GenomeInfoDb); library(GenomicRanges)
  library(GenomicFeatures); library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(MutationalPatterns)
  library(org.Hs.eg.db); library(rtracklayer); library(ComplexHeatmap)
  library(circlize); library(parallel)
})

# ║  PART 1                                                                    ║
case_raw <- read_excel(paste0(INDIR, "Cases_20260117.xlsx"))
clin <- case_raw %>%
  mutate(
    HBV_Group = case_when(
      is.na(Hepatitis)~"Unknown",
      grepl("HBV",Hepatitis,ignore.case=TRUE)&!grepl("HBV-HCV-",Hepatitis,fixed=TRUE)~"HBV_Pos",
      TRUE~"HBV_Neg"),
    Virus_Status = case_when(
      grepl("HBV",Hepatitis)&grepl("HCV",Hepatitis)~"HBV+HCV",
      grepl("HBV",Hepatitis)&!grepl("HBV-",Hepatitis)~"HBV",
      grepl("HCV",Hepatitis)&!grepl("HCV-",Hepatitis)~"HCV",
      TRUE~"NBNC")
  ) %>% dplyr::select(CaseID, HBV_Group, Virus_Status, Gender, everything())
fwrite(as.data.table(clin), paste0(OUTDIR, "clinical_info.csv"))

mut_raw <- read_excel(paste0(INDIR, "Mutations_20260117.xlsx"))
mut_snp <- mut_raw %>% filter(Type=="SNP", nchar(Ref)==1, nchar(Allele)==1)
mut_gr_df <- as.data.frame(mut_snp) %>%
  dplyr::rename(chromosome=Chr,start=Start,end=End,REF=Ref,ALT=Allele) %>%
  mutate(chromosome=paste0("chr",sub("^chr","",chromosome)))
gr_all <- makeGRangesFromDataFrame(mut_gr_df,keep.extra.columns=TRUE,ignore.strand=TRUE,
  seqnames.field="chromosome",start.field="start",end.field="end")
genome(gr_all)<-"hg19"; seqlevelsStyle(gr_all)<-"UCSC"
gr_all <- keepStandardChromosomes(gr_all,pruning.mode="coarse")

# ║  PART 2: APOBEC                                                            ║
vcfs_by_patient <- split(gr_all, gr_all$CaseID)
mut_mat <- mut_matrix(vcfs_by_patient, "BSgenome.Hsapiens.UCSC.hg19")
cosmic_all <- get_known_signatures(source="COSMIC_v3.2")
fit_res <- fit_to_signatures(mut_mat, cosmic_all)
sig_contrib <- fit_res$contribution

# Signature counts
sig_counts <- as.data.frame(t(sig_contrib))
sig_counts$CaseID <- rownames(sig_counts)
fwrite(sig_counts, paste0(OUTDIR, "Signature_Counts.csv"))

# Motif-based APOBEC (HCC  posterior ,  motif)
tc_res <- type_context(gr_all, "BSgenome.Hsapiens.UCSC.hg19")
is_apobec <- (tc_res$types %in% c("C>T","C>G")) &
  (substr(tc_res$context,1,1)=="T") &
  (substr(tc_res$context,3,3) %in% c("A","T"))
message(sprintf("   APOBEC (TCW motif): %d (%.1f%%)", sum(is_apobec), 100*mean(is_apobec)))

gr_clean <- gr_all[!is_apobec]
mut_clean <- as.data.table(as.data.frame(gr_clean))
if ("seqnames" %in% colnames(mut_clean)) setnames(mut_clean,"seqnames","chromosome")

# ║  PART 3:                                                             ║
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; genes_gr <- genes(txdb)
gene_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, genes_gr)
gc <- letterFrequency(gene_seqs, letters="GC", as.prob=TRUE)[,"G|C"]
exons_reduced <- reduce(exonsBy(txdb, by="gene"))
exon_len <- sum(width(exons_reduced)); exon_len <- exon_len[names(genes_gr)]
exon_len[is.na(exon_len)] <- 0; gene_w <- width(genes_gr)

bw_file <- paste0(INDIR, "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig")
if (!file.exists(bw_file)) download.file(paste0(
  "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/",
  "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig"), destfile=bw_file, mode="wb")
rt_track <- import(bw_file)
ov <- findOverlaps(genes_gr, rt_track)
avg_rt <- tapply(score(rt_track)[subjectHits(ov)], queryHits(ov), mean)
rt_vec <- rep(0, length(genes_gr)); rt_vec[as.numeric(names(avg_rt))] <- avg_rt

covariates <- data.table(entrez_id=names(genes_gr), width=as.numeric(gene_w),
  exon_width=as.numeric(exon_len), gc=as.numeric(gc),
  exon_frac=as.numeric(exon_len/gene_w), rept=rt_vec)
symbols <- mapIds(org.Hs.eg.db, keys=covariates$entrez_id, column="SYMBOL",
                  keytype="ENTREZID", multiVals="first")
covariates$Gene <- symbols
covariates <- covariates[!is.na(Gene)][order(-width)][!duplicated(Gene)]
fwrite(covariates, paste0(OUTDIR, "gene_covariates.csv"))

# ║  PART 4:                                                           ║
run_coo_regression <- function(mut_dt, covariates, centroids, hvg_list, cell_types, label="") {
  mut_agg <- mut_dt[, .(mutation_count=.N), by=Gene]
  dat <- mut_agg %>% inner_join(covariates,by="Gene") %>%
    inner_join(centroids,by="Gene") %>% filter(Gene %in% hvg_list) %>% as.data.table()
  if (nrow(dat)<50) return(NULL)
  n_mut <- sum(dat$mutation_count)
  
  expr_mat <- as.matrix(dat[, ..cell_types]); avg_expr <- rowMeans(expr_mat)
  base <- data.frame(y=dat$mutation_count, log_avg=log(avg_expr+1),
    gc=dat$gc, rept=dat$rept, log_gene_w=log(dat$width))
  
  results <- data.frame()
  for (ct in cell_types) {
    md <- base; md$log_expr <- log(dat[[ct]]+1)
    fit <- tryCatch(glm.nb(y~log_expr+log_avg+gc+rept+offset(log_gene_w),
      data=md, control=glm.control(maxit=100)), error=function(e) NULL)
    if (is.null(fit)) next
    coefs <- summary(fit)$coefficients
    if (!"log_expr" %in% rownames(coefs)) next
    est <- coefs["log_expr","Estimate"]; se <- coefs["log_expr","Std. Error"]
    results <- rbind(results, data.frame(CellType=ct, Estimate=est, StdError=se,
      PValue=coefs["log_expr","Pr(>|z|)"],
      RelativeRisk=exp(est), RR_Lower=exp(est-1.96*se), RR_Upper=exp(est+1.96*se),
      stringsAsFactors=FALSE))
  }
  if (nrow(results)==0) return(NULL)
  results$FDR <- p.adjust(results$PValue, method="BH")
  
  sig_neg <- results %>% filter(FDR<0.1 & Estimate<0)
  coo_set <- character(0)
  if (nrow(sig_neg)>0) {
    top <- sig_neg[which.min(sig_neg$Estimate),]
    coo_set <- top$CellType
    tw <- top$RR_Upper - top$RR_Lower
    for (i in seq_len(nrow(sig_neg))) {
      if (sig_neg$CellType[i]==top$CellType) next
      ov_w <- max(0, min(top$RR_Upper,sig_neg$RR_Upper[i])-max(top$RR_Lower,sig_neg$RR_Lower[i]))
      if (tw>0 && ov_w/tw>=0.6) coo_set <- c(coo_set, sig_neg$CellType[i])
    }
  }
  results$IsCOO <- results$CellType %in% coo_set
  results$Significance <- case_when(
    results$IsCOO~"Inferred COO", results$FDR<0.1&results$Estimate<0~"Sig. negative",
    results$FDR<0.1&results$Estimate>0~"Sig. positive", TRUE~"Not significant")
  
  message(sprintf("   [%s] COO: %s", label,
    ifelse(length(coo_set)>0, paste(coo_set,collapse=", "), "none")))
  attr(results,"label")<-label; attr(results,"n_mutations")<-n_mut; attr(results,"n_genes")<-nrow(dat)
  return(results)
}

# ║  PART 5: Aggregate                                                         ║
message("\n>> PART 5: Aggregate COO")
cent_atlas <- fread(paste0(OUTDIR,"liver_centroids_atlas.csv"))
hvg_atlas <- fread(paste0(OUTDIR,"hvg_atlas.csv"))$Gene
ct_atlas <- setdiff(colnames(cent_atlas),"Gene")

cent_hepcom <- fread(paste0(OUTDIR,"liver_centroids_hepcom.csv"))
hvg_hepcom <- fread(paste0(OUTDIR,"hvg_hepcom.csv"))$Gene
ct_hepcom <- setdiff(colnames(cent_hepcom),"Gene")

message("\n--- Atlas ---")
res_atlas <- run_coo_regression(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,"Atlas")
message("\n--- HepCom ---")
res_hepcom <- run_coo_regression(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,"HepCom")

for (res in list(res_atlas,res_hepcom)) {
  if (!is.null(res)) write.csv(res, paste0(OUTDIR,"COO_Result_",attr(res,"label"),".csv"), row.names=FALSE)
}

# ║  PART 6: Per-Patient                                                       ║
message("\n>> PART 6: Per-Patient")
run_per_patient <- function(mut_dt, covariates, centroids, hvg_list, cell_types, min_snv=80, label="") {
  shared <- Reduce(intersect, list(hvg_list, covariates$Gene, centroids$Gene))
  ref <- merge(covariates[Gene %in% shared], centroids[Gene %in% shared], by="Gene")
  per_gene <- mut_dt[Gene %in% shared, .(N=.N), by=.(CaseID,Gene)]
  valid_pid <- per_gene[,.(Tot=sum(N)),by=CaseID][Tot>=min_snv]$CaseID
  if (length(valid_pid)<5) return(NULL)
  run_one <- function(pid) {
    df <- copy(ref)
    df <- merge(df, per_gene[CaseID==pid,.(Gene,N)], by="Gene", all.x=TRUE); df[is.na(N),N:=0]
    avg_e <- rowMeans(as.matrix(df[,..cell_types]))
    df$log_avg <- log(avg_e+1); df$log_gene_w <- log(df$width)
    beta <- setNames(rep(NA_real_,length(cell_types)),cell_types)
    for (ct in cell_types) {
      df$log_expr <- log(df[[ct]]+1)
      fit <- tryCatch(glm.nb(N~log_expr+log_avg+gc+rept+offset(log_gene_w),
        data=df,control=glm.control(maxit=50)),error=function(e) NULL)
      if (!is.null(fit)&&"log_expr" %in% names(coef(fit))) beta[ct] <- coef(fit)["log_expr"]
    }; return(beta)
  }
  nc <- min(detectCores()-2,8)
  blist <- mclapply(valid_pid, run_one, mc.cores=nc)
  bmat <- do.call(rbind,blist); rownames(bmat) <- valid_pid
  bmat <- bmat[complete.cases(bmat),,drop=FALSE]
}
beta_atlas <- run_per_patient(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,80,"Atlas")
beta_hepcom <- run_per_patient(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,80,"HepCom")
save(beta_atlas,beta_hepcom,sig_counts,file=paste0(OUTDIR,"PerPatient_Betas.RData"))

# ║  PART 7:  (Nature Genetics )                                       ║
clin <- fread(paste0(OUTDIR,"clinical_info.csv"))
sig_counts_dt <- fread(paste0(OUTDIR,"Signature_Counts.csv"))

# 7A. Aggregate Barplot — Relative Risk ( 1 ,  Fig2 )
plot_coo_bar <- function(res, fname) {
  if (is.null(res)) return(invisible(NULL))
  lbl <- attr(res,"label")
  
  cols <- c("Inferred COO"="#C0392B", "Sig. negative"="#E74C3C",
            "Sig. positive"="#2980B9", "Not significant"="#BDC3C7")
  
  res$LogRR    <- res$Estimate
  res$CI_lower <- res$Estimate - 1.96 * res$StdError
  res$CI_upper <- res$Estimate + 1.96 * res$StdError
  
  max_abs <- max(abs(c(res$CI_lower, res$CI_upper)), na.rm=TRUE) * 1.1
  step <- if (max_abs > 0.5) 0.2 else if (max_abs > 0.25) 0.1 else 0.05
  brks <- seq(-ceiling(max_abs/step)*step, ceiling(max_abs/step)*step, by=step)
  
  p <- ggplot(res, aes(x=reorder(CellType, LogRR), y=LogRR, fill=Significance)) +
    geom_bar(stat="identity", width=0.7, alpha=0.9) +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.25, linewidth=0.5, color="grey30") +
    scale_fill_manual(values=cols, drop=FALSE) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.6) +
    coord_flip() +
    scale_y_continuous(breaks=brks, expand=expansion(mult=0.02)) +
    theme_bw(base_size=12) +
    labs(title=sprintf("Cell-of-Origin: %s", lbl),
         subtitle=sprintf("n = %s mutations | %d genes | APOBEC excluded",
                          format(attr(res,"n_mutations"),big.mark=","), attr(res,"n_genes")),
         x="", y="Log relative risk",
         caption="< 0 = TCR protection (potential origin) | > 0 = positive correlation") +
    theme(axis.text.y=element_text(size=10, face="bold", color="black"),
          axis.text.x=element_text(size=9),
          plot.title=element_text(hjust=0.5, face="bold", size=14),
          panel.grid.major.y=element_blank(),
          legend.position="bottom")
  
  ggsave(fname, p, width=9, height=max(6, length(unique(res$CellType))*0.4))
  message(sprintf("   %s", fname))
}

plot_coo_bar(res_atlas,  paste0(OUTDIR,"Fig_COO_Atlas.pdf"))
plot_coo_bar(res_hepcom, paste0(OUTDIR,"Fig_COO_HepCom.pdf"))

# 7B. Per-Patient Heatmap — Nature Fig 3a  ()
plot_heatmap <- function(bmat, lbl, fname) {
  if (is.null(bmat)||nrow(bmat)<5) { message(sprintf("   [%s] skipping heatmap",lbl)); return(invisible(NULL)) }
  
  # beta (= log RR),  0
  pmat <- t(bmat)
  
  adf <- data.frame(CaseID=colnames(pmat), stringsAsFactors=FALSE) %>%
    left_join(clin, by="CaseID") %>%
    left_join(sig_counts_dt, by="CaseID")
  
  tmb_per_patient <- mut_clean[, .N, by=CaseID]
  adf <- adf %>% left_join(as.data.frame(tmb_per_patient) %>% dplyr::rename(TMB=N), by="CaseID")
  adf$TMB[is.na(adf$TMB)] <- 0
  
  # Signature : Aging(SBS1+SBS5), APOBEC(SBS2+SBS13), AA(SBS22), Aflatoxin(SBS24)
  adf$Aging   <- ifelse("SBS1" %in% colnames(adf) & "SBS5" %in% colnames(adf),
                        adf$SBS1 + adf$SBS5, 0)
  adf$APOBEC  <- ifelse("SBS2" %in% colnames(adf) & "SBS13" %in% colnames(adf),
                        adf$SBS2 + adf$SBS13, 0)
  adf$AA      <- ifelse("SBS22" %in% colnames(adf), adf$SBS22, 0)
  adf$Aflatoxin <- ifelse("SBS24" %in% colnames(adf), adf$SBS24, 0)
  
  col_virus  <- c("HBV"="#d95f02","HCV"="#7570b3","HBV+HCV"="#e7298a","NBNC"="#66a61e")
  col_gender <- c("Male"="#1f78b4","Female"="#fb9a99")
  col_hbv    <- c("HBV_Pos"="#d95f02","HBV_Neg"="#66a61e","Unknown"="grey85")
  
  col_sig <- c("AA"="#984ea3", "APOBEC"="#e41a1c", "Aflatoxin"="#ffff33", "Aging"="#bababa")
  
  ha_top <- HeatmapAnnotation(
    "TMB" = anno_barplot(adf$TMB, gp=gpar(fill="black",col=NA),
                         height=unit(1.2,"cm"), border=FALSE,
                         axis_param=list(side="right",gp=gpar(fontsize=7))),
    "Signatures" = anno_barplot(
      as.matrix(adf[, c("AA","Aflatoxin","APOBEC","Aging")]),
      gp=gpar(fill=col_sig, col=NA), height=unit(1.5,"cm"), border=FALSE),
    "Virus"  = adf$Virus_Status,
    "HBV"    = adf$HBV_Group,
    "Gender" = adf$Gender,
    col = list("Virus"=col_virus, "HBV"=col_hbv, "Gender"=col_gender),
    na_col = "grey90",
    gap = unit(c(1,2,1,1), "mm"),
    annotation_name_side = "left",
    annotation_label = c("TMB","Signatures","Virus","HBV","Gender")
  )
  
  # log(RR)  0  (=teal=COO, =brown)
  lim <- max(abs(quantile(pmat, c(0.02,0.98), na.rm=TRUE)))
  col_fn <- colorRamp2(c(-lim, 0, lim), c("#01665e","#f7f7f7","#8c510a"))
  
  opt_k <- 4
  if (requireNamespace("factoextra",quietly=TRUE)) {
    library(factoextra)
    nb <- fviz_nbclust(as.data.frame(t(pmat)),kmeans,method="silhouette",
                        k.max=min(8,ncol(pmat)-1))
    opt_k <- which.max(nb$data$y)
  }
  
  # legend
  lg_breaks <- round(seq(-lim, lim, length.out=5), 2)
  
  ht <- Heatmap(pmat, name="Log relative\nrisk", col=col_fn,
    top_annotation=ha_top,
    cluster_columns=TRUE, clustering_distance_columns="pearson",
    clustering_method_columns="ward.D2",
    cluster_rows=TRUE, clustering_distance_rows="pearson",
    show_column_names=FALSE,
    row_names_gp=gpar(fontsize=9, fontface="plain"),
    row_names_side="left",
    column_title=sprintf("%s: Per-Patient COO (n=%d)", lbl, ncol(pmat)),
    column_split=opt_k, border=TRUE,
    heatmap_legend_param=list(title="Log relative\nrisk", at=lg_breaks)
  )
  
  pdf(fname, width=max(14, ncol(pmat)*0.15), height=max(8, nrow(pmat)*0.4))
  draw(ht, merge_legend=TRUE, padding=unit(c(2,2,2,10),"mm"))
  dev.off()
  message(sprintf("   %s", fname))
}

plot_heatmap(beta_atlas,  "Atlas",  paste0(OUTDIR,"Fig_Heatmap_Atlas.pdf"))
plot_heatmap(beta_hepcom, "HepCom", paste0(OUTDIR,"Fig_Heatmap_HepCom.pdf"))

# 7C. Sensitivity
message("\n>> Sensitivity")
run_sensitivity <- function(mut_dt, covariates, centroids, hvg_list, cell_types, label) {
  mut_agg <- mut_dt[,.(mutation_count=.N),by=Gene]
  dat <- mut_agg %>% inner_join(covariates,by="Gene") %>%
    inner_join(centroids,by="Gene") %>% filter(Gene %in% hvg_list) %>% as.data.table()
  expr_mat <- as.matrix(dat[,..cell_types]); avg_expr <- rowMeans(expr_mat)
  out <- list()
  for (otype in c("gene_width","exon_width")) {
    bd <- data.frame(y=dat$mutation_count, log_avg=log(avg_expr+1), gc=dat$gc, rept=dat$rept)
    if (otype=="gene_width") {
      bd$log_off <- log(dat$width)
      fml <- y~log_expr+log_avg+gc+rept+offset(log_off)
    } else {
      bd$log_off <- log(dat$exon_width+1); bd$log_ef <- log(dat$exon_frac+1e-6)
      fml <- y~log_expr+log_avg+gc+rept+offset(log_off)
    }
    res <- data.frame()
    for (ct in cell_types) {
      md <- bd; md$log_expr <- log(dat[[ct]]+1)
      fit <- tryCatch(glm.nb(fml,data=md,control=glm.control(maxit=100)),error=function(e) NULL)
      if (!is.null(fit)&&"log_expr" %in% rownames(summary(fit)$coefficients))
        res <- rbind(res, data.frame(CellType=ct,
          LogRR=summary(fit)$coefficients["log_expr","Estimate"], Offset=otype))
    }; out[[otype]] <- res
  }
  combined <- do.call(rbind, out)
  r1 <- combined %>% filter(Offset=="exon_width") %>% arrange(LogRR) %>% mutate(rk1=row_number())
  r2 <- combined %>% filter(Offset=="gene_width") %>% arrange(LogRR) %>% mutate(rk2=row_number())
  rr <- inner_join(r1 %>% dplyr::select(CellType,rk1,LogRR1=LogRR),
                   r2 %>% dplyr::select(CellType,rk2,LogRR2=LogRR), by="CellType")
  rho <- cor(rr$rk1, rr$rk2, method="spearman")
  message(sprintf("   [%s] rho: %.3f", label, rho))
  
  p <- ggplot(combined, aes(x=reorder(CellType,LogRR), y=LogRR, fill=Offset)) +
    geom_bar(stat="identity", position=position_dodge(0.8), width=0.7, alpha=0.85) +
    scale_fill_manual(values=c("gene_width"="#2980B9","exon_width"="#C0392B"),
                      labels=c("gene_width"="Gene width (primary)","exon_width"="Exon width (sensitivity)")) +
    geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme_bw() +
    labs(title=sprintf("Sensitivity: %s (rho=%.3f)", label, rho), x="", y="Log relative risk") +
    theme(axis.text.y=element_text(size=9,face="bold"), panel.grid.major.y=element_blank(),
          legend.position="bottom")
  ggsave(paste0(OUTDIR,"Fig_Sensitivity_",label,".pdf"), p, width=10, height=7)
  write.csv(rr, paste0(OUTDIR,"Sensitivity_",label,".csv"), row.names=FALSE)
}
run_sensitivity(mut_clean,covariates,cent_atlas,hvg_atlas,ct_atlas,"Atlas")
run_sensitivity(mut_clean,covariates,cent_hepcom,hvg_hepcom,ct_hepcom,"HepCom")

# ║  PART 8:                                                               ║
message("\n",paste(rep("=",70),collapse=""))
for (res in list(res_atlas,res_hepcom)) {
  if (is.null(res)) next; lbl <- attr(res,"label"); coo <- res$CellType[res$IsCOO]
  if (length(coo)>0) {
    message(sprintf("  COO: %s", paste(coo,collapse=", ")))
    for (c in coo) { r <- res[res$CellType==c,]
      message(sprintf("    %s: RR=%.4f [%.4f-%.4f], FDR=%.2e", c, r$RelativeRisk, r$RR_Lower, r$RR_Upper, r$FDR)) }
  } else {
    t3 <- res %>% arrange(RelativeRisk) %>% head(3)
    message("  COO: none. Top 3:")
    for (i in 1:nrow(t3)) message(sprintf("    %s: RR=%.4f, FDR=%.2e", t3$CellType[i],t3$RelativeRisk[i],t3$FDR[i]))
  }
}
