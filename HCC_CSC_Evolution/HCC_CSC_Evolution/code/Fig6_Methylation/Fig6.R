#!/usr/bin/env Rscript
# Figure 6 — Dynamic methylation experiments reveal stable epigenetic memory

# --- 6b: Epigenetic entropy ---

#!/usr/bin/env Rscript
# (Epigenetic Entropy)
# 1.  cgMerged.935K.txt  beta bedgraph
# 2.  /home/root2/methylation/output/
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(data.table)

# ======================== 0.  ==========================================
input_dir   <- "/home/root2/methylation/input_mat"
output_dir  <- "/home/root2/methylation/output"
entropy_dir <- file.path(output_dir, "Epigenetic_Entropy")
dir.create(entropy_dir, showWarnings = FALSE, recursive = TRUE)

MIN_COV <- 5

# ======================== 1.  Shannon  ==============================
calc_cpg_entropy <- function(beta) {
  beta <- pmax(pmin(beta, 1 - 1e-10), 1e-10)
  entropy <- -(beta * log2(beta) + (1 - beta) * log2(1 - beta))
  return(entropy)
}

calc_sample_entropy <- function(beta_vec) {
  valid <- !is.na(beta_vec)
  mean(calc_cpg_entropy(beta_vec[valid]))
}

# ======================== 2.  cgMerged  beta ==================
cg_files  <- list.files(input_dir, pattern = "\\.cgMerged\\.935K\\.txt$", full.names = TRUE)
bed_names <- gsub("\\.cgMerged\\.935K\\.txt$", "", basename(cg_files))

beta_list <- list()
for (i in seq_along(cg_files)) {
  # cgMerged: V1=chr, V2=start, V3=strand, ..., V7=numCs, V8=numTs
  tmp <- fread(cg_files[i], select = c(1, 2, 7, 8),
               col.names = c("chr", "start", "numCs", "numTs"))
  tmp[, coverage := numCs + numTs]
  tmp <- tmp[coverage >= MIN_COV]
  tmp[, beta := numCs / coverage]
  tmp[, site_id := paste0(chr, ":", start)]
  beta_list[[bed_names[i]]] <- as.data.frame(tmp[, .(site_id, beta)])
}

# ======================== 3.  =================================
global_entropy <- data.frame(
  sample = bed_names,
  mean_entropy = sapply(beta_list, function(x) calc_sample_entropy(x$beta)),
  median_entropy = sapply(beta_list, function(x) median(calc_cpg_entropy(x$beta), na.rm = TRUE)),
  n_sites = sapply(beta_list, function(x) sum(!is.na(x$beta))),
  stringsAsFactors = FALSE
)

global_entropy$group <- NA
global_entropy$group[grepl("^N_[1-4]_B$", global_entropy$sample)] <- "FN"
global_entropy$group[grepl("^P_[1-4]_B$", global_entropy$sample)] <- "FP"
global_entropy$group[grepl("^P_N_",        global_entropy$sample)] <- "P_N"
global_entropy$group[grepl("^P_P_",        global_entropy$sample)] <- "P_P"
global_entropy$group[grepl("^N_[1-2]_N$",  global_entropy$sample)] <- "N_N"
global_entropy$group[grepl("^N_[1-2]_P$",  global_entropy$sample)] <- "N_P"

write.csv(global_entropy, file.path(entropy_dir, "Global_Entropy_per_sample.csv"), row.names = FALSE)
print(global_entropy)

# ======================== 4.  ======================================
comparison_list <- list(
  "FN_vs_FP" = list(
    group0 = "FN", group1 = "FP",
    label0 = "FN", label1 = "FP",
    color0 = "#2166AC", color1 = "#B2182B"
  ),
  "P_N_vs_P_P" = list(
    group0 = "P_N", group1 = "P_P",
    label0 = "P_N", label1 = "P_P",
    color0 = "#4393C3", color1 = "#D6604D"
  ),
  "N_N_vs_N_P" = list(
    group0 = "N_N", group1 = "N_P",
    label0 = "N_N", label1 = "N_P",
    color0 = "#92C5DE", color1 = "#F4A582"
  )
)

# ======================== 5.  ===================================

save_plot <- function(p, path_no_ext, w, h) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(path_no_ext, ".png"), p, width = w, height = h, dpi = 300)
}

for (comp_name in names(comparison_list)) {
  comp <- comparison_list[[comp_name]]
  comp_dir <- file.path(entropy_dir, comp_name)
  dir.create(comp_dir, showWarnings = FALSE)
  
  message("\n--- ", comp_name, " ---")
  samples_g0 <- global_entropy$sample[global_entropy$group == comp$group0]
  samples_g1 <- global_entropy$sample[global_entropy$group == comp$group1]
  all_samples <- c(samples_g0, samples_g1)
  
  message("  Group0 (", comp$label0, "): ", paste(samples_g0, collapse = ", "))
  message("  Group1 (", comp$label1, "): ", paste(samples_g1, collapse = ", "))
  
  common_sites <- Reduce(intersect, lapply(all_samples, function(s) beta_list[[s]]$site_id))
  
  beta_matrix <- matrix(NA, nrow = length(common_sites), ncol = length(all_samples),
                        dimnames = list(common_sites, all_samples))
  for (s in all_samples) {
    idx <- match(common_sites, beta_list[[s]]$site_id)
    beta_matrix[, s] <- beta_list[[s]]$beta[idx]
  }
  
  entropy_matrix <- calc_cpg_entropy(beta_matrix)
  sample_entropy <- colMeans(entropy_matrix, na.rm = TRUE)
  
  entropy_df <- data.frame(
    sample = all_samples, entropy = sample_entropy,
    group = factor(ifelse(all_samples %in% samples_g0, comp$label0, comp$label1),
                   levels = c(comp$label0, comp$label1))
  )
  
  p1 <- ggplot(entropy_df, aes(x = group, y = entropy, fill = group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
    scale_fill_manual(values = c(comp$color0, comp$color1)) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x.npc = "center", size = 5) +
    labs(title = paste0("Epigenetic Entropy: ", comp_name),
         y = "Mean Shannon Entropy", x = "") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p1, file.path(comp_dir, "01_Global_Entropy_Boxplot"), 5, 6)
  
  mean_beta_g0 <- rowMeans(beta_matrix[, samples_g0, drop = FALSE], na.rm = TRUE)
  mean_beta_g1 <- rowMeans(beta_matrix[, samples_g1, drop = FALSE], na.rm = TRUE)
  entropy_g0 <- calc_cpg_entropy(mean_beta_g0)
  entropy_g1 <- calc_cpg_entropy(mean_beta_g1)
  
  density_df <- data.frame(
    entropy = c(entropy_g0, entropy_g1),
    group = factor(rep(c(comp$label0, comp$label1), each = length(common_sites)),
                   levels = c(comp$label0, comp$label1))
  )
  
  p2 <- ggplot(density_df, aes(x = entropy, fill = group, color = group)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    scale_fill_manual(values = c(comp$color0, comp$color1)) +
    scale_color_manual(values = c(comp$color0, comp$color1)) +
    labs(title = paste0("CpG-level Entropy Distribution: ", comp_name),
         x = "Shannon Entropy per CpG", y = "Density") +
    theme_bw(base_size = 14) +
    theme(legend.position = c(0.8, 0.85),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p2, file.path(comp_dir, "02_CpG_Entropy_Density"), 8, 5)
  
  # --- 5.3 Epigenetic Erosion ---
  delta_entropy <- entropy_g1 - entropy_g0
  erosion_df <- data.frame(mean_beta_g0 = mean_beta_g0, delta_entropy = delta_entropy)
  set.seed(42)
  erosion_sub <- erosion_df[sample(nrow(erosion_df), min(50000, nrow(erosion_df))), ]
  
  p4 <- ggplot(erosion_sub, aes(x = mean_beta_g0, y = delta_entropy)) +
    geom_point(alpha = 0.05, size = 0.3, color = "grey30") +
    geom_smooth(method = "loess", color = "#B2182B", linewidth = 1.2, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = paste0("Epigenetic Erosion: ", comp_name),
         subtitle = paste0("\u0394 Entropy = ", comp$label1, " \u2212 ", comp$label0),
         x = paste0("Mean \u03B2 in ", comp$label0),
         y = "\u0394 Shannon Entropy") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p4, file.path(comp_dir, "04_Entropy_Erosion_Scatter"), 8, 6)
}

# ======================== 6.  ======================================
plot_df <- global_entropy[!is.na(global_entropy$group), ]
plot_df$group <- factor(plot_df$group,
                        levels = c("FN", "FP", "P_N", "P_P", "N_N", "N_P"))

group_colors <- c(
  "FN"  = "#2166AC", "FP"  = "#B2182B",
  "P_N" = "#4393C3", "P_P" = "#D6604D",
  "N_N" = "#92C5DE", "N_P" = "#F4A582"
)

p_all <- ggplot(plot_df, aes(x = group, y = mean_entropy, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(c("FN", "FP"), c("P_N", "P_P"), c("N_N", "N_P")),
    method = "wilcox.test", label = "p.format", size = 4
  ) +
  labs(title = "Epigenetic Entropy Across All Groups",
       y = "Mean Shannon Entropy", x = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

save_plot(p_all, file.path(entropy_dir, "All_Groups_Entropy_Comparison"), 10, 7)

# ======================== 7.  ====================================
cat("\n============================================================\n")
cat("============================================================\n\n")

group_summary <- plot_df %>%
  group_by(group) %>%
  summarise(n = n(),
            mean_entropy = round(mean(mean_entropy), 6),
            sd_entropy = round(sd(mean_entropy), 6),
            .groups = "drop")
print(as.data.frame(group_summary))

for (comp_name in names(comparison_list)) {
  comp <- comparison_list[[comp_name]]
  g0_vals <- plot_df$mean_entropy[plot_df$group == comp$group0]
  g1_vals <- plot_df$mean_entropy[plot_df$group == comp$group1]
  if (length(g0_vals) >= 2 & length(g1_vals) >= 2) {
    wt <- wilcox.test(g0_vals, g1_vals)
    cat(sprintf("  %s: %s (%.6f) vs %s (%.6f), p = %.4f\n",
                comp_name, comp$label0, mean(g0_vals),
                comp$label1, mean(g1_vals), wt$p.value))
  } else {
  }
}

# --- 6c: Global methylation levels ---

#!/usr/bin/env Rscript
# 4.  PDF + PNG

library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ======================== 1.  ========================
input_dir  <- "/home/root2/methylation/input_mat"
output_dir <- "/home/root2/methylation/output/global/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_COV <- 5

# ======================== 2.  ========================

cg_files <- list.files(input_dir, pattern = "\\.cgMerged\\.935K\\.txt$", full.names = TRUE)

if (length(cg_files) == 0) {
  input_dir <- "/home/root2/methylation/input_bed/"
  cg_files  <- list.files(input_dir, pattern = "\\.CG_5x\\.bedgraph\\.gz$", full.names = TRUE)
  use_bedgraph <- TRUE
} else {
  use_bedgraph <- FALSE
}

if (length(cg_files) == 0) stop("No input files found, check path!")

results <- list()

for (f in cg_files) {
  if (use_bedgraph) {
    sample_id <- gsub("\\.CG_5x\\.bedgraph\\.gz$", "", basename(f))
    dt <- fread(cmd = paste("zcat", f), header = FALSE, select = 4, col.names = "beta")
    mean_meth   <- mean(dt$beta, na.rm = TRUE)
    median_meth <- median(dt$beta, na.rm = TRUE)
    cpg_count   <- nrow(dt)
  } else {
    sample_id <- gsub("\\.cgMerged\\.935K\\.txt$", "", basename(f))
    dt <- fread(f, select = c(7, 8), col.names = c("numCs", "numTs"))
    dt[, coverage := numCs + numTs]
    dt <- dt[coverage >= MIN_COV]
    dt[, beta := numCs / coverage]
    mean_meth   <- mean(dt$beta, na.rm = TRUE)
    median_meth <- median(dt$beta, na.rm = TRUE)
    cpg_count   <- nrow(dt)
  }

  results[[sample_id]] <- data.frame(
    Sample = sample_id,
    Mean_Methylation = mean_meth,
    Median_Methylation = median_meth,
    CpG_Count = cpg_count,
    stringsAsFactors = FALSE
  )
}

df_summary <- bind_rows(results)

# ======================== 3.  ========================
df_summary <- df_summary %>%
  mutate(group = case_when(
    grepl("^N_[1-4]_B$", Sample) ~ "FN",
    grepl("^P_[1-4]_B$", Sample) ~ "FP",
    grepl("^P_N_",        Sample) ~ "P_N",
    grepl("^P_P_",        Sample) ~ "P_P",
    grepl("^N_[1-2]_N$",  Sample) ~ "N_N",
    grepl("^N_[1-2]_P$",  Sample) ~ "N_P",
    TRUE ~ "Other"
  ))

df_summary$group <- factor(df_summary$group,
                           levels = c("FN", "FP", "P_N", "P_P", "N_N", "N_P"))

print(df_summary)
write.csv(df_summary, file.path(output_dir, "Global_Methylation_Summary.csv"), row.names = FALSE)

# ======================== 4.  ========================

group_colors <- c(
  "FN"  = "#2166AC",
  "FP"  = "#B2182B",
  "P_N" = "#4393C3",
  "P_P" = "#D6604D",
  "N_N" = "#92C5DE",
  "N_P" = "#F4A582"
)

my_comparisons <- list(c("FN", "FP"), c("P_N", "P_P"), c("N_N", "N_P"))

p <- ggplot(df_summary, aes(x = group, y = Mean_Methylation, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", label = "p.format", size = 4
  ) +
  labs(title = "Global DNA Methylation Levels Across All Groups",
       y = "Mean Methylation Level (\u03B2)", x = "") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# PDF + PNG
ggsave(file.path(output_dir, "Global_Methylation_Boxplot.pdf"),
       p, width = 10, height = 7, device = cairo_pdf)
#ggsave(file.path(output_dir, "Global_Methylation_Boxplot.png"),
#       p, width = 10, height = 7, dpi = 300)

# ======================== 5.  ========================
group_stats <- df_summary %>%
  group_by(group) %>%
  summarise(n = n(),
            mean_beta = round(mean(Mean_Methylation), 4),
            sd_beta   = round(sd(Mean_Methylation), 4),
            .groups = "drop")
print(as.data.frame(group_stats))

for (comp in my_comparisons) {
  g0 <- df_summary$Mean_Methylation[df_summary$group == comp[1]]
  g1 <- df_summary$Mean_Methylation[df_summary$group == comp[2]]
  if (length(g0) >= 2 & length(g1) >= 2) {
    wt <- wilcox.test(g0, g1)
    cat(sprintf("  %s vs %s: p = %.4f\n", comp[1], comp[2], wt$p.value))
  }
}

# --- 6d: DMC enrichment in ChromHMM chromatin states ---

# ║  PART 2: EPIC × ChromHMM 18-state (root2)        ║
if (RUN_PART %in% c("2", "all")) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("  PART 2: EPIC Methylation Entropy × ChromHMM 18-state\n")
  cat(strrep("=", 60), "\n")
  
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(data.table)
    library(patchwork)
  })
  
  EPIC_DIR   <- "/home/root2/methylation/output/dma/"
  CHROM_FILE <- file.path(EPIC_DIR, "E008_18_core_K27ac_mnemonics.bed.gz")
  CHAIN_FILE <- file.path(EPIC_DIR, "hg19ToHg38.over.chain")
  OUT_DIR    <- "/home/root2/methylation/output/entropy_chromhmm/"
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(CHAIN_FILE)) {
    message(">>> Downloading liftOver chain...")
    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                  paste0(CHAIN_FILE, ".gz"))
    system(paste("gunzip", paste0(CHAIN_FILE, ".gz")))
  }
  
  # --- ChromHMM ---
  message(">>> Loading ChromHMM 18-state (E008 H9 ESC)...")
  chain <- import.chain(CHAIN_FILE)
  chrom_df <- fread(cmd = paste("zcat", CHROM_FILE), header = FALSE,
                    col.names = c("chr", "start", "end", "state"))
  chrom_gr_hg19 <- GRanges(chrom_df$chr, IRanges(chrom_df$start + 1, chrom_df$end),
                           state = chrom_df$state)
  chrom_gr <- unlist(liftOver(chrom_gr_hg19, chain))
  
  STATE_ORDER <- c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD",
                   "5_Tx","6_TxWk","7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk",
                   "12_ZNF/Rpts","13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies")
  
  STATE_CLASS <- c(
    "1_TssA"="Promoter","2_TssFlnk"="Promoter","3_TssFlnkU"="Promoter","4_TssFlnkD"="Promoter",
    "5_Tx"="Transcription","6_TxWk"="Transcription",
    "7_EnhG1"="Enhancer","8_EnhG2"="Enhancer","9_EnhA1"="Enhancer","10_EnhA2"="Enhancer","11_EnhWk"="Enhancer",
    "12_ZNF/Rpts"="Heterochromatin","13_Het"="Heterochromatin",
    "14_TssBiv"="Bivalent/Polycomb","15_EnhBiv"="Bivalent/Polycomb",
    "16_ReprPC"="Bivalent/Polycomb","17_ReprPCWk"="Bivalent/Polycomb",
    "18_Quies"="Quiescent"
  )
  
  # --- EPIC ---
  epic_files <- list.files(EPIC_DIR, pattern = "All_tested_sites.csv",
                           recursive = TRUE, full.names = TRUE)
  message(sprintf(">>> Found %d EPIC comparison files", length(epic_files)))
  
  all_entropy_df <- data.frame()
  
  for (f in epic_files) {
    comp_name <- basename(dirname(f))
    message(sprintf("   -> %s", comp_name))
    
    df <- fread(f)
    if (!"chr" %in% names(df) || !"start" %in% names(df)) {
      message("     SKIP: missing chr/start columns")
      next
    }
    
    # beta (treatmentcontrolbeta)
    # methylKitAll_tested_sites.csv: chr, start, end, strand, pvalue, qvalue, meth.diff
    # beta。meth.diff = treatment - control
    # meth.diff : treatment_beta ≈ 50 + meth.diff/2 ()
    
    beta_cols_t <- grep("^numCs\\.", names(df), value = TRUE)
    beta_cols_n <- grep("^numTs\\.", names(df), value = TRUE)
    
    if (length(beta_cols_t) >= 2 && length(beta_cols_n) >= 2) {
      # numCs/numTs -> beta
      # treatment (.1), control (.2)
      cov_t <- df[[beta_cols_t[1]]] + df[[beta_cols_n[1]]]
      cov_c <- df[[beta_cols_t[2]]] + df[[beta_cols_n[2]]]
      beta_t <- df[[beta_cols_t[1]]] / pmax(cov_t, 1)
      beta_c <- df[[beta_cols_t[2]]] / pmax(cov_c, 1)
    } else {
      # control ≈ 0.5, treatment ≈ 0.5 + meth.diff/100
      message("     Using meth.diff approximation for beta")
      beta_c <- rep(0.5, nrow(df))
      beta_t <- pmin(pmax(0.5 + df$meth.diff / 100, 0.01), 0.99)
    }
    
    # CpGShannon: H = -(β·log2(β) + (1-β)·log2(1-β))
    safe_entropy <- function(b) {
      b <- pmin(pmax(b, 0.001), 0.999)
      -(b * log2(b) + (1 - b) * log2(1 - b))
    }
    
    entropy_t <- safe_entropy(beta_t)
    entropy_c <- safe_entropy(beta_c)
    
    # GRanges
    gr <- GRanges(df$chr, IRanges(df$start, df$end))
    
    # ChromHMM
    hits <- findOverlaps(gr, chrom_gr)
    state_annot <- rep(NA_character_, length(gr))
    state_annot[queryHits(hits)] <- chrom_gr$state[subjectHits(hits)]
    dup_q <- duplicated(queryHits(hits))
    
    tmp <- data.frame(
      comp = comp_name,
      state = state_annot,
      entropy_treatment = entropy_t,
      entropy_control = entropy_c,
      meth_diff = df$meth.diff,
      stringsAsFactors = FALSE
    )
    tmp <- tmp[!is.na(tmp$state), ]
    all_entropy_df <- rbind(all_entropy_df, tmp)
  }
  
  if (nrow(all_entropy_df) == 0) {
    stop("No data after ChromHMM annotation! Check EPIC file format.")
  }
  
  message(sprintf(">>> Total annotated CpGs: %d", nrow(all_entropy_df)))
  all_entropy_df$state <- factor(all_entropy_df$state, levels = STATE_ORDER)
  all_entropy_df$state_class <- STATE_CLASS[as.character(all_entropy_df$state)]
  
  fwrite(all_entropy_df, file.path(OUT_DIR, "EPIC_Entropy_by_ChromHMM.csv"))
  
  message(">>> Generating main figure...")
  
  summary_df <- all_entropy_df %>%
    group_by(comp, state, state_class) %>%
    summarise(
      mean_entropy_t = mean(entropy_treatment, na.rm = TRUE),
      mean_entropy_c = mean(entropy_control, na.rm = TRUE),
      median_entropy_t = median(entropy_treatment, na.rm = TRUE),
      median_entropy_c = median(entropy_control, na.rm = TRUE),
      n_cpg = n(),
      .groups = "drop"
    )
  
  # --- Panel A: Boxplot  by ChromHMM state (treatment) ---
  # Bivalent/Polycomb
  state_colors <- c(
    "1_TssA"="#FF4500","2_TssFlnk"="#FF6347","3_TssFlnkU"="#FF6347","4_TssFlnkD"="#FF6347",
    "5_Tx"="#00B000","6_TxWk"="#00B000",
    "7_EnhG1"="#FFC34D","8_EnhG2"="#FFC34D","9_EnhA1"="#FFC34D","10_EnhA2"="#FFC34D","11_EnhWk"="#FFC34D",
    "12_ZNF/Rpts"="#8A91D0","13_Het"="#8A91D0",
    "14_TssBiv"="#B71C1C","15_EnhBiv"="#B71C1C","16_ReprPC"="#4A148C","17_ReprPCWk"="#4A148C",
    "18_Quies"="#C0C0C0"
  )
  
  p_a <- ggplot(all_entropy_df, aes(x = state, y = entropy_treatment, fill = state_class)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.3) +
    facet_wrap(~comp, ncol = 2) +
    scale_fill_manual(values = c(
      "Promoter" = "#FF6347",
      "Transcription" = "#00B000",
      "Enhancer" = "#FFC34D",
      "Heterochromatin" = "#8A91D0",
      "Bivalent/Polycomb" = "#C62828",
      "Quiescent" = "#C0C0C0"
    ), name = "Chromatin Class") +
    labs(title = "A. CpG Methylation Entropy by ChromHMM State",
         subtitle = "Low entropy = ordered methylation; High entropy = disordered",
         y = "Methylation Entropy H(β)", x = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "top"
    )
  
  # --- Panel B:  - Bivalent vs  (barplot) ---
  all_entropy_df$is_bivalent <- ifelse(
    all_entropy_df$state %in% c("14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk"),
    "Bivalent/Polycomb", "Other states")
  
  biv_summary <- all_entropy_df %>%
    group_by(comp, is_bivalent) %>%
    summarise(
      mean_H = mean(entropy_treatment, na.rm = TRUE),
      se_H = sd(entropy_treatment, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  # Wilcoxon
  biv_tests <- all_entropy_df %>%
    group_by(comp) %>%
    summarise(
      p_val = wilcox.test(
        entropy_treatment[is_bivalent == "Bivalent/Polycomb"],
        entropy_treatment[is_bivalent == "Other states"],
        alternative = "less"  # Bivalent < Other
      )$p.value,
      .groups = "drop"
    )
  biv_tests$label <- ifelse(biv_tests$p_val < 0.001, "***",
                            ifelse(biv_tests$p_val < 0.01, "**",
                                   ifelse(biv_tests$p_val < 0.05, "*", "ns")))
  
  p_b <- ggplot(biv_summary, aes(x = comp, y = mean_H, fill = is_bivalent)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
    geom_errorbar(aes(ymin = mean_H - se_H, ymax = mean_H + se_H),
                  position = position_dodge(0.7), width = 0.2) +
    scale_fill_manual(values = c("Bivalent/Polycomb" = "#C62828", "Other states" = "#90A4AE"),
                      name = NULL) +
    labs(title = "B. Bivalent/Polycomb Regions Show Lower Methylation Entropy",
         subtitle = paste(biv_tests$comp, biv_tests$label, collapse = "  |  "),
         y = "Mean Methylation Entropy H(β)", x = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
  
  # --- Panel C: 4Bivalent/Polycomb ---
  biv_detail <- all_entropy_df %>%
    filter(state %in% c("14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk"))
  
  p_c <- ggplot(biv_detail, aes(x = state, y = entropy_treatment, fill = state)) +
    geom_violin(alpha = 0.7, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
    facet_wrap(~comp, ncol = 2) +
    scale_fill_manual(values = c(
      "14_TssBiv" = "#D32F2F", "15_EnhBiv" = "#F57C00",
      "16_ReprPC" = "#7B1FA2", "17_ReprPCWk" = "#9C27B0"
    )) +
    labs(title = "C. Methylation Entropy in Bivalent/Polycomb Sub-states",
         subtitle = "TssBiv and ReprPC show lowest entropy → most ordered methylation",
         y = "Methylation Entropy H(β)", x = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
  
  p_final <- p_a / (p_b | p_c) + plot_layout(heights = c(1.2, 1))
  
  out_pdf <- file.path(OUT_DIR, "Fig_MethylationEntropy_ChromHMM.pdf")
  ggsave(out_pdf, p_final, width = 14, height = 16, device = cairo_pdf)
  message(">>> PDF: ", out_pdf)
  
  out_png <- file.path(OUT_DIR, "Fig_MethylationEntropy_ChromHMM.png")
  ggsave(out_png, p_final, width = 14, height = 16, dpi = 300)
  message(">>> PNG: ", out_png)
  
  cat("\n--- Methylation Entropy Summary ---\n")
  state_means <- all_entropy_df %>%
    group_by(state, state_class) %>%
    summarise(mean_H = mean(entropy_treatment, na.rm = TRUE),
              median_H = median(entropy_treatment, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    arrange(mean_H)
  print(as.data.frame(state_means))
  
  cat("\n--- Bivalent vs Other: Wilcoxon Tests ---\n")
  print(as.data.frame(biv_tests))
  
  cat("\n--- Key Result ---\n")
  biv_mean <- mean(all_entropy_df$entropy_treatment[all_entropy_df$is_bivalent == "Bivalent/Polycomb"], na.rm = TRUE)
  oth_mean <- mean(all_entropy_df$entropy_treatment[all_entropy_df$is_bivalent != "Bivalent/Polycomb"], na.rm = TRUE)
  cat(sprintf("  Bivalent/Polycomb mean H = %.4f\n", biv_mean))
  cat(sprintf("  Other states mean H      = %.4f\n", oth_mean))
  cat(sprintf("  Difference               = %.4f\n", biv_mean - oth_mean))
  
  message("\n>>> Part 2 complete. Output: ", OUT_DIR)
  
} # end Part 2

# --- 6e,g: Per-site methylation correlation + heterogeneity ---

#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

input_dir   <- "/home/root2/methylation/input_mat"
output_dir  <- "/home/root2/methylation/output"
hetero_dir  <- file.path(output_dir, "Heterogeneity")
dir.create(hetero_dir, showWarnings = FALSE, recursive = TRUE)

MIN_COV <- 5

save_plot <- function(p, path_no_ext, w, h) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(path_no_ext, ".png"), p, width = w, height = h, dpi = 300)
}

cg_files  <- list.files(input_dir, pattern = "\\.cgMerged\\.935K\\.txt$", full.names = TRUE)
bed_names <- gsub("\\.cgMerged\\.935K\\.txt$", "", basename(cg_files))

beta_list <- list()
for (i in seq_along(cg_files)) {
  tmp <- fread(cg_files[i], select = c(1, 2, 7, 8),
               col.names = c("chr", "start", "numCs", "numTs"))
  tmp[, coverage := numCs + numTs]
  tmp <- tmp[coverage >= MIN_COV]
  tmp[, beta := numCs / coverage]
  tmp[, site_id := paste0(chr, ":", start)]
  beta_list[[bed_names[i]]] <- tmp
}

assign_group <- function(s) {
  case_when(
    grepl("^N_[1-4]_B$", s) ~ "FN",
    grepl("^P_[1-4]_B$", s) ~ "FP",
    grepl("^P_N_",        s) ~ "P_N",
    grepl("^P_P_",        s) ~ "P_P",
    grepl("^N_[1-2]_N$",  s) ~ "N_N",
    grepl("^N_[1-2]_P$",  s) ~ "N_P",
    TRUE ~ NA_character_
  )
}

sample_info <- data.frame(
  sample = bed_names,
  group  = assign_group(bed_names),
  stringsAsFactors = FALSE
)

comparison_list <- list(
  "FN_vs_FP" = list(g0 = "FN", g1 = "FP", c0 = "#2166AC", c1 = "#B2182B"),
  "P_N_vs_P_P" = list(g0 = "P_N", g1 = "P_P", c0 = "#4393C3", c1 = "#D6604D"),
  "N_N_vs_N_P" = list(g0 = "N_N", g1 = "N_P", c0 = "#92C5DE", c1 = "#F4A582")
)

group_colors <- c(
  "FN" = "#2166AC", "FP" = "#B2182B",
  "P_N" = "#4393C3", "P_P" = "#D6604D",
  "N_N" = "#92C5DE", "N_P" = "#F4A582"
)

# Part A:  Delta Entropy

calc_cpg_entropy <- function(beta) {
  beta <- pmax(pmin(beta, 1 - 1e-10), 1e-10)
  -(beta * log2(beta) + (1 - beta) * log2(1 - beta))
}

for (comp_name in names(comparison_list)) {
  comp <- comparison_list[[comp_name]]
  comp_dir <- file.path(hetero_dir, "DeltaEntropy", comp_name)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  s0 <- sample_info$sample[sample_info$group == comp$g0]
  s1 <- sample_info$sample[sample_info$group == comp$g1]
  all_s <- c(s0, s1)

  common <- Reduce(intersect, lapply(all_s, function(s) beta_list[[s]]$site_id))
  message(comp_name, ": ", length(common), " shared sites")

  bmat <- matrix(NA, length(common), length(all_s), dimnames = list(common, all_s))
  for (s in all_s) {
    idx <- match(common, beta_list[[s]]$site_id)
    bmat[, s] <- beta_list[[s]]$beta[idx]
  }

  mean_b0 <- rowMeans(bmat[, s0, drop = FALSE], na.rm = TRUE)
  mean_b1 <- rowMeans(bmat[, s1, drop = FALSE], na.rm = TRUE)
  ent_0 <- calc_cpg_entropy(mean_b0)
  ent_1 <- calc_cpg_entropy(mean_b1)
  delta_ent <- ent_1 - ent_0

  site_df <- data.frame(
    site_id       = common,
    beta_g0       = mean_b0,
    beta_g1       = mean_b1,
    entropy_g0    = ent_0,
    entropy_g1    = ent_1,
    delta_entropy = delta_ent
  )

  site_df$beta_bin <- cut(site_df$beta_g0,
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
    include.lowest = TRUE
  )

  set.seed(42)
  sub_df <- site_df[sample(nrow(site_df), min(80000, nrow(site_df))), ]

  p1 <- ggplot(sub_df, aes(x = beta_g0, y = delta_entropy)) +
    geom_point(aes(color = beta_bin), alpha = 0.08, size = 0.3) +
    geom_smooth(method = "loess", color = "black", linewidth = 1.2, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      "0-0.2" = "#313695", "0.2-0.4" = "#4575B4",
      "0.4-0.6" = "#ABD9E9", "0.6-0.8" = "#FDAE61", "0.8-1" = "#A50026"
    ), name = expression(beta~"bin")) +
    labs(
      title = paste0("Epigenetic Erosion: ", comp_name),
      subtitle = paste0("\u0394Entropy = ", comp$g1, " \u2212 ", comp$g0,
                        "   (>0 = entropy gain = disorder increase)"),
      x = bquote("Mean " ~ beta ~ " in " ~ .(comp$g0)),
      y = "\u0394 Shannon Entropy"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p1, file.path(comp_dir, "01_DeltaEntropy_colored_scatter"), 9, 6)

  # --- A2:  beta  delta entropy  ---
  p2 <- ggplot(site_df, aes(x = beta_bin, y = delta_entropy, fill = beta_bin)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c(
      "0-0.2" = "#313695", "0.2-0.4" = "#4575B4",
      "0.4-0.6" = "#ABD9E9", "0.6-0.8" = "#FDAE61", "0.8-1" = "#A50026"
    )) +
    coord_cartesian(ylim = quantile(site_df$delta_entropy, c(0.01, 0.99))) +
    labs(
      title = paste0("\u0394Entropy by Baseline Methylation: ", comp_name),
      x = bquote("Baseline " ~ beta ~ " bin (" ~ .(comp$g0) ~ ")"),
      y = paste0("\u0394 Entropy (", comp$g1, " \u2212 ", comp$g0, ")")
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p2, file.path(comp_dir, "02_DeltaEntropy_by_betaBin"), 8, 6)

  bin_stats <- site_df %>%
    group_by(beta_bin) %>%
    summarise(
      n_sites      = n(),
      mean_delta   = round(mean(delta_entropy), 6),
      median_delta = round(median(delta_entropy), 6),
      pct_gain     = round(100 * mean(delta_entropy > 0), 1),
      pct_loss     = round(100 * mean(delta_entropy < 0), 1),
      p_value      = wilcox.test(delta_entropy, mu = 0)$p.value,
      .groups = "drop"
    )
  write.csv(bin_stats, file.path(comp_dir, "DeltaEntropy_binStats.csv"), row.names = FALSE)
  cat("\n", comp_name, "delta entropy by interval:\n")
  print(as.data.frame(bin_stats))
}

# Part B: PIM (Proportion of Intermediate Methylation)

pim_dir <- file.path(hetero_dir, "PIM")
dir.create(pim_dir, showWarnings = FALSE, recursive = TRUE)

pim_results <- data.frame()

for (s in bed_names) {
  b <- beta_list[[s]]$beta
  b <- b[!is.na(b)]

  n_total <- length(b)
  n_low   <- sum(b < 0.2)
  n_mid   <- sum(b >= 0.2 & b <= 0.8)
  n_high  <- sum(b > 0.8)

  pim_results <- rbind(pim_results, data.frame(
    sample    = s,
    group     = assign_group(s),
    n_total   = n_total,
    n_low     = n_low,
    n_mid     = n_mid,
    n_high    = n_high,
    PIM       = n_mid / n_total,
    pct_low   = n_low / n_total,
    pct_high  = n_high / n_total,
    mean_beta = mean(b),
    var_beta  = var(b),
    stringsAsFactors = FALSE
  ))
}

pim_results$group <- factor(pim_results$group,
                            levels = c("FN", "FP", "P_N", "P_P", "N_N", "N_P"))
write.csv(pim_results, file.path(pim_dir, "PIM_per_sample.csv"), row.names = FALSE)

print(pim_results[, c("sample", "group", "PIM", "mean_beta", "var_beta")])

# --- B1: PIM  ---
p_pim <- ggplot(pim_results, aes(x = group, y = PIM, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(c("FN", "FP"), c("P_N", "P_P"), c("N_N", "N_P")),
    method = "wilcox.test", label = "p.format", size = 4
  ) +
  labs(
    title = "Proportion of Intermediate Methylation (PIM)",
    subtitle = "Fraction of CpGs with 0.2 < \u03B2 < 0.8 \u2014 measures cell population heterogeneity",
    y = "PIM", x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
save_plot(p_pim, file.path(pim_dir, "PIM_Boxplot_AllGroups"), 10, 7)

# --- B2: Beta  ---
p_var <- ggplot(pim_results, aes(x = group, y = var_beta, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(c("FN", "FP"), c("P_N", "P_P"), c("N_N", "N_P")),
    method = "wilcox.test", label = "p.format", size = 4
  ) +
  labs(
    title = "Methylation Variance Across Groups",
    subtitle = "Higher variance = more diverse methylation states across CpG sites",
    y = "Variance of \u03B2", x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
save_plot(p_var, file.path(pim_dir, "Beta_Variance_Boxplot"), 10, 7)

# --- B3: PIM vs Mean Beta  ---
p_pim_vs_beta <- ggplot(pim_results, aes(x = mean_beta, y = PIM, color = group)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), vjust = -1, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PIM vs Mean Methylation",
    subtitle = "If PIM varies independently of mean \u03B2, it captures unique heterogeneity info",
    x = "Mean \u03B2", y = "PIM"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_plot(p_pim_vs_beta, file.path(pim_dir, "PIM_vs_MeanBeta_scatter"), 8, 6)

density_list <- list()
for (s in bed_names) {
  g <- assign_group(s)
  if (is.na(g)) next
  set.seed(42)
  b <- beta_list[[s]]$beta
  b <- b[!is.na(b)]
  if (length(b) > 200000) b <- sample(b, 200000)
  density_list[[s]] <- data.frame(beta = b, sample = s, group = g)
}
density_all <- bind_rows(density_list)
density_all$group <- factor(density_all$group,
                            levels = c("FN", "FP", "P_N", "P_P", "N_N", "N_P"))

p_dens <- ggplot(density_all, aes(x = beta, color = group)) +
  geom_density(linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = c(0.2, 0.8), linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 0.2, xmax = 0.8, ymin = 0, ymax = Inf,
           alpha = 0.08, fill = "red") +
  annotate("text", x = 0.5, y = Inf, vjust = 1.5, label = "Intermediate zone",
           color = "red", fontface = "italic", size = 4) +
  labs(
    title = "Beta Value Distribution by Group",
    subtitle = "Red zone = intermediate methylation (heterogeneous)",
    x = "\u03B2 value", y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_plot(p_dens, file.path(pim_dir, "Beta_Distribution_AllGroups"), 10, 6)

cat("\n============================================================\n")
cat("============================================================\n")
cat("\nPart A (Delta Entropy): ", file.path(hetero_dir, "DeltaEntropy"), "\n")
cat("Part B (PIM + Variance):", pim_dir, "\n")
message("\n>>> Done!")

# --- 6h: DMC classification + enrichment (core imprint / stable memory / no memory) ---

# 09_MyData_Targeted_Enrichment_Manual_v2.R

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(tidyverse)
# library(gground)  # NOTE: removed, not a valid R package 

# ======================== 1.  () ========================
render_thesis_plot <- function(gene_ids, title_prefix) {
  if(length(gene_ids) < 5) return(NULL)
  
  ego <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "ALL", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  ekegg <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.05)
  
  GO_df <- as.data.frame(ego)
  KEGG_df <- as.data.frame(ekegg)
  
  if (nrow(GO_df) > 0) {
    if (!"ONTOLOGY" %in% colnames(GO_df)) { GO_df$ONTOLOGY <- ego@ontology }
    go_subset <- GO_df %>% group_by(ONTOLOGY) %>% arrange(p.adjust) %>% slice_head(n = 5) %>% ungroup()
  } else { go_subset <- data.frame() }
  
  if (nrow(KEGG_df) > 0) {
    KEGG_df$ONTOLOGY <- "KEGG"
    kegg_subset <- KEGG_df %>% arrange(p.adjust) %>% slice_head(n = 5)
  } else { kegg_subset <- data.frame() }
  
  plot_data <- bind_rows(go_subset, kegg_subset)
  if(nrow(plot_data) == 0) return(NULL)
  
  plot_data <- plot_data %>%
    separate(GeneRatio, into = c("num", "den"), sep = "/") %>%
    mutate(RichFactor = as.numeric(num) / as.numeric(den)) %>%
    mutate(Description = str_wrap(Description, width = 45)) %>% 
    mutate(Description = fct_reorder(Description, RichFactor)) 
  
  plot_data$ONTOLOGY <- factor(plot_data$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))
  
  p <- ggplot(plot_data, aes(x = RichFactor, y = Description)) +
    
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = "grey30", stroke = 0.6, alpha = 0.85) +
    
    scale_fill_gradientn(colours = c("#b2182b", "#d6604d", "#f4a582", "#92c5de", "#2166ac"),
                         name = "p.adjust", guide = guide_colorbar(reverse = TRUE)) +
    
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    labs(title = title_prefix, x = "Rich Factor", y = NULL) +
    
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey75"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 10, lineheight = 0.8),
      axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 15)),
      strip.background = element_rect(fill = "#f4f4f4", color = "grey40", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      panel.border = element_rect(color = "grey40", linewidth = 0.8),
      legend.position = "right",
      legend.key = element_blank()
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  
  return(p)
}
# ======================== 2.  () ========================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
promoters_gr <- promoters(genes_gr, upstream = 3000, downstream = 3000)

get_native_annotation <- function(query_gr, target_genes, target_promoters) {
  if(length(query_gr) == 0) return(list(all_entrez=c(), prom_entrez=c()))
  hits_all <- distanceToNearest(query_gr, target_genes)
  all_entrez <- names(target_genes)[subjectHits(hits_all)]
  hits_promoter <- findOverlaps(query_gr, target_promoters)
  prom_entrez <- names(target_promoters)[subjectHits(hits_promoter)]
  return(list(all_entrez = unique(na.omit(all_entrez)), prom_entrez = unique(na.omit(prom_entrez))))
}

# ======================== 3.  () ========================
my_dmc_dir <- "/home/root2/methylation/output/dma"
my_dmr_dir <- "/home/root2/methylation/output/Batch_DMR_Final_v5"
all_plots <- list()

process_my_data <- function(file_path, type_name) {
  comp <- basename(dirname(file_path)) 
  message(paste0("\n=== Analyzing comparison: ", comp, " [", type_name, "] ==="))
  df <- read.csv(file_path)
  
  hyper_df <- df %>% filter(meth.diff > 0)
  hypo_df  <- df %>% filter(meth.diff < 0)
  
  run_enrich <- function(sub_df, direction) {
    if(nrow(sub_df) == 0) return()
    gr <- makeGRangesFromDataFrame(sub_df, keep.extra.columns = TRUE)
    res <- get_native_annotation(gr, genes_gr, promoters_gr)
    
    # Global
    title_g <- paste0(comp, " | ", type_name, " | ", direction, " (Global)")
    p_g <- render_thesis_plot(res$all_entrez, title_g)
    if(!is.null(p_g)) all_plots[[title_g]] <<- p_g
    
    # Promoter
    title_p <- paste0(comp, " | ", type_name, " | ", direction, " (Promoter)")
    p_p <- render_thesis_plot(res$prom_entrez, title_p)
    if(!is.null(p_p)) all_plots[[title_p]] <<- p_p
  }
  run_enrich(hyper_df, "Hyper"); run_enrich(hypo_df, "Hypo")
}

# ======================== 4.  ========================
dmc_files <- list.files(my_dmc_dir, pattern = "Significant_DMCs_all.csv", recursive = TRUE, full.names = TRUE)
for(f in dmc_files) { process_my_data(f, "DMC") }

dmr_files <- list.files(my_dmr_dir, pattern = "DMR_Significant.csv", recursive = TRUE, full.names = TRUE)
for(f in dmr_files) {
  if(!grepl("_vs_", basename(dirname(f)))) next
  process_my_data(f, "DMR")
}

# 09b_Export_Enrichment_Tables.R

out_table_dir <- "/home/root2/methylation/output/Enrichment_Tables"
if(!dir.exists(out_table_dir)) dir.create(out_table_dir, recursive = TRUE)

export_enrichment_tables <- function(file_path, type_name) {
  comp <- basename(dirname(file_path)) 
  df <- read.csv(file_path)
  
  hyper_df <- df %>% filter(meth.diff > 0)
  hypo_df  <- df %>% filter(meth.diff < 0)
  
  run_export <- function(sub_df, direction) {
    if(nrow(sub_df) == 0) return()
    
    gr <- makeGRangesFromDataFrame(sub_df, keep.extra.columns = TRUE)
    res <- get_native_annotation(gr, genes_gr, promoters_gr)
    
    save_tables <- function(gene_ids, region_level) {
      if(length(gene_ids) < 5) return()
      
      ego <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "ALL", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
      ekegg <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.05)
      if (!is.null(ekegg) && nrow(ekegg) > 0) ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      GO_df <- as.data.frame(ego)
      KEGG_df <- as.data.frame(ekegg)
      
      prefix <- paste0(comp, "_", type_name, "_", direction, "_", region_level)
      
      if (nrow(GO_df) > 0) {
        write.csv(GO_df, file.path(out_table_dir, paste0(prefix, "_GO.csv")), row.names = FALSE)
      }
      
      if (nrow(KEGG_df) > 0) {
        write.csv(KEGG_df, file.path(out_table_dir, paste0(prefix, "_KEGG.csv")), row.names = FALSE)
      }
    }
    
    save_tables(res$all_entrez, "Global")
    save_tables(res$prom_entrez, "Promoter")
  }
  
  run_export(hyper_df, "Hyper")
  run_export(hypo_df, "Hypo")
}

for(f in list.files(my_dmc_dir, pattern = "Significant_DMCs_all.csv", recursive = TRUE, full.names = TRUE)) {
  export_enrichment_tables(f, "DMC")
}

for(f in list.files(my_dmr_dir, pattern = "DMR_Significant.csv", recursive = TRUE, full.names = TRUE)) {
  if(!grepl("_vs_", basename(dirname(f)))) next
  export_enrichment_tables(f, "DMR")
}

message(out_table_dir)

# --- 6i,j: KM survival curves for core imprint score ---
# --- 6k: Time-dependent AUC ---

#!/usr/bin/env Rscript
# CSC Core Imprint Complete Pipeline (Final Version)
# mean_diff scoring, AUC-focused, 5 cohorts, early/late recurrence

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(timeROC)
  library(pROC)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# KM Plot Style (matching reference figure)
# Colors: steel blue for High-risk, golden yellow for Low-risk
KM_COL <- c("group=High" = "#4472C4", "group=Low" = "#DAA520")
KM_PAL <- c("#4472C4", "#DAA520")

plot_km <- function(fit, df, title_str, time_var, event_var, hr_val, p_val, 
                    fname, xlim_val = NULL, break_by = NULL, time_unit = "months") {
  
  # palette order must match factor levels: Low first, High second
  pal <- c("#DAA520", "#4472C4")  # Low=gold, High=blue
  
  p <- ggsurvplot(fit, data = df,
    title = title_str,
    palette = pal,
    conf.int = TRUE,
    conf.int.alpha = 0.2,
    surv.median.line = "hv",
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.fontsize = 3.5,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    legend.labs = c("Low", "High"),
    legend.title = "Strata",
    xlab = "Time (months)",
    ylab = "Survival Probability",
    ggtheme = theme_classic(base_size = 12),
    font.title = c(13, "bold"),
    font.x = c(11), font.y = c(11),
    xlim = xlim_val,
    break.time.by = break_by,
    xscale = if (time_unit == "days") 30.44 else 1
  )
  
  # Add HR and p-value annotation inside the plot area
  hr_text <- sprintf("HR=%.2f  p=%s", hr_val, signif(p_val, 3))
  # Use relative coordinates via annotation_custom to stay inside plot
  p$plot <- p$plot + 
    annotate("text", x = -Inf, y = 0.1, label = hr_text, 
             size = 4, fontface = "italic", hjust = -0.1, vjust = 0)
  
  pdf(fname, width = 5, height = 5)
  print(p)
  dev.off()
}

OUT <- "/home/download/csc_article/fig5/Recurrence/Final/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
HCCDB <- "/data1/hccbulk/HCCDB_all"

# PART 1: Signature Definition
cat("========== PART 1: Core Imprint Signature ==========\n")
memory <- read.csv("/home/download/csc_article/fig5/Epigenetic_Integration/Table_Gene_Memory.csv",
                   stringsAsFactors = FALSE)
erosion <- read.csv("/home/download/csc_article/fig5/Drug_Screening/Table_Erosion_Guided_Targets.csv",
                    stringsAsFactors = FALSE)

ci <- memory[memory$pct_core == 100 & memory$is_deg == TRUE, ]
ci2 <- merge(ci, erosion[, c("gene", "rna_log2FC", "rna_padj")], by = "gene", all.x = TRUE)
ci2$meth_dir <- ifelse(ci2$mean_baseline > 0, "Hyper", "Hypo")
ci2$rna_dir  <- ifelse(is.na(ci2$rna_log2FC) | ci2$rna_log2FC <= 0, "DOWN", "UP")
has_padj <- !is.na(ci2$rna_padj) & ci2$rna_padj < 0.05

sil <- ci2$gene[ci2$meth_dir == "Hyper" & ci2$rna_dir == "DOWN" & has_padj]
act <- ci2$gene[ci2$meth_dir == "Hypo"  & ci2$rna_dir == "UP"   & has_padj]
cat("  CSC Silenced (Hyper+DOWN):", length(sil), "->", paste(sil, collapse = ", "), "\n")
cat("  CSC Activated (Hypo+UP):", length(act), "->", paste(act, collapse = ", "), "\n")
cat("  Total:", length(sil) + length(act), "genes\n\n")

# Save gene table
sig_table <- ci2[has_padj & ((ci2$meth_dir == "Hyper" & ci2$rna_dir == "DOWN") |
                              (ci2$meth_dir == "Hypo"  & ci2$rna_dir == "UP")), ]
sig_table$Group <- ifelse(sig_table$meth_dir == "Hyper", "CSC_Silenced", "CSC_Activated")
write.csv(sig_table, file.path(OUT, "Table_CoreImprint_16genes.csv"), row.names = FALSE)

compute_risk <- function(expr) {
  u <- intersect(act, rownames(expr))
  d <- intersect(sil, rownames(expr))
  if (length(u) < 1 || length(d) < 1) return(list(risk = rep(NA, ncol(expr)), nu = length(u), nd = length(d)))
  ez <- t(scale(t(expr)))
  list(risk = colMeans(ez[u, , drop = FALSE], na.rm = TRUE) -
              colMeans(ez[d, , drop = FALSE], na.rm = TRUE),
       nu = length(u), nd = length(d))
}

# PART 2: Load Cohorts
cat("========== PART 2: Load Cohorts ==========\n")

load_hccdb <- function(expr_file, sample_file, patient_file) {
  expr <- read.delim(expr_file, stringsAsFactors = FALSE, check.names = FALSE)
  mat <- as.matrix(expr[, -(1:2)])
  rownames(mat) <- expr$Symbol
  bad <- duplicated(expr$Symbol) | is.na(expr$Symbol) | expr$Symbol == ""
  mat <- mat[!bad, ]
  sam <- read.delim(sample_file, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
  sam_t <- as.data.frame(t(sam[, -1]), stringsAsFactors = FALSE)
  colnames(sam_t) <- sam[, 1]
  pat <- read.delim(patient_file, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
  pat_t <- as.data.frame(t(pat[, -1]), stringsAsFactors = FALSE)
  colnames(pat_t) <- pat[, 1]
  list(mat = mat, sam = sam_t, pat = pat_t)
}

# --- TCGA expression matrix (shared) ---
cat("  TCGA expression...\n")
tcga_raw <- read.csv("/data1/WGCNA/CellType_Specific_Results/TCGA_Survival/LIHC_fpkm.csv",
                     stringsAsFactors = FALSE, check.names = FALSE)
gc <- tcga_raw[, 1]; tm <- as.matrix(tcga_raw[, -1]); rownames(tm) <- gc
dup <- duplicated(gc)
if (any(dup)) {
  mn <- rowMeans(tm, na.rm = TRUE); kp <- rep(TRUE, nrow(tm))
  for (g in unique(gc[dup])) { i <- which(gc == g); b <- i[which.max(mn[i])]; kp[i] <- FALSE; kp[b] <- TRUE }
  tm <- tm[kp, ]; rownames(tm) <- gc[kp]
}
tc <- grep("-01[A-Z]", colnames(tm), value = TRUE)
tm <- tm[, tc]; colnames(tm) <- substr(colnames(tm), 1, 12); tm <- tm[, !duplicated(colnames(tm))]

# --- TCGA OS ---
cat("  TCGA OS...\n")
cl <- read.csv("/data1/WGCNA/CellType_Specific_Results/TCGA_Survival/LIHC_clinical.csv",
               stringsAsFactors = FALSE)
cl$pid <- substr(cl$barcode, 1, 12)
cl$time <- as.numeric(ifelse(cl$vital_status == "Dead", cl$days_to_death, cl$days_to_last_follow_up))
cl$event <- ifelse(cl$vital_status == "Dead", 1, 0)
cm <- intersect(colnames(tm), cl$pid); cl_m <- cl[match(cm, cl$pid), ]
expr_tcga_os <- tm[, cm]
ok <- !is.na(cl_m$time) & cl_m$time > 0 & !is.na(cl_m$event)
expr_tcga_os <- expr_tcga_os[, ok]; time_tcga_os <- cl_m$time[ok]; ev_tcga_os <- cl_m$event[ok]
risk_tcga_os <- compute_risk(expr_tcga_os)
cat("    n=", ncol(expr_tcga_os), " events=", sum(ev_tcga_os), " genes:", risk_tcga_os$nu, "+", risk_tcga_os$nd, "\n")

# --- TCGA DFS ---
cat("  TCGA DFS...\n")
dfs <- read.delim("~/JWD_LIHC/lihc_tcga_pan_can_atlas_2018_clinical_data.tsv",
                  stringsAsFactors = FALSE, check.names = FALSE)
dfs$dm <- as.numeric(dfs[["Disease Free (Months)"]])
dfs$de <- ifelse(grepl("Recurred", dfs[["Disease Free Status"]]), 1, 0)
if ("Patient ID" %in% colnames(dfs)) dfs$pid <- dfs[["Patient ID"]] else {
  for (cc in colnames(dfs)) if (any(grepl("^TCGA-", dfs[[cc]]))) { dfs$pid <- substr(dfs[[cc]], 1, 12); break }
}
dfs <- dfs[!duplicated(dfs$pid), ]
cm2 <- intersect(colnames(tm), dfs$pid); dfs_m <- dfs[match(cm2, dfs$pid), ]
ok2 <- !is.na(dfs_m$dm) & dfs_m$dm > 0 & !is.na(dfs_m$de)
expr_tcga_dfs <- tm[, cm2[ok2]]; rfs_m_tcga <- dfs_m$dm[ok2]; rfs_e_tcga <- dfs_m$de[ok2]
risk_tcga_dfs <- compute_risk(expr_tcga_dfs)
cat("    n=", ncol(expr_tcga_dfs), " events=", sum(rfs_e_tcga), "\n")

# --- GSE14520 RFS ---
cat("  GSE14520...\n")
expr_gse <- readRDS("/home/download/csc_article/fig5/Recurrence/GSE14520/GSE14520_expr.rds")
supp <- read.delim("/home/download/csc_article/fig6/GSE14520_Extra_Supplement.txt", stringsAsFactors = FALSE)
cg <- intersect(colnames(expr_gse), supp$Affy_GSM)
supp <- supp[match(cg, supp$Affy_GSM), ]; expr_gse <- expr_gse[, cg]
rfs_m_gse <- as.numeric(supp$Recurr.months); rfs_e_gse <- as.numeric(supp$Recurr.status)
ok3 <- !is.na(rfs_m_gse) & rfs_m_gse >= 0 & !is.na(rfs_e_gse)
expr_gse <- expr_gse[, ok3]; rfs_m_gse <- rfs_m_gse[ok3]; rfs_e_gse <- rfs_e_gse[ok3]
risk_gse <- compute_risk(expr_gse)
cat("    n=", ncol(expr_gse), " events=", sum(rfs_e_gse), " genes:", risk_gse$nu, "+", risk_gse$nd, "\n")

# --- HCCDB18 / ICGC LIRI-JP OS ---
cat("  HCCDB18/ICGC...\n")
d18 <- load_hccdb(file.path(HCCDB, "HCCDB18_mRNA_level3.txt"),
                  file.path(HCCDB, "sample/HCCDB18.sample.txt"),
                  file.path(HCCDB, "patient/HCCDB18.patient.txt"))
hcc18 <- d18$sam$TYPE == "HCC"; sid18 <- d18$sam$SAMPLE_ID[hcc18]; pid18 <- d18$sam$PATIENT_ID[hcc18]
expr_18 <- d18$mat[, intersect(colnames(d18$mat), sid18)]
colnames(expr_18) <- setNames(pid18, sid18)[colnames(expr_18)]
expr_18 <- expr_18[, !duplicated(colnames(expr_18))]
ri18 <- match(colnames(expr_18), d18$pat$PATIENT_ID)
time_18 <- as.numeric(d18$pat$SUR[ri18]) * 30.44
ev_18 <- ifelse(d18$pat$STATUS[ri18] == "Dead" | d18$pat$STATUS[ri18] == "Deceased", 1, 0)
ok4 <- !is.na(time_18) & time_18 > 0 & !is.na(ev_18)
expr_18 <- expr_18[, ok4]; time_18 <- time_18[ok4]; ev_18 <- ev_18[ok4]
risk_18 <- compute_risk(expr_18)
cat("    n=", ncol(expr_18), " events=", sum(ev_18), " genes:", risk_18$nu, "+", risk_18$nd, "\n")

# --- HCCDB17 / GSE76427 RFS ---
cat("  HCCDB17...\n")
d17 <- load_hccdb(file.path(HCCDB, "HCCDB17_mRNA_level3.txt"),
                  file.path(HCCDB, "sample/HCCDB17.sample.txt"),
                  file.path(HCCDB, "patient/HCCDB17.patient.txt"))
hcc17 <- d17$sam$TYPE == "HCC"; sid17 <- d17$sam$SAMPLE_ID[hcc17]; pid17 <- d17$sam$PATIENT_ID[hcc17]
expr_17 <- d17$mat[, intersect(colnames(d17$mat), sid17)]
colnames(expr_17) <- setNames(pid17, sid17)[colnames(expr_17)]
ri17 <- match(colnames(expr_17), d17$pat$PATIENT_ID)
rfs_m_17 <- as.numeric(d17$pat$RFS_YEAR[ri17]) * 12
rfs_e_17 <- ifelse(d17$pat$RFS_STATUS[ri17] == "Yes", 1, 0)
ok5 <- !is.na(rfs_m_17) & rfs_m_17 > 0 & !is.na(rfs_e_17)
expr_17 <- expr_17[, ok5]; rfs_m_17 <- rfs_m_17[ok5]; rfs_e_17 <- rfs_e_17[ok5]
risk_17 <- compute_risk(expr_17)
cat("    n=", ncol(expr_17), " events=", sum(rfs_e_17), " genes:", risk_17$nu, "+", risk_17$nd, "\n")

# PART 3: OS Validation (TCGA + HCCDB18)
cat("\n========== PART 3: Overall Survival ==========\n")

validate_os <- function(risk, time, event, name) {
  df <- data.frame(time = time, event = event, risk = risk)
  df <- df[!is.na(df$risk), ]
  cx <- coxph(Surv(time, event) ~ risk, data = df)
  hr <- exp(coef(cx)); ci_v <- exp(confint(cx)); pv <- summary(cx)$coefficients[, "Pr(>|z|)"]
  # KM
  df$grp <- factor(ifelse(df$risk > median(df$risk), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ grp, data = df)
  lr <- survdiff(Surv(time, event) ~ grp, data = df); p_lr <- 1 - pchisq(lr$chisq, 1)
  hr_g <- exp(coef(coxph(Surv(time, event) ~ grp, data = df)))  # High vs Low, HR > 1
  # AUC
  tp <- c(365, 730, 1095, 1825)
  troc <- tryCatch(timeROC(T = df$time, delta = df$event, marker = df$risk,
                            cause = 1, times = tp, iid = FALSE), error = function(e) NULL)
  auc <- if (!is.null(troc)) round(troc$AUC, 3) else rep(NA, 4)
  
  cat(sprintf("  %s: n=%d ev=%d HR=%.3f(%.3f-%.3f) p=%.1e p_lr=%.1e AUC1y=%.3f 3y=%.3f\n",
              name, nrow(df), sum(df$event), hr, ci_v[1], ci_v[2], pv, p_lr, auc[1], auc[3]))
  
  # KM plot
  plot_km(fit, df, paste0(name, ": Overall Survival"), 
          "time", "event", hr_g, p_lr,
          file.path(OUT, paste0("KM_OS_", gsub("[^A-Za-z0-9]", "_", name), ".pdf")),
          time_unit = "days")
  
  data.frame(Cohort = name, Endpoint = "OS", N = nrow(df), Events = sum(df$event),
             HR = round(hr, 3), CI_lo = round(ci_v[1], 3), CI_hi = round(ci_v[2], 3),
             P = signif(pv, 3), HR_grp = round(hr_g, 3), P_lr = signif(p_lr, 3),
             AUC_1y = auc[1], AUC_2y = auc[2], AUC_3y = auc[3], AUC_5y = auc[4])
}

os_res <- list()
os_res[[1]] <- validate_os(risk_tcga_os$risk, time_tcga_os, ev_tcga_os, "TCGA-LIHC")
os_res[[2]] <- validate_os(risk_18$risk, time_18, ev_18, "ICGC LIRI-JP")

# PART 4: RFS Validation (GSE14520 + TCGA_DFS + HCCDB17)
cat("\n========== PART 4: Recurrence-Free Survival ==========\n")

validate_rfs <- function(risk, rfs_months, rfs_event, name) {
  rfs_days <- rfs_months * 30.44
  df <- data.frame(m = rfs_months, d = rfs_days, ev = rfs_event, risk = risk)
  df <- df[!is.na(df$risk), ]
  
  # Overall RFS
  cx <- coxph(Surv(d, ev) ~ risk, data = df)
  hr <- exp(coef(cx)); ci_v <- exp(confint(cx)); pv <- summary(cx)$coefficients[, "Pr(>|z|)"]
  df$grp <- factor(ifelse(df$risk > median(df$risk), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(d, ev) ~ grp, data = df)
  lr <- survdiff(Surv(d, ev) ~ grp, data = df); p_lr <- 1 - pchisq(lr$chisq, 1)
  hr_g <- exp(coef(coxph(Surv(d, ev) ~ grp, data = df)))  # High vs Low
  
  tp <- c(365, 730, 1095, 1825)
  troc <- tryCatch(timeROC(T = df$d, delta = df$ev, marker = df$risk,
                            cause = 1, times = tp, iid = FALSE), error = function(e) NULL)
  auc <- if (!is.null(troc)) round(troc$AUC, 3) else rep(NA, 4)
  
  cat(sprintf("  %s: n=%d ev=%d HR=%.3f(%.3f-%.3f) p=%.1e p_lr=%.1e AUC1y=%.3f 3y=%.3f\n",
              name, nrow(df), sum(df$ev), hr, ci_v[1], ci_v[2], pv, p_lr, auc[1], auc[3]))
  
  # KM
  plot_km(fit, df, paste0(name, ": Recurrence-Free Survival"),
          "d", "ev", hr_g, p_lr,
          file.path(OUT, paste0("KM_RFS_", gsub("[^A-Za-z0-9]", "_", name), ".pdf")),
          time_unit = "days")
  
  list(overall = data.frame(Cohort = name, Endpoint = "RFS", N = nrow(df), Events = sum(df$ev),
                             HR = round(hr, 3), CI_lo = round(ci_v[1], 3), CI_hi = round(ci_v[2], 3),
                             P = signif(pv, 3), HR_grp = round(hr_g, 3), P_lr = signif(p_lr, 3),
                             AUC_1y = auc[1], AUC_2y = auc[2], AUC_3y = auc[3], AUC_5y = auc[4]),
       df = df)
}

rfs_res <- list()
rfs_res[["GSE14520"]]  <- validate_rfs(risk_gse$risk, rfs_m_gse, rfs_e_gse, "GSE14520")
rfs_res[["TCGA_DFS"]]  <- validate_rfs(risk_tcga_dfs$risk, rfs_m_tcga, rfs_e_tcga, "TCGA-DFS")
rfs_res[["GSE76427"]]   <- validate_rfs(risk_17$risk, rfs_m_17, rfs_e_17, "GSE76427")

# PART 5: Early vs Late Recurrence
cat("\n========== PART 5: Early vs Late Recurrence ==========\n")

early_late <- function(df, name) {
  cat(sprintf("\n  --- %s (n=%d, events=%d) ---\n", name, nrow(df), sum(df$ev)))
  
  # Event counts
  cat("    Events by time:\n")
  for (m in c(3, 6, 9, 12, 18, 24, 36, 48, 60)) {
    n_r <- sum(df$ev == 1 & df$m <= m)
    if (n_r > 0) cat(sprintf("      %2dm: %3d (%.1f%%)\n", m, n_r, n_r / nrow(df) * 100))
  }
  
  # Early (<=24m)
  df_e <- df
  df_e$t <- pmin(df_e$d, 24 * 30.44)
  df_e$e <- ifelse(df_e$ev == 1 & df_e$m <= 24, 1, 0)
  ne <- sum(df_e$e)
  hr_e <- NA; p_e <- NA; c_e <- NA
  if (ne >= 5) {
    cx_e <- coxph(Surv(t, e) ~ risk, data = df_e)
    hr_e <- exp(coef(cx_e)); ci_e <- exp(confint(cx_e)); p_e <- summary(cx_e)$coefficients[, "Pr(>|z|)"]
    c_e <- concordance(cx_e)$concordance
    cat(sprintf("    Early(<=2y): ev=%d HR=%.3f(%.3f-%.3f) p=%.1e C=%.3f\n", ne, hr_e, ci_e[1], ci_e[2], p_e, c_e))
    
    # KM early
    df_e$grp <- factor(ifelse(df_e$risk > median(df_e$risk), "High", "Low"), levels = c("Low", "High"))
    fit <- survfit(Surv(t, e) ~ grp, data = df_e)
    lr <- survdiff(Surv(t, e) ~ grp, data = df_e); p_lr <- 1 - pchisq(lr$chisq, 1)
    hr_g <- exp(coef(coxph(Surv(t, e) ~ grp, data = df_e)))  # High vs Low
    plot_km(fit, df_e, paste0(name, ": Early Recurrence (<=2y)"),
            "t", "e", hr_g, p_lr,
            file.path(OUT, paste0("KM_Early_", gsub("[^A-Za-z0-9]", "_", name), ".pdf")),
            xlim_val = c(0, 24 * 30.44), break_by = 180, time_unit = "days")
  }
  
  # Late (>24m)
  df_l <- df[df$m > 24, ]; df_l$t <- df_l$d - 24 * 30.44; nl <- sum(df_l$ev)
  hr_l <- NA; p_l <- NA; c_l <- NA
  if (nl >= 5) {
    cx_l <- coxph(Surv(t, ev) ~ risk, data = df_l)
    hr_l <- exp(coef(cx_l)); ci_l <- exp(confint(cx_l)); p_l <- summary(cx_l)$coefficients[, "Pr(>|z|)"]
    c_l <- concordance(cx_l)$concordance
    cat(sprintf("    Late(>2y):  ev=%d HR=%.3f(%.3f-%.3f) p=%.1e C=%.3f\n", nl, hr_l, ci_l[1], ci_l[2], p_l, c_l))
    
    df_l$grp <- factor(ifelse(df_l$risk > median(df$risk), "High", "Low"), levels = c("Low", "High"))
    fit_l <- survfit(Surv(t, ev) ~ grp, data = df_l)
    lr_l <- survdiff(Surv(t, ev) ~ grp, data = df_l); p_lr_l <- 1 - pchisq(lr_l$chisq, 1)
    hr_g_l <- exp(coef(coxph(Surv(t, ev) ~ grp, data = df_l)))  # High vs Low
    plot_km(fit_l, df_l, paste0(name, ": Late Recurrence (>2y)"),
            "t", "ev", hr_g_l, p_lr_l,
            file.path(OUT, paste0("KM_Late_", gsub("[^A-Za-z0-9]", "_", name), ".pdf")),
            time_unit = "days")
  }
  
  # Time-dependent AUC (3-month intervals)
  months_grid <- seq(3, min(60, floor(max(df$m[df$ev == 1]))), by = 3)
  tp <- months_grid * 30.44
  troc <- tryCatch(timeROC(T = df$d, delta = df$ev, marker = df$risk,
                            cause = 1, times = tp, iid = TRUE), error = function(e) NULL)
  auc_rows <- list()
  if (!is.null(troc)) {
    cat("    Time-dependent AUC:\n")
    cat("      Month  Events  AUC    95%CI\n")
    for (i in seq_along(tp)) {
      m <- months_grid[i]
      n_ev <- sum(df$ev == 1 & df$d <= tp[i])
      se <- sqrt(troc$inference$vect_sd_1[i])
      cat(sprintf("      %4dm   %3d     %.3f  (%.3f-%.3f)\n",
                  m, n_ev, troc$AUC[i], troc$AUC[i] - 1.96 * se, troc$AUC[i] + 1.96 * se))
      auc_rows[[i]] <- data.frame(Cohort = name, Month = m, Events = n_ev,
                                   AUC = round(troc$AUC[i], 3), SE = round(se, 3))
    }
  }
  
  list(
    early = data.frame(Cohort = name, N = nrow(df), Events = ne,
                        HR = round(hr_e, 3), CI_lo = round(ci_e[1], 3), CI_hi = round(ci_e[2], 3),
                        P = p_e, C = round(c_e, 3)),
    late  = data.frame(Cohort = name, N = ifelse(nl >= 5, nrow(df_l), NA), Events = nl,
                        HR = round(hr_l, 3), CI_lo = round(ci_l[1], 3), CI_hi = round(ci_l[2], 3),
                        P = p_l, C = round(c_l, 3)),
    auc   = if (length(auc_rows) > 0) do.call(rbind, auc_rows) else NULL
  )
}

el_res <- list()
for (nm in names(rfs_res)) {
  cohort_name <- rfs_res[[nm]]$overall$Cohort
  el_res[[nm]] <- early_late(rfs_res[[nm]]$df, cohort_name)
}

# PART 6: Summary Tables
cat("\n\n========== PART 6: Summary Tables ==========\n")

# OS + RFS overview
cat("\n--- Multi-cohort Overview ---\n")
overview <- rbind(
  do.call(rbind, os_res),
  do.call(rbind, lapply(rfs_res, function(x) x$overall))
)
print(overview, row.names = FALSE)
write.csv(overview, file.path(OUT, "Table_MultiCohort_Overview.csv"), row.names = FALSE)

# Early vs Late
cat("\n--- Early vs Late Recurrence ---\n")
early_df <- do.call(rbind, lapply(el_res, function(x) x$early))
late_df  <- do.call(rbind, lapply(el_res, function(x) x$late))
cat("Early:\n"); print(early_df, row.names = FALSE)
cat("Late:\n");  print(late_df, row.names = FALSE)
write.csv(rbind(cbind(Phase = "Early", early_df), cbind(Phase = "Late", late_df)),
          file.path(OUT, "Table_EarlyLate.csv"), row.names = FALSE)

# AUC over time
auc_all <- do.call(rbind, lapply(el_res, function(x) x$auc))
write.csv(auc_all, file.path(OUT, "Table_TimeAUC.csv"), row.names = FALSE)

# PART 7: Meta-analyses
cat("\n========== PART 7: Meta-analyses ==========\n")

do_meta <- function(df_sub, label) {
  df_sub <- df_sub[!is.na(df_sub$HR) & !is.na(df_sub$P) & df_sub$HR > 0 & df_sub$P < 1, ]
  if (nrow(df_sub) < 2) { cat("  ", label, ": not enough cohorts\n"); return(NULL) }
  lhr <- log(df_sub$HR); z <- qnorm(1 - df_sub$P / 2); se <- abs(lhr / z)
  w <- 1 / se^2; pl <- sum(w * lhr) / sum(w); ps <- sqrt(1 / sum(w))
  Q <- sum(w * (lhr - pl)^2); dff <- nrow(df_sub) - 1
  I2 <- max(0, (Q - dff) / Q * 100)
  pp <- 2 * pnorm(-abs(pl / ps))
  cat(sprintf("  %s: pooled HR=%.3f(%.3f-%.3f) p=%.2e I2=%.0f%% N=%d\n",
              label, exp(pl), exp(pl - 1.96 * ps), exp(pl + 1.96 * ps), pp, I2, sum(df_sub$N)))
  list(HR = exp(pl), CI_lo = exp(pl - 1.96 * ps), CI_hi = exp(pl + 1.96 * ps), P = pp, I2 = I2)
}

# OS meta (TCGA + HCCDB18)
os_overview <- do.call(rbind, os_res)
meta_os <- do_meta(os_overview, "Overall Survival")

# RFS meta (all 3)
rfs_overview <- do.call(rbind, lapply(rfs_res, function(x) x$overall))
meta_rfs <- do_meta(rfs_overview, "Recurrence-Free Survival")

# All 5 cohorts
meta_all <- do_meta(overview, "All cohorts combined")

# Early recurrence meta
meta_early <- do_meta(early_df, "Early Recurrence")

# Late recurrence meta
meta_late <- do_meta(late_df, "Late Recurrence")

# PART 8: Figures
cat("\n========== PART 8: Publication Figures ==========\n")

# --- Figure: AUC Decay Curve ---
auc_all$CI_lo <- pmax(0.5, auc_all$AUC - 1.96 * auc_all$SE)
auc_all$CI_hi <- pmin(1.0, auc_all$AUC + 1.96 * auc_all$SE)

cols <- c("GSE14520" = "#D62728", "TCGA-DFS" = "#2166AC", "GSE76427" = "#4DAF4A")

p_auc <- ggplot(auc_all, aes(x = Month, y = AUC, color = Cohort)) +
  annotate("rect", xmin = 0, xmax = 24, ymin = 0.4, ymax = 0.85, fill = "#FFE5E5", alpha = 0.3) +
  annotate("rect", xmin = 24, xmax = 60, ymin = 0.4, ymax = 0.85, fill = "#E5E5FF", alpha = 0.3) +
  annotate("text", x = 12, y = 0.83, label = "Early Recurrence", color = "#B22222", size = 3.5, fontface = "bold") +
  annotate("text", x = 42, y = 0.83, label = "Late Recurrence", color = "#00008B", size = 3.5, fontface = "bold") +
  geom_vline(xintercept = 24, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_lo, ymax = CI_hi, fill = Cohort), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, 60, 6)) +
  scale_y_continuous(limits = c(0.4, 0.85), breaks = seq(0.4, 0.85, 0.05)) +
  labs(x = "Time (months)", y = "Time-dependent AUC",
       title = "CSC Core Imprint: Predictive Accuracy Over Time") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"),
        plot.title = element_text(face = "bold", size = 13))

pdf(file.path(OUT, "Fig_AUC_Decay.pdf"), width = 9, height = 5.5)
print(p_auc); dev.off()
cat("  -> Fig_AUC_Decay.pdf\n")

# --- Figure: Forest Plot Early vs Late ---
forest_rows <- list()
# Early
for (i in 1:nrow(early_df)) {
  forest_rows[[length(forest_rows) + 1]] <- data.frame(
    Label = paste0(early_df$Cohort[i], " (n=", early_df$N[i], ", ev=", early_df$Events[i], ")"),
    HR = early_df$HR[i], CI_lo = early_df$CI_lo[i], CI_hi = early_df$CI_hi[i],
    P = early_df$P[i], Phase = "Early (<=2y)", Type = "Study")
}
if (!is.null(meta_early)) {
  forest_rows[[length(forest_rows) + 1]] <- data.frame(
    Label = paste0("Pooled (I2=", round(meta_early$I2), "%)"),
    HR = meta_early$HR, CI_lo = meta_early$CI_lo, CI_hi = meta_early$CI_hi,
    P = meta_early$P, Phase = "Early (<=2y)", Type = "Pooled")
}
# Late
for (i in 1:nrow(late_df)) {
  if (!is.na(late_df$HR[i])) {
    forest_rows[[length(forest_rows) + 1]] <- data.frame(
      Label = paste0(late_df$Cohort[i], " (n=", late_df$N[i], ", ev=", late_df$Events[i], ")"),
      HR = late_df$HR[i], CI_lo = late_df$CI_lo[i], CI_hi = late_df$CI_hi[i],
      P = late_df$P[i], Phase = "Late (>2y)", Type = "Study")
  }
}
if (!is.null(meta_late)) {
  forest_rows[[length(forest_rows) + 1]] <- data.frame(
    Label = "Pooled",
    HR = meta_late$HR, CI_lo = meta_late$CI_lo, CI_hi = meta_late$CI_hi,
    P = meta_late$P, Phase = "Late (>2y)", Type = "Pooled")
}

forest_df <- do.call(rbind, forest_rows)
forest_df$y <- nrow(forest_df):1
forest_df$HR_text <- sprintf("%.2f (%.2f-%.2f)", forest_df$HR, forest_df$CI_lo, forest_df$CI_hi)
forest_df$CI_lo_clip <- pmax(forest_df$CI_lo, 0.1)
forest_df$CI_hi_clip <- pmin(forest_df$CI_hi, 5)

p_forest <- ggplot(forest_df, aes(y = y, x = HR)) +
  annotate("rect", xmin = 0.1, xmax = 5,
           ymin = min(forest_df$y[forest_df$Phase == "Early (<=2y)"]) - 0.4,
           ymax = max(forest_df$y[forest_df$Phase == "Early (<=2y)"]) + 0.4,
           fill = "#FFE5E5", alpha = 0.3) +
  annotate("rect", xmin = 0.1, xmax = 5,
           ymin = min(forest_df$y[forest_df$Phase == "Late (>2y)"]) - 0.4,
           ymax = max(forest_df$y[forest_df$Phase == "Late (>2y)"]) + 0.4,
           fill = "#E5E5FF", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lo_clip, xmax = CI_hi_clip), height = 0.2, linewidth = 0.7) +
  geom_point(aes(shape = Type, size = Type, color = Phase)) +
  scale_shape_manual(values = c(Study = 16, Pooled = 18)) +
  scale_size_manual(values = c(Study = 3, Pooled = 5)) +
  scale_color_manual(values = c("Early (<=2y)" = "#D62728", "Late (>2y)" = "#2166AC")) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5), limits = c(0.1, 5)) +
  scale_y_continuous(breaks = forest_df$y, labels = forest_df$Label) +
  geom_text(aes(x = 4.5, label = HR_text), hjust = 1, size = 3, color = "grey20") +
  labs(x = "Hazard Ratio (log scale)", y = "",
       title = "Early vs Late Recurrence") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

pdf(file.path(OUT, "Fig_Forest_EarlyLate.pdf"), width = 9, height = 5)
print(p_forest); dev.off()
cat("  -> Fig_Forest_EarlyLate.pdf\n")

# Combined panel
pdf(file.path(OUT, "Fig_Recurrence_Panel.pdf"), width = 14, height = 5.5)
grid.arrange(
  ggplotGrob(p_auc), ggplotGrob(p_forest),
  ncol = 2, widths = c(1.1, 0.9),
  top = textGrob("CSC Core Imprint Recurrence Prediction",
                 gp = gpar(fontface = "bold", fontsize = 14))
)
dev.off()
cat("  -> Fig_Recurrence_Panel.pdf\n")

cat("\n========== COMPLETE ==========\n")
cat("All outputs:", OUT, "\n")

# --- 6k supplement: Erosion-guided recurrence prediction ---

#!/usr/bin/env Rscript
#  Erosion-Guided Gene Signature Predicts HCC Recurrence
#
#  Hypothesis: Genes with low erosion rate (actively maintained by CSC)
#              should predict tumor recurrence better than high-erosion genes
#
#  Design:
#    Signature A: Low erosion (erosion < 0.3) = CSC actively maintained
#    Signature B: High erosion (erosion > 0.7) = passively different (control)
#    Signature C: Random genes (same size as A) = null control
#
#  Cohorts:
#    Training:    TCGA-LIHC (n=371, RFS)
#    Validation1: GSE14520  (n=242, RFS)
#    Validation2: ICGC-LIRI (n=243, OS+RFS)
#
#  Output: KM curves, Cox HR, ROC AUC comparing A vs B vs C

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(GSVA)         # ssGSEA scoring
  library(GEOquery)
  library(data.table)
  library(ggplot2)
  library(pROC)
  library(dplyr)
})

OUT_DIR <- "/home/download/csc_article/fig5/Recurrence/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# PART 0: Load erosion-guided targets
message("===== PART 0: Load erosion data =====")

erosion <- read.csv("/home/download/csc_article/fig5/Drug_Screening/Table_Erosion_Guided_Targets.csv",
                    stringsAsFactors = FALSE)

message("  Total genes: ", nrow(erosion))
message("  With erosion data: ", sum(!is.na(erosion$mean_erosion)))

# Define five directional signatures
# Need RNA direction info: Hyper+DOWN vs Hypo+UP
erosion$meth_dir <- ifelse(erosion$mean_baseline > 0, "Hyper", "Hypo")
erosion$rna_dir  <- ifelse(!is.na(erosion$rna_log2FC) & erosion$rna_log2FC > 0, "UP", "DOWN")
erosion$concordant <- (erosion$meth_dir == "Hyper" & erosion$rna_dir == "DOWN") |
                      (erosion$meth_dir == "Hypo"  & erosion$rna_dir == "UP")

message("  Direction concordant genes: ", sum(erosion$concordant, na.rm = TRUE))

# Use erosion < 0.5 and |baseline| > 5 for active signatures (more inclusive)
# Use erosion > 0.7 for passive controls
# Cap each signature at N_SIG genes by priority score
N_SIG <- 20

# A1: Low erosion + Hyper + DOWN = CSC actively silenced tumor suppressors
sig_silenced <- erosion %>%
  filter(!is.na(mean_erosion), mean_erosion < 0.5,
         mean_baseline > 5,
         meth_dir == "Hyper", rna_dir == "DOWN",
         !is.na(rna_padj), rna_padj < 0.05) %>%
  arrange(mean_erosion) %>%
  head(N_SIG) %>%
  pull(gene)

# A2: Low erosion + Hypo + UP = CSC actively activated oncogenes
sig_activated <- erosion %>%
  filter(!is.na(mean_erosion), mean_erosion < 0.5,
         mean_baseline < -5,
         meth_dir == "Hypo", rna_dir == "UP",
         !is.na(rna_padj), rna_padj < 0.05) %>%
  arrange(mean_erosion) %>%
  head(N_SIG) %>%
  pull(gene)

# B1: High erosion + Hyper + DOWN (passive control)
sig_passive_hyper <- erosion %>%
  filter(!is.na(mean_erosion), mean_erosion > 0.7,
         mean_baseline > 5,
         meth_dir == "Hyper", rna_dir == "DOWN",
         !is.na(rna_padj), rna_padj < 0.05) %>%
  arrange(desc(mean_erosion)) %>%
  head(N_SIG) %>%
  pull(gene)

# B2: High erosion + Hypo + UP (passive control)
sig_passive_hypo <- erosion %>%
  filter(!is.na(mean_erosion), mean_erosion > 0.7,
         mean_baseline < -5,
         meth_dir == "Hypo", rna_dir == "UP",
         !is.na(rna_padj), rna_padj < 0.05) %>%
  arrange(desc(mean_erosion)) %>%
  head(N_SIG) %>%
  pull(gene)

# C: Random genes (same size as N_SIG)
set.seed(42)
all_genes_with_erosion <- erosion %>%
  filter(!is.na(mean_erosion), !is.na(rna_padj), rna_padj < 0.05) %>%
  pull(gene)
used_genes <- c(sig_silenced, sig_activated, sig_passive_hyper, sig_passive_hypo)
sig_random <- sample(setdiff(all_genes_with_erosion, used_genes), size = N_SIG)

message("  CSC_Silenced  (low erosion, Hyper+DOWN): ", length(sig_silenced), " genes")
message("  CSC_Activated (low erosion, Hypo+UP):    ", length(sig_activated), " genes")
message("  Passive_Hyper (high erosion, Hyper+DOWN):", length(sig_passive_hyper), " genes")
message("  Passive_Hypo  (high erosion, Hypo+UP):   ", length(sig_passive_hypo), " genes")
message("  Random control:                          ", length(sig_random), " genes")

# Save gene lists
all_sig_info <- rbind(
  data.frame(gene = sig_silenced,     signature = "CSC_Silenced"),
  data.frame(gene = sig_activated,    signature = "CSC_Activated"),
  data.frame(gene = sig_passive_hyper, signature = "Passive_Hyper"),
  data.frame(gene = sig_passive_hypo,  signature = "Passive_Hypo"),
  data.frame(gene = sig_random,       signature = "Random")
)
write.csv(all_sig_info, file.path(OUT_DIR, "Signature_GeneLists.csv"), row.names = FALSE)

signatures <- list(
  CSC_Silenced  = sig_silenced,
  CSC_Activated = sig_activated,
  Passive_Hyper = sig_passive_hyper,
  Passive_Hypo  = sig_passive_hypo,
  Random        = sig_random
)

# PART 1: TCGA-LIHC (Training)
message("\n===== PART 1: TCGA-LIHC =====")

# ---- 1a. Download expression ----
tcga_expr_file <- "/data1/WGCNA/CellType_Specific_Results/TCGA_Survival/LIHC_fpkm.csv"
tcga_clin_file <- "/data1/WGCNA/CellType_Specific_Results/TCGA_Survival/LIHC_clinical.csv"

message("  Loading local TCGA-LIHC data...")
expr_df <- read.csv(tcga_expr_file, check.names = FALSE)
gene_names <- expr_df[[1]]
expr_df <- expr_df[, -1]

# Handle duplicate gene names: keep highest mean expression
expr_df$gene <- gene_names
expr_df <- expr_df %>%
  group_by(gene) %>%
  filter(row_number() == which.max(rowMeans(pick(where(is.numeric)), na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame()
rownames(expr_df) <- expr_df$gene
expr_df$gene <- NULL
expr_mat <- as.matrix(expr_df)

clin <- read.csv(tcga_clin_file, stringsAsFactors = FALSE)

message("  Expression: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
message("  Clinical: ", nrow(clin), " patients")

# ---- Keep tumor samples only (01A) ----
tumor_cols <- grep("-01[A-Z]", colnames(expr_mat), value = TRUE)
expr_tcga <- expr_mat[, tumor_cols]
colnames(expr_tcga) <- substr(colnames(expr_tcga), 1, 12)
expr_tcga <- expr_tcga[, !duplicated(colnames(expr_tcga))]
message("  Tumor samples: ", ncol(expr_tcga))

# ---- Prepare survival from clinical ----
clin$patient_id <- substr(clin$barcode, 1, 12)

surv_tcga <- data.frame(
  sample = clin$patient_id,
  OS_time = as.numeric(ifelse(clin$vital_status == "Dead",
                               clin$days_to_death,
                               clin$days_to_last_follow_up)),
  OS_event = ifelse(clin$vital_status == "Dead", 1, 0),
  stringsAsFactors = FALSE
)

# Remove NAs and keep tumor patients only
surv_tcga <- surv_tcga[surv_tcga$sample %in% colnames(expr_tcga), ]
surv_tcga <- surv_tcga[!duplicated(surv_tcga$sample), ]
surv_tcga <- surv_tcga[!is.na(surv_tcga$OS_time) & surv_tcga$OS_time > 0, ]
message("  Patients with OS data: ", nrow(surv_tcga))

# Use OS as endpoint (TCGA-LIHC clinical doesn't have clean RFS columns)
time_col <- "OS_time"
event_col <- "OS_event"
message("  Using OS (Overall Survival)")

score_signatures <- function(expr, sigs, method = "ssgsea") {
  # Filter genes present in expression matrix
  sigs_filtered <- lapply(sigs, function(g) intersect(g, rownames(expr)))
  
  n_found <- sapply(sigs_filtered, length)
  message("  Genes found in expression: ", paste(names(n_found), n_found, sep = "=", collapse = ", "))
  
  # Remove empty signatures
  sigs_filtered <- sigs_filtered[n_found >= 5]
  
  if (length(sigs_filtered) == 0) {
    warning("No signatures with >= 5 genes found!")
    return(NULL)
  }
  
  # ssGSEA (new GSVA API)
  param <- ssgseaParam(exprData = as.matrix(expr), geneSets = sigs_filtered, normalize = TRUE)
  scores <- gsva(param, verbose = FALSE)
  
  return(scores)
}

# ---- Score TCGA ----
message("  Scoring TCGA-LIHC...")
scores_tcga <- score_signatures(expr_tcga, signatures)

if (!is.null(scores_tcga)) {
  message("  Score dimensions: ", nrow(scores_tcga), " signatures x ", ncol(scores_tcga), " samples")
}

# PART 3: Survival analysis function
message("\n===== PART 3: Survival analysis =====")

run_survival <- function(scores, surv_data, cohort_name, time_col, event_col) {
  
  results <- list()
  plots <- list()
  
  # Match samples
  common <- intersect(colnames(scores), surv_data$sample)
  if (length(common) < 20) {
    message("  WARNING: Only ", length(common), " matched samples in ", cohort_name)
    return(NULL)
  }
  message("  ", cohort_name, ": ", length(common), " patients matched")
  
  surv_sub <- surv_data[match(common, surv_data$sample), ]
  
  for (sig_name in rownames(scores)) {
    sig_scores <- scores[sig_name, common]
    
    # Median split
    median_score <- median(sig_scores, na.rm = TRUE)
    group <- ifelse(sig_scores > median_score, "High", "Low")
    
    df <- data.frame(
      time = as.numeric(surv_sub[[time_col]]),
      event = as.numeric(surv_sub[[event_col]]),
      group = factor(group, levels = c("Low", "High")),
      score = sig_scores,
      stringsAsFactors = FALSE
    )
    
    # Remove NA
    df <- df[complete.cases(df), ]
    df <- df[df$time > 0, ]
    
    if (nrow(df) < 20) next
    
    # Convert days to months
    df$time_months <- df$time / 30.44
    
    # ---- KM ----
    fit <- survfit(Surv(time_months, event) ~ group, data = df)
    
    # Log-rank test
    lr <- survdiff(Surv(time_months, event) ~ group, data = df)
    p_logrank <- 1 - pchisq(lr$chisq, 1)
    
    # ---- Cox (use median-split group for interpretable HR) ----
    cox <- coxph(Surv(time_months, event) ~ group, data = df)
    cox_sum <- summary(cox)
    hr <- cox_sum$conf.int[1, 1]       # HR for High vs Low
    hr_lower <- cox_sum$conf.int[1, 3]
    hr_upper <- cox_sum$conf.int[1, 4]
    p_cox <- cox_sum$coefficients[1, 5]
    
    # ---- Time-dependent ROC (3-year) ----
    auc_3yr <- NA
    tryCatch({
      library(timeROC)
      roc <- timeROC(T = df$time_months, delta = df$event,
                     marker = scale(df$score)[,1], cause = 1,
                     times = c(12, 24, 36),
                     iid = FALSE)
      auc_3yr <- roc$AUC[3]  # 36-month AUC
      # Correct direction: if AUC < 0.5, flip
      if (!is.na(auc_3yr) && auc_3yr < 0.5) auc_3yr <- 1 - auc_3yr
    }, error = function(e) {})
    
    results[[sig_name]] <- data.frame(
      cohort = cohort_name,
      signature = sig_name,
      n = nrow(df),
      n_events = sum(df$event),
      HR = round(hr, 3),
      HR_lower = round(hr_lower, 3),
      HR_upper = round(hr_upper, 3),
      p_cox = signif(p_cox, 3),
      p_logrank = signif(p_logrank, 3),
      AUC_3yr = round(auc_3yr, 3),
      stringsAsFactors = FALSE
    )
    
    # ---- KM plot ----
    p <- ggsurvplot(fit, data = df,
                    title = paste0(cohort_name, ": ", sig_name, " Signature"),
                    subtitle = paste0("HR=", round(hr, 2),
                                     " (", round(hr_lower, 2), "-", round(hr_upper, 2), ")",
                                     ", p=", signif(p_logrank, 3)),
                    xlab = "Time (months)",
                    ylab = ifelse(grepl("DFS|PFI|RFS", time_col),
                                  "Recurrence-Free Survival", "Overall Survival"),
                    palette = c("#2166AC", "#B2182B"),
                    risk.table = TRUE,
                    risk.table.height = 0.25,
                    conf.int = TRUE,
                    pval = FALSE,  # already in subtitle
                    legend.labs = c("Low score", "High score"))
    
    plots[[sig_name]] <- p
  }
  
  return(list(results = do.call(rbind, results), plots = plots))
}

# ---- Run TCGA ----
res_tcga <- run_survival(scores_tcga, surv_tcga, "TCGA-LIHC", time_col, event_col)

if (!is.null(res_tcga)) {
  message("\n  TCGA Results:")
  print(res_tcga$results)
  
  # Save plots
  for (sig_name in names(res_tcga$plots)) {
    pdf(file.path(OUT_DIR, paste0("KM_TCGA_", sig_name, ".pdf")), width = 7, height = 6)
    print(res_tcga$plots[[sig_name]])
    dev.off()
  }
}

# PART 4: GSE14520 (Validation 1)
message("\n===== PART 4: GSE14520 =====")

gse14520_dir <- file.path(OUT_DIR, "GSE14520")
dir.create(gse14520_dir, showWarnings = FALSE)

gse14520_expr_file <- file.path(gse14520_dir, "GSE14520_expr.rds")
gse14520_clin_file <- file.path(gse14520_dir, "GSE14520_clin.rds")

if (!file.exists(gse14520_expr_file)) {
  message("  Downloading GSE14520...")
  
  # Download the GPL3921 series matrix file
  gse_file <- file.path(gse14520_dir, "GSE14520-GPL3921_series_matrix.txt.gz")
  if (!file.exists(gse_file)) {
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz",
      destfile = gse_file, mode = "wb", quiet = FALSE
    )
  }
  
  # Manual parsing - this file has two data tables, we want the second (larger) one
  lines <- readLines(gzfile(gse_file))
  
  # Find all data table boundaries
  begin_idx <- which(grepl("^!series_matrix_table_begin", lines))
  end_idx <- which(grepl("^!series_matrix_table_end", lines))
  message("  Found ", length(begin_idx), " data tables")
  
  # Use the largest table (most rows)
  table_sizes <- end_idx - begin_idx
  best <- which.max(table_sizes)
  message("  Using table ", best, " (", table_sizes[best], " rows)")
  
  # Extract expression data
  data_lines <- lines[(begin_idx[best] + 1):(end_idx[best] - 1)]
  expr_con <- textConnection(data_lines)
  expr_raw <- read.delim(expr_con, row.names = 1, check.names = FALSE)
  close(expr_con)
  expr_gse <- as.matrix(expr_raw)
  message("  Raw expression: ", nrow(expr_gse), " probes x ", ncol(expr_gse), " samples")
  
  # Extract phenotype data from header
  meta_lines <- lines[1:(begin_idx[1] - 1)]
  
  # Parse sample characteristics
  sample_lines <- meta_lines[grepl("^!Sample_", meta_lines)]
  
  # Get sample IDs
  title_line <- meta_lines[grepl("^!Sample_geo_accession", meta_lines)][1]
  sample_ids <- unlist(strsplit(gsub("^!Sample_geo_accession\t", "", title_line), "\t"))
  sample_ids <- gsub('"', '', sample_ids)
  
  # Build phenotype data frame
  pheno <- data.frame(row.names = sample_ids)
  
  char_lines <- meta_lines[grepl("^!Sample_characteristics_ch1", meta_lines)]
  for (cl in char_lines) {
    vals <- unlist(strsplit(gsub("^!Sample_characteristics_ch1\t", "", cl), "\t"))
    vals <- gsub('"', '', vals)
    if (length(vals) != length(sample_ids)) next
    # Extract key:value
    key <- sub(":.*", "", vals[1])
    vals_clean <- sub("^[^:]+: *", "", vals)
    pheno[[key]] <- vals_clean
  }
  
  message("  Pheno columns: ", paste(colnames(pheno), collapse = ", "))
  
  # Download GPL3921 annotation directly
  message("  Downloading GPL3921 annotation...")
  gpl_file <- file.path(gse14520_dir, "GPL3921.annot.gz")
  if (!file.exists(gpl_file)) {
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL3nnn/GPL3921/annot/GPL3921.annot.gz",
      destfile = gpl_file, mode = "wb", quiet = FALSE
    )
  }
  
  # Parse annotation: find table start after !platform_table_begin
  gpl_lines <- readLines(gzfile(gpl_file))
  table_start <- which(grepl("^!platform_table_begin", gpl_lines))
  table_end <- which(grepl("^!platform_table_end", gpl_lines))
  if (length(table_end) == 0) table_end <- length(gpl_lines)
  
  data_lines <- gpl_lines[(table_start + 1):(table_end - 1)]
  
  # Only read first 3 columns: ID, Gene title, Gene symbol
  # Split each line by tab and extract columns 1 and 3
  parsed <- strsplit(data_lines, "\t")
  header <- parsed[[1]]
  message("  GPL header: ", paste(header[1:min(5, length(header))], collapse = " | "))
  
  probe_ids <- sapply(parsed[-1], function(x) if(length(x) >= 1) x[1] else NA)
  gene_syms <- sapply(parsed[-1], function(x) if(length(x) >= 3) x[3] else NA)
  
  probe_genes <- data.frame(ID = probe_ids, Gene_Symbol = gene_syms, stringsAsFactors = FALSE)
  probe_genes <- probe_genes[!is.na(probe_genes$Gene_Symbol) & 
                              probe_genes$Gene_Symbol != "" & 
                              probe_genes$Gene_Symbol != "---", ]
  message("  Probes with gene symbols: ", nrow(probe_genes))
  
  # Match probes
  common_probes <- intersect(rownames(expr_gse), probe_genes$ID)
  expr_gse <- expr_gse[common_probes, ]
  gene_symbols <- probe_genes$Gene_Symbol[match(common_probes, probe_genes$ID)]
  gene_symbols <- sapply(strsplit(gene_symbols, " */// *"), `[`, 1)
  
  # Remove empty
  keep <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  expr_gse <- expr_gse[keep, ]
  gene_symbols <- gene_symbols[keep]
  
  # Average duplicates
  expr_gse <- as.data.frame(expr_gse)
  expr_gse$gene <- gene_symbols
  expr_gse <- expr_gse %>%
    group_by(gene) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()
  rownames(expr_gse) <- expr_gse$gene
  expr_gse$gene <- NULL
  expr_gse <- as.matrix(expr_gse)
  
  message("  Final: ", nrow(expr_gse), " genes x ", ncol(expr_gse), " samples")
  
  saveRDS(expr_gse, gse14520_expr_file)
  saveRDS(pheno, gse14520_clin_file)
} else {
  message("  Loading cached GSE14520...")
  expr_gse <- readRDS(gse14520_expr_file)
  pheno <- readRDS(gse14520_clin_file)
}

message("  GSE14520: ", nrow(expr_gse), " genes x ", ncol(expr_gse), " samples")
message("  Clinical columns: ", paste(head(colnames(pheno), 20), collapse = ", "))

# ---- Extract survival ----
# GSE14520 pheno columns vary; look for recurrence info
rfs_candidates <- grep("recurr|relapse|rfs|dfs|time|survival|status",
                       colnames(pheno), ignore.case = TRUE, value = TRUE)
message("  Survival-related columns: ", paste(rfs_candidates, collapse = ", "))

# Typical GSE14520 columns:
#   "Recurrence Free Survival (Months):ch1" and "Recurrence:ch1"
#   or "Time to recurrence (months):ch1" and "Recurrence status:ch1"
surv_gse <- data.frame(sample = rownames(pheno), stringsAsFactors = FALSE)

# Try to find time and event columns
for (col in colnames(pheno)) {
  col_lower <- tolower(col)
  vals <- pheno[[col]]
  
  if (grepl("recurr.*free.*survival.*month|time.*recurr.*month|rfs.*month", col_lower)) {
    surv_gse$RFS_time <- as.numeric(vals) * 30.44  # convert to days
    message("  RFS time column: ", col)
  }
  if (grepl("recurr.*status|recurrence:ch1$|recurr$", col_lower)) {
    surv_gse$RFS_event <- ifelse(tolower(vals) %in% c("yes", "1", "recurrence", "recurred"), 1, 0)
    message("  RFS event column: ", col)
  }
  if (grepl("overall.*survival.*month|os.*month", col_lower)) {
    surv_gse$OS_time <- as.numeric(vals) * 30.44
    message("  OS time column: ", col)
  }
  if (grepl("survival.*status|vital.*status|os.*status", col_lower)) {
    surv_gse$OS_event <- ifelse(tolower(vals) %in% c("dead", "1", "deceased"), 1, 0)
    message("  OS event column: ", col)
  }
}

# ---- Keep tumor samples only ----
tissue_col <- grep("tissue|type|source", colnames(pheno), ignore.case = TRUE, value = TRUE)
if (length(tissue_col) > 0) {
  for (tc in tissue_col) {
    if (any(grepl("tumor|tumour|hcc|cancer", pheno[[tc]], ignore.case = TRUE))) {
      tumor_mask <- grepl("tumor|tumour|hcc|cancer", pheno[[tc]], ignore.case = TRUE)
      expr_gse <- expr_gse[, tumor_mask]
      surv_gse <- surv_gse[tumor_mask, ]
      message("  Filtered to tumor samples: ", ncol(expr_gse))
      break
    }
  }
}

# ---- Score ----
message("  Scoring GSE14520...")
scores_gse <- score_signatures(expr_gse, signatures)

# ---- Survival ----
if ("RFS_time" %in% names(surv_gse)) {
  res_gse <- run_survival(scores_gse, surv_gse, "GSE14520", "RFS_time", "RFS_event")
} else if ("OS_time" %in% names(surv_gse)) {
  res_gse <- run_survival(scores_gse, surv_gse, "GSE14520", "OS_time", "OS_event")
} else {
  message("  WARNING: No survival data found in GSE14520 phenodata")
  res_gse <- NULL
}

if (!is.null(res_gse)) {
  message("\n  GSE14520 Results:")
  print(res_gse$results)
  
  for (sig_name in names(res_gse$plots)) {
    pdf(file.path(OUT_DIR, paste0("KM_GSE14520_", sig_name, ".pdf")), width = 7, height = 6)
    print(res_gse$plots[[sig_name]])
    dev.off()
  }
}

# PART 5: ICGC-LIRI (Validation 2)
message("\n===== PART 5: ICGC-LIRI =====")

# ICGC data needs to be downloaded from https://dcc.icgc.org/
# If not available locally, try to download
icgc_dir <- file.path(OUT_DIR, "ICGC")
dir.create(icgc_dir, showWarnings = FALSE)

icgc_expr_file <- file.path(icgc_dir, "ICGC_LIRI_expr.rds")

if (!file.exists(icgc_expr_file)) {
  message("  ICGC-LIRI data not found locally")
  message("  Please download from: https://dcc.icgc.org/projects/LIRI-JP")
  message("  Required files: exp_seq.tsv.gz + donor.tsv.gz")
  message("  Skipping ICGC validation...")
  res_icgc <- NULL
} else {
  message("  Loading ICGC-LIRI...")
  expr_icgc <- readRDS(icgc_expr_file)
  clin_icgc <- readRDS(file.path(icgc_dir, "ICGC_LIRI_clin.rds"))
  
  scores_icgc <- score_signatures(expr_icgc, signatures)
  res_icgc <- run_survival(scores_icgc, clin_icgc, "ICGC-LIRI", "OS_time", "OS_event")
}

# PART 6: Combined results + comparative figure
message("\n===== PART 6: Summary =====")

all_results <- rbind(
  if (!is.null(res_tcga)) res_tcga$results,
  if (!is.null(res_gse)) res_gse$results,
  if (!is.null(res_icgc)) res_icgc$results
)

write.csv(all_results, file.path(OUT_DIR, "Table_Recurrence_Prediction.csv"), row.names = FALSE)

message("\n  ====== SUMMARY TABLE ======")
print(all_results)

# ---- Forest plot: HR comparison across signatures and cohorts ----
if (nrow(all_results) >= 3) {
  
  all_results$label <- paste0(all_results$cohort, "\n", all_results$signature)
  all_results$sig_type <- factor(all_results$signature,
                                  levels = c("CSC_Silenced", "CSC_Activated",
                                             "Passive_Hyper", "Passive_Hypo", "Random"))
  
  p_forest <- ggplot(all_results, aes(x = HR, y = reorder(label, HR),
                                       color = sig_type)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("CSC_Silenced"  = "#2166AC",
                                   "CSC_Activated" = "#D62728",
                                   "Passive_Hyper" = "#6BAED6",
                                   "Passive_Hypo"  = "#FC9272",
                                   "Random"        = "#7F7F7F"),
                       name = "Signature") +
    labs(x = "Hazard Ratio (95% CI)", y = "",
         title = "Erosion-Guided Signature Predicts HCC Recurrence",
         subtitle = "Low erosion = CSC actively maintained = strongest predictor") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
  
  # Add significance stars
  p_forest <- p_forest +
    geom_text(aes(x = HR_upper + 0.1,
                  label = ifelse(p_cox < 0.001, "***",
                                 ifelse(p_cox < 0.01, "**",
                                        ifelse(p_cox < 0.05, "*", "ns")))),
              hjust = 0, size = 4, show.legend = FALSE)
  
  ggsave(file.path(OUT_DIR, "Fig_Forest_ErosionVsRecurrence.pdf"),
         p_forest, width = 8, height = max(4, nrow(all_results) * 0.8), dpi = 300)
  message("  -> Fig_Forest_ErosionVsRecurrence.pdf")
}

# ---- Bar plot: AUC comparison ----
auc_data <- all_results[!is.na(all_results$AUC_3yr), ]
if (nrow(auc_data) > 0) {
  p_auc <- ggplot(auc_data, aes(x = signature, y = AUC_3yr, fill = signature)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_text(aes(label = round(AUC_3yr, 3)), vjust = -0.5, size = 3.5) +
    facet_wrap(~cohort) +
    scale_fill_manual(values = c("CSC_Silenced"  = "#2166AC",
                                  "CSC_Activated" = "#D62728",
                                  "Passive_Hyper" = "#6BAED6",
                                  "Passive_Hypo"  = "#FC9272",
                                  "Random"        = "#7F7F7F")) +
    labs(y = "3-Year AUC", x = "", title = "Predictive Accuracy for HCC Recurrence") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    ylim(0, 1)
  
  ggsave(file.path(OUT_DIR, "Fig_AUC_ErosionVsRecurrence.pdf"),
         p_auc, width = 8, height = 5, dpi = 300)
  message("  -> Fig_AUC_ErosionVsRecurrence.pdf")
}

# PART 7: Individual gene survival (top erosion-guided genes)
message("\n===== PART 7: Top gene individual survival =====")

# Top 10 low-erosion genes with RNA data
top_genes <- erosion %>%
  filter(!is.na(mean_erosion), mean_erosion < 0.3,
         abs(mean_baseline) > 10,
         !is.na(rna_padj), rna_padj < 0.05) %>%
  arrange(mean_erosion) %>%
  head(10) %>%
  pull(gene)

message("  Testing individual genes: ", paste(top_genes, collapse = ", "))

# Use TCGA expression
gene_results <- list()

for (g in top_genes) {
  if (!(g %in% rownames(expr_tcga))) next
  
  gexpr <- as.numeric(expr_tcga[g, ])
  names(gexpr) <- colnames(expr_tcga)
  
  common <- intersect(names(gexpr), surv_tcga$sample)
  if (length(common) < 30) next
  
  df <- data.frame(
    time = as.numeric(surv_tcga[[time_col]][match(common, surv_tcga$sample)]),
    event = as.numeric(surv_tcga[[event_col]][match(common, surv_tcga$sample)]),
    expr = gexpr[common],
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df) & df$time > 0, ]
  if (nrow(df) < 30) next
  
  df$time_months <- df$time / 30.44
  
  cox <- coxph(Surv(time_months, event) ~ expr, data = df)
  s <- summary(cox)
  
  # Get erosion info
  erosion_val <- erosion$mean_erosion[erosion$gene == g]
  baseline_val <- erosion$mean_baseline[erosion$gene == g]
  
  gene_results[[g]] <- data.frame(
    gene = g,
    erosion = round(erosion_val, 3),
    baseline_diff = round(baseline_val, 1),
    HR = round(s$conf.int[1, 1], 3),
    p = signif(s$coefficients[1, 5], 3),
    stringsAsFactors = FALSE
  )
}

if (length(gene_results) > 0) {
  gene_df <- do.call(rbind, gene_results)
  gene_df <- gene_df[order(gene_df$p), ]
  write.csv(gene_df, file.path(OUT_DIR, "Table_TopGenes_IndividualSurvival.csv"), row.names = FALSE)
  
  message("\n  Individual gene survival (TCGA RFS):")
  print(gene_df)
}

message("\n", paste(rep("=", 60), collapse = ""))
message("  ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("\n  Output: ", OUT_DIR)
message("\n  Key question answered:")
message("  If Low_Erosion HR > 1 and p < 0.05,")
message("    AND High_Erosion / Random are NOT significant,")
