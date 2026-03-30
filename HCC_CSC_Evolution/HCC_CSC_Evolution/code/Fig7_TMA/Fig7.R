# NOTE: Chinese column names (点位置, 组织类型, 生存状态, 生存期, 癌/癌旁) match
# the clinical data Excel file and cannot be renamed without modifying input data.
#!/usr/bin/env Rscript
# Figure 7 — Tissue microarray immunofluorescence
# 7h,i: 5hmC positivity rate comparisons
# 7j: KM survival (CD13+CD133+ 5hmC)
# S8: Adjacent tissue IF + prognostic significance

# --- 7h,i,j + S8k: MIF quantification + KM survival ---

##############################################################################
##############################################################################

# ---- 0.  ----
suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(survival)
  library(survminer)
  library(dplyr)
})

load_and_merge <- function(file, clin_file = "clin.xlsx") {
  dat <- read_excel(file)
  
  b <- paste(str_extract(dat$`Annotation ID`, "(?<=,)[^,]+"))
  a <- str_match(dat$`Annotation ID`, ",[^,]*,([^\\]]*)\\]")[, 2]
  ID <- character(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    if (nchar(b[i]) == 1) {
      ID[i] <- paste(a[i], b[i], sep = "0")
    } else {
      ID[i] <- paste0(a[i], b[i])
    }
  }
  dat$ID <- ID
  
  clin <- read_excel(clin_file)
  clintumor <- clin
  clinpara  <- clin
  clintumor$点位置1 <- str_sub(clintumor$点位置, 1, 3)
  clintumor$组织类型 <- str_sub(clintumor$组织类型, 1, 1)
  clinpara$点位置1  <- str_sub(clinpara$点位置, 5, 7)
  clinpara <- clinpara[clinpara$点位置1 != "", ]
  clinpara$组织类型 <- str_sub(clinpara$组织类型, 3, 4)
  clin_all <- rbind(clintumor, clinpara)
  colnames(clin_all)[22] <- "ID"
  
  merged <- merge(dat, clin_all, by = "ID")
  
  if (all(c("EpCAM+/5hmC+", "EpCAM+") %in% colnames(merged))) {
    merged$EpCAM_5hmc <- merged$`EpCAM+/5hmC+` / merged$`EpCAM+`
  }
  if (all(c("CD133+/5hmC+", "CD133+") %in% colnames(merged))) {
    merged$cd133_5hmc <- merged$`CD133+/5hmC+` / merged$`CD133+`
  }
  if (all(c("CD13+/CD133+/5hmC+", "CD13+/CD133+") %in% colnames(merged))) {
    merged$cd13_133_5hmc <- merged$`CD13+/CD133+/5hmC+` / merged$`CD13+/CD133+`
  }
  if (all(c("CD133+/EpCAM+", "CD133+") %in% colnames(merged))) {
    merged$cd133_EpCAM <- merged$`CD133+/EpCAM+` / merged$`CD133+`
  }
  if (all(c("EpCAM+/5hmC+", "CD133+") %in% colnames(merged))) {
    merged$EpCAM5hmc_cd133 <- merged$`EpCAM+/5hmC+` / merged$`CD133+`
  }
  single_markers <- intersect(
    c("EpCAM+", "CD133+", "5hmC+", "CD13+/CD133+", 
      "EpCAM+/5hmC+", "CD133+/5hmC+", "CD13+/CD133+/5hmC+"),
    colnames(merged)
  )
  
  merged$生存状态 <- as.numeric(str_sub(merged$生存状态, 1, 1))
  
  return(merged)
}

# ---- 2. cutoff ----
get_cutoffs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  
  cuts <- list()
  cuts[["median"]]  <- median(x, na.rm = TRUE)
  cuts[["mean"]]    <- mean(x, na.rm = TRUE)
  cuts[["Q1"]]      <- quantile(x, 0.25, na.rm = TRUE)
  cuts[["Q3"]]      <- quantile(x, 0.75, na.rm = TRUE)
  cuts[["tertile_low"]]  <- quantile(x, 1/3, na.rm = TRUE)
  cuts[["tertile_high"]] <- quantile(x, 2/3, na.rm = TRUE)
  
  # surv_cutpoint (maxstat )
  return(cuts)
}

run_single_km <- function(data, marker, cutoff_value, cutoff_name,
                          surv_time_col = "生存期2019/12",
                          surv_status_col = "生存状态") {
  
  df <- data[, c(surv_time_col, surv_status_col, marker), drop = FALSE]
  colnames(df) <- c("time", "status", "marker")
  df <- df[is.finite(df$marker) & is.finite(df$time) & is.finite(df$status), ]
  
  if (nrow(df) < 10) return(NULL)
  
  df$group <- ifelse(df$marker > cutoff_value, "High", "Low")
  
  tbl <- table(df$group)
  if (length(tbl) < 2 || any(tbl < 3)) return(NULL)
  
  km_fit  <- survfit(Surv(time, status) ~ group, data = df)
  
  # log-rank
  logrank <- survdiff(Surv(time, status) ~ group, data = df)
  pval    <- 1 - pchisq(logrank$chisq, df = 1)
  
  return(list(
    km_fit  = km_fit,
    data    = df,
    pval    = pval,
    n_high  = tbl["High"],
    n_low   = tbl["Low"],
    cutoff  = cutoff_value
  ))
}

# ---- 4. Maxstat  ----
run_maxstat_km <- function(data, marker,
                           surv_time_col = "生存期2019/12",
                           surv_status_col = "生存状态") {
  df <- data[, c(surv_time_col, surv_status_col, marker), drop = FALSE]
  colnames(df) <- c("time", "status", "marker")
  df <- df[is.finite(df$marker) & is.finite(df$time) & is.finite(df$status), ]
  if (nrow(df) < 10) return(NULL)
  
  tryCatch({
    cut_res <- surv_cutpoint(df, time = "time", event = "status", variables = "marker")
    cp <- cut_res$cutpoint$cutpoint
    run_single_km(data, marker, cp, "maxstat", surv_time_col, surv_status_col)
  }, error = function(e) NULL)
}

# ---- 5.  ----
run_all_prognosis <- function(data_files = c("cellcount.xlsx", "celldensity.xlsx", "cellpercent.xlsx"),
                              clin_file = "clin.xlsx",
                              p_threshold = 0.05,
                              output_pdf = "significant_KM_curves.pdf") {
  
  ratio_markers <- c("EpCAM_5hmc", "cd133_5hmc", "cd13_133_5hmc", 
                     "cd133_EpCAM", "EpCAM5hmc_cd133")
  single_markers <- c("EpCAM+", "CD133+", "5hmC+", "CD13+/CD133+",
                      "EpCAM+/5hmC+", "CD133+/5hmC+", "CD13+/CD133+/5hmC+")
  all_markers <- c(ratio_markers, single_markers)
  
  tissue_types <- c("癌", "癌旁")  # "Cancer", "Para-cancer" (matches Excel column values)
  
  cutoff_methods <- c("maxstat")
  
  results <- list()
  sig_results <- list()
  idx <- 0
  
  cat("========================================================\n")
  cat("========================================================\n\n")
  
  for (file in data_files) {
    
    if (!file.exists(file)) {
      next
    }
    
    data_name <- gsub("\\.xlsx$", "", file)
    
    merged <- tryCatch(load_and_merge(file, clin_file), error = function(e) {
      return(NULL)
    })
    if (is.null(merged)) next
    
    for (tissue in tissue_types) {
      
      if (tissue == "癌") {
        sub_data <- merged[merged$组织类型 == "癌", ]
      } else {
        sub_data <- merged[merged$组织类型 != "癌", ]
      }
      
      if (nrow(sub_data) < 10) next
      
      available_markers <- intersect(all_markers, colnames(sub_data))
      available_markers <- available_markers[
        sapply(available_markers, function(m) sum(is.finite(sub_data[[m]])) >= 10)
      ]
      
      for (marker in available_markers) {
        
        res <- run_maxstat_km(sub_data, marker)
        if (is.null(res)) next
        
        idx <- idx + 1
        
        entry <- list(
          dataset  = data_name,
          tissue   = tissue,
          marker   = marker,
          method   = "maxstat",
          cutoff   = res$cutoff,
          pval     = res$pval,
          n_high   = res$n_high,
          n_low    = res$n_low,
          km_fit   = res$km_fit,
          data     = res$data
        )
        
        results[[idx]] <- entry
        
        if (res$pval < p_threshold) {
          sig_results[[length(sig_results) + 1]] <- entry
        }
      }
    }
  }
  
  # ---- 6.  ----
  cat("\n========================================================\n")
  cat("========================================================\n\n")
  
  if (length(sig_results) == 0) {
    return(invisible(NULL))
  }
  
  pvals <- sapply(sig_results, function(x) x$pval)
  sig_results <- sig_results[order(pvals)]
  
  cat(sprintf("%-15s %-6s %-20s %-12s %-10s %-8s %-8s %s\n",
              "数据集", "组织", "标志物", "Cutoff方法", "Cutoff值",
              "N_High", "N_Low", "P值"))
  cat(paste(rep("-", 100), collapse = ""), "\n")
  
  for (r in sig_results) {
    stars <- ifelse(r$pval < 0.001, "***",
                    ifelse(r$pval < 0.01,  "**",
                           ifelse(r$pval < 0.05,  "*", "")))
    cat(sprintf("%-15s %-6s %-20s %-12s %-10.4f %-8d %-8d %.4f %s\n",
                r$dataset, r$tissue, r$marker, r$method,
                r$cutoff, r$n_high, r$n_low, r$pval, stars))
  }
  
  pdf(output_pdf, width = 8, height = 7)
  
  for (i in seq_along(sig_results)) {
    r <- sig_results[[i]]
    
    cox_fit <- coxph(Surv(time, status) ~ group, data = r$data)
    hr <- round(exp(coef(cox_fit)), 2)
    pval_fmt <- formatC(r$pval, format = "g", digits = 3)
    hr_label <- paste0("italic(HR)==", hr, "~italic(P)==", pval_fmt)
    
    tissue_label <- ifelse(r$tissue == "癌", "Tumor", "Para-tumor")
    title_str <- sprintf("%s [%s] %s (cutoff=%.4f)",
                         tissue_label, r$dataset, r$marker, r$cutoff)
    
    p <- ggsurvplot(r$km_fit, data = r$data,
                    pval = FALSE,
                    conf.int = TRUE,
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    risk.table.y.text = FALSE,
                    surv.median.line = "hv",
                    palette = c("#E5A727", "#4A90D9"),
                    legend.labs = c("High", "Low"),
                    legend.title = "Strata",
                    xlab = "Time (months)",
                    ylab = "Survival Probability",
                    title = title_str,
                    ggtheme = theme_bw() + theme(
                      panel.grid = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
                      legend.position = c(0.85, 0.9)
                    ))
    
    p$plot <- p$plot + 
      annotate("text", x = 0, y = 0.15, label = hr_label,
               parse = TRUE, hjust = 0, size = 4.5, color = "#333333")
    
    print(p)
  }
  
  dev.off()
  
  # ---- 8. data.frame ----
  result_df <- do.call(rbind, lapply(sig_results, function(r) {
    data.frame(
      数据集     = r$dataset,
      组织类型   = r$tissue,
      标志物     = r$marker,
      Cutoff方法 = r$method,
      Cutoff值   = round(r$cutoff, 6),
      N_High     = r$n_high,
      N_Low      = r$n_low,
      P值        = signif(r$pval, 4),
      stringsAsFactors = FALSE
    )
  }))
  
  return(result_df)
}

sig_table <- run_all_prognosis(
  data_files = c("cellcount.xlsx", "celldensity.xlsx", "cellpercent.xlsx", "MFI.xlsx"),
  clin_file  = "clin.xlsx",
  p_threshold = 0.05,
  output_pdf  = "significant_KM_curves.pdf"
)

print(sig_table)

##############################################################################
# Part A: cellcountratiocount/density/percent
# {Marker+, Marker-} × {Tumor, Adjacent}3panel
#
# {Marker+, Marker-} × {Tumor, Adjacent}3panel
#
##############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
})

process_one_dataset <- function(file, clin_file = "clin.xlsx") {
  
  dat <- read_excel(file)
  
  b <- paste(str_extract(dat$`Annotation ID`, "(?<=,)[^,]+"))
  a <- str_match(dat$`Annotation ID`, ",[^,]*,([^\\]]*)\\]")[, 2]
  ID <- character(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    if (nchar(b[i]) == 1) {
      ID[i] <- paste(a[i], b[i], sep = "0")
    } else {
      ID[i] <- paste0(a[i], b[i])
    }
  }
  dat$ID <- ID
  
  clin <- read_excel(clin_file)
  clintumor <- clin
  clinpara  <- clin
  clintumor$点位置1 <- str_sub(clintumor$点位置, 1, 3)
  clintumor$组织类型 <- str_sub(clintumor$组织类型, 1, 1)
  clinpara$点位置1  <- str_sub(clinpara$点位置, 5, 7)
  clinpara <- clinpara[clinpara$点位置1 != "", ]
  clinpara$组织类型 <- str_sub(clinpara$组织类型, 3, 4)
  clin_all <- rbind(clintumor, clinpara)
  colnames(clin_all)[22] <- "ID"
  
  merged <- merge(dat, clin_all, by = "ID")
  
  d <- merged[nchar(merged$点位置) == 7, ]
  d <- d[!d$点位置 %in% c("A09/A10", "B17/B18", "B01/B02", "A03/A04", "E05/E06"), ]
  d <- d[d$点位置 %in% d$点位置[duplicated(d$点位置) | duplicated(d$点位置, fromLast = TRUE)], ]
  
  return(d)
}

plot_four_box <- function(df_long, title_str, ylab_str) {
  
  if (is.null(df_long) || nrow(df_long) < 10) return(NULL)
  
  grp_levels <- levels(df_long$group)
  if (length(grp_levels) < 2) return(NULL)
  comparisons <- combn(grp_levels, 2, simplify = FALSE)
  
  fill_colors <- c("#8FBBD9", "#A8D5A2", "#D9C0A0", "#F2C97D")
  
  p <- ggplot(df_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2, color = "grey30") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01,
                       step.increase = 0.08, size = 3.5) +
    labs(title = title_str, x = NULL, y = ylab_str) +
    scale_fill_manual(values = fill_colors[1:length(grp_levels)]) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 13))
  
  return(p)
}

run_ratio_four_groups <- function(file = "cellcount.xlsx",
                                  clin_file = "clin.xlsx",
                                  output_pdf = "ratio_four_group_boxplots.pdf") {
  
  cat("========================================================\n")
  cat("========================================================\n\n")
  
  paired <- process_one_dataset(file, clin_file)
  
  safe_div <- function(a, b) ifelse(b <= 0, NA_real_, a / b)
  
  paired$EpCAM_pos_rate     <- safe_div(paired$`EpCAM+/5hmC+`, paired$`EpCAM+`)
  paired$EpCAM_neg_rate     <- safe_div(paired$`5hmC+` - paired$`EpCAM+/5hmC+`,
                                        paired$`Total Cells` - paired$`EpCAM+`)
  paired$CD133_pos_rate     <- safe_div(paired$`CD133+/5hmC+`, paired$`CD133+`)
  paired$CD133_neg_rate     <- safe_div(paired$`5hmC+` - paired$`CD133+/5hmC+`,
                                        paired$`Total Cells` - paired$`CD133+`)
  paired$CD13CD133_pos_rate <- safe_div(paired$`CD13+/CD133+/5hmC+`, paired$`CD13+/CD133+`)
  paired$CD13CD133_neg_rate <- safe_div(
    paired$`5hmC+` - paired$`CD13+/CD133+/5hmC+`,
    paired$`Total Cells` - paired$`CD13+` - paired$`CD133+` + paired$`CD13+/CD133+`
  )
  
  panels <- list(
    list(pos = "EpCAM_pos_rate",     neg = "EpCAM_neg_rate",
         pos_lab = "EpCAM+",         neg_lab = "EpCAM-",
         title = "EpCAM 5hmC rate"),
    list(pos = "CD133_pos_rate",     neg = "CD133_neg_rate",
         pos_lab = "CD133+",         neg_lab = "CD133-",
         title = "CD133 5hmC rate"),
    list(pos = "CD13CD133_pos_rate", neg = "CD13CD133_neg_rate",
         pos_lab = "CD13+CD133+",    neg_lab = "CD13-CD133-",
         title = "CD13/CD133 5hmC rate")
  )
  
  all_plots <- list()
  
  for (pdef in panels) {
    tissue_en <- ifelse(paired$组织类型 == "癌", "Tumor", "Adjacent")
    
    df_pos <- data.frame(点位置 = paired$点位置,
                         group = paste0(pdef$pos_lab, "\n", tissue_en),
                         value = as.numeric(paired[[pdef$pos]]))
    df_neg <- data.frame(点位置 = paired$点位置,
                         group = paste0(pdef$neg_lab, "\n", tissue_en),
                         value = as.numeric(paired[[pdef$neg]]))
    df_long <- rbind(df_pos, df_neg)
    df_long <- df_long[!is.na(df_long$value) & is.finite(df_long$value), ]
    
    lvls <- c(paste0(pdef$pos_lab, "\nTumor"),  paste0(pdef$pos_lab, "\nAdjacent"),
              paste0(pdef$neg_lab, "\nTumor"),   paste0(pdef$neg_lab, "\nAdjacent"))
    lvls <- intersect(lvls, unique(df_long$group))
    df_long$group <- factor(df_long$group, levels = lvls)
    
    if (pdef$title == "CD13/CD133 5hmC rate") {
      Q1 <- quantile(df_long$value, 0.25, na.rm = TRUE)
      Q3 <- quantile(df_long$value, 0.75, na.rm = TRUE)
      upper <- Q3 + 1.5 * (Q3 - Q1)
      n_before <- nrow(df_long)
      df_long <- df_long[df_long$value <= upper, ]
    }
    
    p <- plot_four_box(df_long, pdef$title, "5hmC positive rate")
    if (!is.null(p)) {
      all_plots[[length(all_plots) + 1]] <- p
    }
  }
  
  pdf(output_pdf, width = 7, height = 6)
  for (p in all_plots) print(p)
  dev.off()
}

run_mfi_four_groups <- function(file = "MFI.xlsx",
                                clin_file = "clin.xlsx",
                                output_pdf = "MFI_four_group_boxplots.pdf") {
  
  cat("========================================================\n")
  cat("========================================================\n\n")
  
  paired <- process_one_dataset(file, clin_file)
  
  # Marker+MFI = EpCAM+/5hmC+
  # Marker-MFI = 5hmC+ - EpCAM+/5hmC+5hmC MFI
  
  panels <- list(
    list(pos_col = "EpCAM+/5hmC+",         neg_col = "5hmC+",
         pos_lab = "EpCAM+",               neg_lab = "EpCAM-",
         neg_sub = "EpCAM+/5hmC+",
         title   = "EpCAM 5hmC MFI"),
    list(pos_col = "CD133+/5hmC+",         neg_col = "5hmC+",
         pos_lab = "CD133+",               neg_lab = "CD133-",
         neg_sub = "CD133+/5hmC+",
         title   = "CD133 5hmC MFI"),
    list(pos_col = "CD13+/CD133+/5hmC+",   neg_col = "5hmC+",
         pos_lab = "CD13+CD133+",          neg_lab = "CD13-CD133-",
         neg_sub = "CD13+/CD133+/5hmC+",
         title   = "CD13/CD133 5hmC MFI")
  )
  
  all_plots <- list()
  
  for (pdef in panels) {
    tissue_en <- ifelse(paired$组织类型 == "癌", "Tumor", "Adjacent")
    
    pos_val <- as.numeric(paired[[pdef$pos_col]])
    neg_col_name <- gsub("\\+", "-", pdef$pos_lab)  # e.g. "EpCAM-"
    
    df_pos <- data.frame(点位置 = paired$点位置,
                         group = paste0(pdef$pos_lab, "\n", tissue_en),
                         value = pos_val)
    
    if (neg_col_name %in% colnames(paired)) {
      neg_val <- as.numeric(paired[[neg_col_name]])
    } else {
      neg_val <- as.numeric(paired[["5hmC+"]]) - pos_val
    }
    
    df_neg <- data.frame(点位置 = paired$点位置,
                         group = paste0(pdef$neg_lab, "\n", tissue_en),
                         value = neg_val)
    
    df_long <- rbind(df_pos, df_neg)
    df_long <- df_long[!is.na(df_long$value) & is.finite(df_long$value), ]
    
    lvls <- c(paste0(pdef$pos_lab, "\nTumor"),  paste0(pdef$pos_lab, "\nAdjacent"),
              paste0(pdef$neg_lab, "\nTumor"),   paste0(pdef$neg_lab, "\nAdjacent"))
    lvls <- intersect(lvls, unique(df_long$group))
    df_long$group <- factor(df_long$group, levels = lvls)
    
    if (grepl("CD13", pdef$title)) {
      Q1 <- quantile(df_long$value, 0.25, na.rm = TRUE)
      Q3 <- quantile(df_long$value, 0.75, na.rm = TRUE)
      upper <- Q3 + 1.5 * (Q3 - Q1)
      n_before <- nrow(df_long)
      df_long <- df_long[df_long$value <= upper, ]
    }
    
    p <- plot_four_box(df_long, pdef$title, "5hmC MFI")
    if (!is.null(p)) {
      all_plots[[length(all_plots) + 1]] <- p
    }
  }
  
  pdf(output_pdf, width = 7, height = 6)
  for (p in all_plots) print(p)
  dev.off()
}

run_mfi_paired <- function(file = "MFI.xlsx",
                           clin_file = "clin.xlsx",
                           output_pdf = "MFI_paired_boxplots.pdf") {
  
  cat("========================================================\n")
  cat("========================================================\n\n")
  
  paired <- process_one_dataset(file, clin_file)
  
  mfi_cols <- c("EpCAM+/5hmC+", "EpCAM+", "CD13+/CD133+", "CD13+/CD133+/5hmC+",
                "CD133+/5hmC+", "CD133+", "EpCAM+/CD133+", "5hmC+", "CD13+", "Total Cells")
  mfi_cols <- intersect(mfi_cols, colnames(paired))
  
  all_plots <- list()
  
  for (col in mfi_cols) {
    
    df <- data.frame(
      tissue_label = ifelse(paired$组织类型 == "癌", "Tumor", "Adjacent"),
      点位置       = paired$点位置,
      value        = as.numeric(paired[[col]])
    )
    df$tissue_label <- factor(df$tissue_label, levels = c("Tumor", "Adjacent"))
    df <- df[!is.na(df$value) & is.finite(df$value), ]
    if (nrow(df) < 6 || length(unique(df$tissue_label)) < 2) next
    
    pval <- tryCatch({
      res <- df %>% wilcox_test(value ~ tissue_label, paired = TRUE)
      res$p
    }, error = function(e) NA)
    
    stars <- ifelse(is.na(pval), "?",
                    ifelse(pval < 0.001, "***",
                           ifelse(pval < 0.01,  "**",
                                  ifelse(pval < 0.05,  "*", "ns"))))
    cat(sprintf("  %-25s p=%.4f %s\n", col, ifelse(is.na(pval), 999, pval), stars))
    
    p <- ggplot(df, aes(x = tissue_label, y = value, fill = tissue_label)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
      stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.signif") +
      labs(title = sprintf("MFI — %s", col), x = NULL, y = "5hmC MFI") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA)) +
      scale_fill_manual(values = c("#00AFBB", "#E7B800"))
    
    all_plots[[length(all_plots) + 1]] <- p
  }
  
  if (length(all_plots) > 0) {
    pdf(output_pdf, width = 6, height = 5.5)
    for (p in all_plots) print(p)
    dev.off()
  }
}

run_ratio_four_groups("cellcount.xlsx", "clin.xlsx", "ratio_four_group_boxplots.pdf")

run_mfi_four_groups("MFI.xlsx", "clin.xlsx", "MFI_four_group_boxplots.pdf")

run_mfi_paired("MFI.xlsx", "clin.xlsx", "MFI_paired_boxplots.pdf")

# --- 7j + S8k: Survival analysis for MIF markers ---

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(survival)
  library(survminer)
  library(dplyr)
})

load_and_merge <- function(file, clin_file = "clin.xlsx") {
  dat <- read_excel(file)
  b <- paste(str_extract(dat$`Annotation ID`, "(?<=,)[^,]+"))
  a <- str_match(dat$`Annotation ID`, ",[^,]*,([^\\]]*)\\]")[, 2]
  ID <- character(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    if (nchar(b[i]) == 1) {
      ID[i] <- paste(a[i], b[i], sep = "0")
    } else {
      ID[i] <- paste0(a[i], b[i])
    }
  }
  dat$ID <- ID
  clin <- read_excel(clin_file)
  clintumor <- clin
  clinpara <- clin
  clintumor$dw1 <- str_sub(clintumor[["\u70b9\u4f4d\u7f6e"]], 1, 3)
  clintumor[["\u7ec4\u7ec7\u7c7b\u578b"]] <- str_sub(clintumor[["\u7ec4\u7ec7\u7c7b\u578b"]], 1, 1)
  clinpara$dw1 <- str_sub(clinpara[["\u70b9\u4f4d\u7f6e"]], 5, 7)
  clinpara <- clinpara[clinpara$dw1 != "", ]
  clinpara[["\u7ec4\u7ec7\u7c7b\u578b"]] <- str_sub(clinpara[["\u7ec4\u7ec7\u7c7b\u578b"]], 3, 4)
  clin_all <- rbind(clintumor, clinpara)
  colnames(clin_all)[22] <- "ID"
  merged <- merge(dat, clin_all, by = "ID")
  
  safe_div <- function(a, b) ifelse(b <= 0, NA_real_, a / b)
  
  num_cols <- c("EpCAM+/5hmC+", "EpCAM+", "CD13+/CD133+", "CD13+/CD133+/5hmC+",
                "CD133+/5hmC+", "CD133+", "5hmC+", "CD13+", "Total Cells")
  for (nc in intersect(num_cols, colnames(merged))) {
    if (is.character(merged[[nc]])) {
      merged[[nc]] <- as.numeric(gsub("%", "", merged[[nc]])) / 100
    }
  }
  
  if (all(c("EpCAM+/5hmC+", "EpCAM+") %in% colnames(merged)))
    merged$EpCAM_5hmc <- safe_div(merged$`EpCAM+/5hmC+`, merged$`EpCAM+`)
  if (all(c("CD133+/5hmC+", "CD133+") %in% colnames(merged)))
    merged$cd133_5hmc <- safe_div(merged$`CD133+/5hmC+`, merged$`CD133+`)
  if (all(c("CD13+/CD133+/5hmC+", "CD13+/CD133+") %in% colnames(merged)))
    merged$cd13_133_5hmc <- safe_div(merged$`CD13+/CD133+/5hmC+`, merged$`CD13+/CD133+`)
  
  merged[["\u751f\u5b58\u72b6\u6001"]] <- as.numeric(str_sub(merged[["\u751f\u5b58\u72b6\u6001"]], 1, 1))
  return(merged)
}

run_single_km <- function(data, marker, cutoff_value,
                          surv_time_col = "\u751f\u5b58\u671f2019/12",
                          surv_status_col = "\u751f\u5b58\u72b6\u6001") {
  df <- data[, c(surv_time_col, surv_status_col, marker), drop = FALSE]
  colnames(df) <- c("time", "status", "marker")
  df <- df[is.finite(df$marker) & is.finite(df$time) & is.finite(df$status), ]
  if (nrow(df) < 10) return(NULL)
  df$group <- ifelse(df$marker > cutoff_value, "High", "Low")
  tbl <- table(df$group)
  if (length(tbl) < 2 || any(tbl < 3)) return(NULL)
  km_fit <- survfit(Surv(time, status) ~ group, data = df)
  logrank <- survdiff(Surv(time, status) ~ group, data = df)
  pval <- 1 - pchisq(logrank$chisq, df = 1)
  return(list(km_fit = km_fit, data = df, pval = pval,
              n_high = tbl["High"], n_low = tbl["Low"], cutoff = cutoff_value))
}

run_maxstat_km <- function(data, marker,
                           surv_time_col = "\u751f\u5b58\u671f2019/12",
                           surv_status_col = "\u751f\u5b58\u72b6\u6001") {
  df <- data[, c(surv_time_col, surv_status_col, marker), drop = FALSE]
  colnames(df) <- c("time", "status", "marker")
  df <- df[is.finite(df$marker) & is.finite(df$time) & is.finite(df$status), ]
  if (nrow(df) < 10) return(NULL)
  tryCatch({
    cut_res <- surv_cutpoint(df, time = "time", event = "status", variables = "marker")
    cp <- cut_res$cutpoint$cutpoint
    run_single_km(data, marker, cp, surv_time_col, surv_status_col)
  }, error = function(e) NULL)
}

run_all_prognosis <- function(
    ratio_files = c("cellcount.xlsx", "celldensity.xlsx", "cellpercent.xlsx"),
    mfi_file = "MFI.xlsx",
    clin_file = "clin.xlsx",
    p_threshold = 0.05,
    output_pdf = "significant_KM_curves.pdf") {
  
  ratio_markers <- c("EpCAM_5hmc", "cd133_5hmc", "cd13_133_5hmc")
  mfi_markers <- c("EpCAM+", "CD133+", "CD13+/CD133+")
  
  tissue_types <- list(
    list(label = "Tumor", filter_val = "\u764c", negate = FALSE),
    list(label = "Para-tumor", filter_val = "\u764c", negate = TRUE)
  )
  
  zz_col <- "\u7ec4\u7ec7\u7c7b\u578b"
  results <- list()
  sig_results <- list()
  idx <- 0
  
  cat("========================================================\n")
  cat("  Prognosis: ratio + MFI\n")
  cat("  Ratio: EpCAM_5hmc, cd133_5hmc, cd13_133_5hmc\n")
  cat("  MFI:   EpCAM+, CD133+, CD13+/CD133+\n")
  cat("========================================================\n\n")
  
  # --- ratio files ---
  for (file in ratio_files) {
    if (!file.exists(file)) next
    data_name <- gsub("\\.xlsx$", "", file)
    cat(sprintf(">>> %s (ratio)\n", data_name))
    
    merged <- tryCatch(load_and_merge(file, clin_file), error = function(e) {
      cat(sprintf("  [error] %s\n", e$message)); NULL
    })
    if (is.null(merged)) next
    
    for (tt in tissue_types) {
      if (isTRUE(tt$negate)) {
        sub_data <- merged[merged[[zz_col]] != tt$filter_val, ]
      } else {
        sub_data <- merged[merged[[zz_col]] == tt$filter_val, ]
      }
      if (nrow(sub_data) < 10) next
      
      for (marker in ratio_markers) {
        if (!marker %in% colnames(sub_data)) next
        n_valid <- sum(is.finite(sub_data[[marker]]))
        cat(sprintf("  %s | %s | %s: n_valid=%d ", data_name, tt$label, marker, n_valid))
        if (n_valid < 10) { cat("skip\n"); next }
        
        res <- run_maxstat_km(sub_data, marker)
        if (is.null(res)) { cat("maxstat failed\n"); next }
        cat(sprintf("p=%.4f\n", res$pval))
        
        idx <- idx + 1
        entry <- list(dataset = data_name, tissue = tt$label, marker = marker,
                      cutoff = res$cutoff, pval = res$pval, n_high = res$n_high,
                      n_low = res$n_low, km_fit = res$km_fit, data = res$data)
        results[[idx]] <- entry
        if (res$pval < p_threshold) sig_results[[length(sig_results) + 1]] <- entry
      }
    }
  }
  
  # --- MFI file ---
  if (file.exists(mfi_file)) {
    cat(sprintf("\n>>> MFI\n"))
    merged <- tryCatch(load_and_merge(mfi_file, clin_file), error = function(e) {
      cat(sprintf("  [error] %s\n", e$message)); NULL
    })
    
    if (!is.null(merged)) {
      for (tt in tissue_types) {
        if (isTRUE(tt$negate)) {
          sub_data <- merged[merged[[zz_col]] != tt$filter_val, ]
        } else {
          sub_data <- merged[merged[[zz_col]] == tt$filter_val, ]
        }
        if (nrow(sub_data) < 10) next
        
        for (marker in mfi_markers) {
          if (!marker %in% colnames(sub_data)) next
          sub_data[[marker]] <- as.numeric(sub_data[[marker]])
          n_valid <- sum(is.finite(sub_data[[marker]]))
          cat(sprintf("  MFI | %s | %s: n_valid=%d ", tt$label, marker, n_valid))
          if (n_valid < 10) { cat("skip\n"); next }
          
          res <- run_maxstat_km(sub_data, marker)
          if (is.null(res)) { cat("maxstat failed\n"); next }
          cat(sprintf("p=%.4f\n", res$pval))
          
          idx <- idx + 1
          entry <- list(dataset = "MFI", tissue = tt$label, marker = marker,
                        cutoff = res$cutoff, pval = res$pval, n_high = res$n_high,
                        n_low = res$n_low, km_fit = res$km_fit, data = res$data)
          results[[idx]] <- entry
          if (res$pval < p_threshold) sig_results[[length(sig_results) + 1]] <- entry
        }
      }
    }
  }
  
  # --- summary ---
  cat("\n========================================================\n")
  cat(sprintf("  Total: %d combinations\n", idx))
  cat(sprintf("  Significant (p < %.2f): %d\n", p_threshold, length(sig_results)))
  cat("========================================================\n\n")
  
  if (length(sig_results) == 0) {
    cat("No significant results.\n")
    return(invisible(NULL))
  }
  
  pvals <- sapply(sig_results, function(x) x$pval)
  sig_results <- sig_results[order(pvals)]
  
  cat(sprintf("%-12s %-12s %-18s %-10s %-8s %-8s %s\n",
              "Dataset", "Tissue", "Marker", "Cutoff", "N_High", "N_Low", "P"))
  cat(paste(rep("-", 80), collapse = ""), "\n")
  for (r in sig_results) {
    stars <- ifelse(r$pval < 0.001, "***", ifelse(r$pval < 0.01, "**", ifelse(r$pval < 0.05, "*", "")))
    cat(sprintf("%-12s %-12s %-18s %-10.4f %-8d %-8d %.4f %s\n",
                r$dataset, r$tissue, r$marker, r$cutoff, r$n_high, r$n_low, r$pval, stars))
  }
  
  # --- KM plots: save ALL results ---
  all_sorted <- results[order(sapply(results, function(x) x$pval))]
  cat(sprintf("\n>>> Saving ALL %d KM curves: %s\n", length(all_sorted), output_pdf))
  pdf(output_pdf, width = 8, height = 7)
  for (r in all_sorted) {
    cox_fit <- coxph(Surv(time, status) ~ group, data = r$data)
    hr <- round(exp(coef(cox_fit)), 2)
    pval_fmt <- formatC(r$pval, format = "g", digits = 3)
    hr_label <- paste0("italic(HR)==", hr, "~italic(P)==", pval_fmt)
    title_str <- sprintf("%s [%s] %s (cutoff=%.4f)", r$tissue, r$dataset, r$marker, r$cutoff)
    
    p <- ggsurvplot(r$km_fit, data = r$data,
                    pval = FALSE, conf.int = TRUE,
                    risk.table = TRUE, risk.table.col = "strata",
                    risk.table.y.text = FALSE,
                    surv.median.line = "hv",
                    palette = c("#E5A727", "#4A90D9"),
                    legend.labs = c("High", "Low"),
                    legend.title = "Strata",
                    xlab = "Time (months)", ylab = "Survival Probability",
                    title = title_str,
                    ggtheme = theme_bw() + theme(
                      panel.grid = element_blank(),
                      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
                      legend.position = c(0.85, 0.9)))
    p$plot <- p$plot +
      annotate("text", x = 0, y = 0.15, label = hr_label,
               parse = TRUE, hjust = 0, size = 4.5, color = "#333333")
    print(p)
  }
  dev.off()
  cat(">>> Done!\n")
  
  result_df <- do.call(rbind, lapply(sig_results, function(r) {
    data.frame(Dataset = r$dataset, Tissue = r$tissue, Marker = r$marker,
               Cutoff = round(r$cutoff, 6), N_High = r$n_high, N_Low = r$n_low,
               P = signif(r$pval, 4), stringsAsFactors = FALSE)
  }))
  return(result_df)
}

sig_table <- run_all_prognosis(
  ratio_files = c("cellcount.xlsx", "celldensity.xlsx", "cellpercent.xlsx"),
  mfi_file = "MFI.xlsx",
  clin_file = "clin.xlsx",
  p_threshold = 0.05,
  output_pdf = "significant_KM_curves.pdf"
)
print(sig_table)
