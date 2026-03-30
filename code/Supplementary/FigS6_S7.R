#!/usr/bin/env Rscript
# Figure S7 — Differential methylation enrichment and clinical relevance

# --- S7a-d: GO enrichment of DMC/DMR genes ---

# 09_Targeted_Enrichment_by_Thesis_v2.R
# 
# 

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tidyverse)
library(ggplot2)

# ======================== 0.  ========================
# methylKit  "F_P_vs_F_N"  "F_N_vs_F_P"

my_dmc_dir <- "/home/root2/methylation/output/dma"
my_dmr_dir <- "/home/root2/methylation/output/Batch_DMR_Final_v5"
out_dir    <- "/home/root2/methylation/output/Thesis_Fig6.2_Enrichment"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# methylKit  treatment=, control= treatment_vs_control
# 
#

# DMC: _vs_ (FN_vs_FP, N_N_vs_N_P, P_N_vs_P_P)
DMC_PARENTAL  <- "FN_vs_FP"
DMC_MAINTAIN  <- "P_N_vs_P_P"  # CSC DMC (6.2E, F)
DMC_ACQUIRE   <- "N_N_vs_N_P"

DMR_ACQUIRE   <- "N_N_vs_N_P"
DMR_PARENTAL  <- "N_B_vs_P_B"

# ======================== 1.  ========================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr     <- genes(txdb)
promoters_gr <- promoters(genes_gr, upstream = 3000, downstream = 3000)

get_gene_annotation <- function(query_gr) {
  if (length(query_gr) == 0) return(list(all_entrez = c(), prom_entrez = c()))
  hits_all    <- distanceToNearest(query_gr, genes_gr)
  all_entrez  <- names(genes_gr)[subjectHits(hits_all)]
  hits_prom   <- findOverlaps(query_gr, promoters_gr)
  prom_entrez <- names(promoters_gr)[subjectHits(hits_prom)]
  list(all_entrez = unique(na.omit(all_entrez)),
       prom_entrez = unique(na.omit(prom_entrez)))
}

# ======================== 2.  +  ========================
read_and_annotate <- function(file_path, direction = "hyper") {
  if (!file.exists(file_path)) {
    return(list(all_entrez = c(), prom_entrez = c()))
  }
  df  <- read.csv(file_path)
  sub <- if (direction == "hyper") dplyr::filter(df, meth.diff > 0) else dplyr::filter(df, meth.diff < 0)
  if (nrow(sub) == 0) return(list(all_entrez = c(), prom_entrez = c()))
  gr <- makeGRangesFromDataFrame(sub, keep.extra.columns = TRUE)
  get_gene_annotation(gr)
}

# ======================== 3.  +  ( GO-only  GO+KEGG) ========================
run_enrichment <- function(gene_ids, title, save_prefix, include_kegg = TRUE) {
  if (length(gene_ids) < 5) {
    return(NULL)
  }
  
  # --- GO ---
  ego   <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  GO_df <- as.data.frame(ego)
  if (nrow(GO_df) > 0) {
    if (!"ONTOLOGY" %in% colnames(GO_df)) GO_df$ONTOLOGY <- ego@ontology
    go_sub <- GO_df %>% group_by(ONTOLOGY) %>% arrange(p.adjust) %>%
      slice_head(n = 5) %>% ungroup()
    write.csv(GO_df, file.path(out_dir, paste0(save_prefix, "_GO.csv")), row.names = FALSE)
  } else { go_sub <- data.frame() }
  
  # --- KEGG () ---
  kegg_sub <- data.frame()
  if (include_kegg) {
    ekegg   <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.05)
    KEGG_df <- as.data.frame(ekegg)
    if (nrow(KEGG_df) > 0) {
      KEGG_df$ONTOLOGY <- "KEGG"
      kegg_sub <- KEGG_df %>% arrange(p.adjust) %>% slice_head(n = 5)
      write.csv(KEGG_df, file.path(out_dir, paste0(save_prefix, "_KEGG.csv")), row.names = FALSE)
    }
  }
  
  plot_data <- bind_rows(go_sub, kegg_sub)
  if (nrow(plot_data) == 0) { message("    no significant enrichment"); return(NULL) }
  
  ont_levels <- if (include_kegg) c("BP", "CC", "MF", "KEGG") else c("BP", "CC", "MF")
  
  plot_data <- plot_data %>%
    separate(GeneRatio, into = c("num", "den"), sep = "/") %>%
    mutate(RichFactor   = as.numeric(num) / as.numeric(den),
           Description  = str_wrap(Description, width = 45),
           Description  = fct_reorder(Description, RichFactor),
           ONTOLOGY     = factor(ONTOLOGY, levels = ont_levels))
  
  p <- ggplot(plot_data, aes(x = RichFactor, y = Description)) +
    geom_point(aes(size = Count, fill = p.adjust),
               shape = 21, color = "grey30", stroke = 0.6, alpha = 0.85) +
    scale_fill_gradientn(
      colours = c("#b2182b", "#d6604d", "#f4a582", "#92c5de", "#2166ac"),
      name = "p.adjust", guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    labs(title = title, x = "Rich Factor", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey75"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text   = element_text(color = "black"),
      axis.text.y = element_text(size = 10, lineheight = 0.8),
      axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14,
                                  margin = margin(b = 15)),
      strip.background = element_rect(fill = "#f4f4f4", color = "grey40", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      panel.border = element_rect(color = "grey40", linewidth = 0.8),
      legend.position = "right", legend.key = element_blank()
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  
  h <- max(6, nrow(plot_data) * 0.45)
  ggsave(file.path(out_dir, paste0(save_prefix, ".pdf")), p,
         width = 10, height = h, limitsize = FALSE)
  return(p)
}

# ======================== 4. 6.2 D-H  ========================

# ====== 6.2D ======
dmc_file <- file.path(my_dmc_dir, DMC_PARENTAL, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hyper")
run_enrichment(res$all_entrez,
               "Fig6.2D: FN vs FP | Hyper DMC (Global)",
               "Fig6.2D_Parental_DMC_Hyper_Global", include_kegg = FALSE)
run_enrichment(res$prom_entrez,
               "Fig6.2D: FN vs FP | Hyper DMC (Promoter)",
               "Fig6.2D_Parental_DMC_Hyper_Promoter", include_kegg = FALSE)

# ====== 6.2E ======
dmc_file <- file.path(my_dmc_dir, DMC_MAINTAIN, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hyper")
run_enrichment(res$all_entrez,
               "Fig6.2E: P_N vs P_P | Hyper DMC (Global)",
               "Fig6.2E_Maintain_DMC_Hyper_Global", include_kegg = TRUE)
run_enrichment(res$prom_entrez,
               "Fig6.2E: P_N vs P_P | Hyper DMC (Promoter)",
               "Fig6.2E_Maintain_DMC_Hyper_Promoter", include_kegg = TRUE)

# ====== 6.2F ======
dmc_file <- file.path(my_dmc_dir, DMC_MAINTAIN, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hypo")
run_enrichment(res$all_entrez,
               "Fig6.2F: P_N vs P_P | Hypo DMC (Global)",
               "Fig6.2F_Maintain_DMC_Hypo_Global", include_kegg = TRUE)
run_enrichment(res$prom_entrez,
               "Fig6.2F: P_N vs P_P | Hypo DMC (Promoter)",
               "Fig6.2F_Maintain_DMC_Hypo_Promoter", include_kegg = TRUE)

# ====== 6.2G ======
dmr_file <- file.path(my_dmr_dir, DMR_ACQUIRE, "DMR_Significant.csv")
res <- read_and_annotate(dmr_file, "hypo")
run_enrichment(res$all_entrez,
               "Fig6.2G: N_N vs N_P | Hypo DMR (Global)",
               "Fig6.2G_Acquire_DMR_Hypo_Global", include_kegg = TRUE)
run_enrichment(res$prom_entrez,
               "Fig6.2G: N_N vs N_P | Hypo DMR (Promoter)",
               "Fig6.2G_Acquire_DMR_Hypo_Promoter", include_kegg = TRUE)

# ====== 6.2H ======
dmr_file <- file.path(my_dmr_dir, DMR_PARENTAL, "DMR_Significant.csv")
res <- read_and_annotate(dmr_file, "hypo")
run_enrichment(res$all_entrez,
               "Fig6.2H: N_B vs P_B | Hypo DMR (Global)",
               "Fig6.2H_Parental_DMR_Hypo_Global", include_kegg = FALSE)
run_enrichment(res$prom_entrez,
               "Fig6.2H: N_B vs P_B | Hypo DMR (Promoter)",
               "Fig6.2H_Parental_DMR_Hypo_Promoter", include_kegg = FALSE)

# ======================== 5.  ========================

dmc_file <- file.path(my_dmc_dir, DMC_ACQUIRE, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hypo")
run_enrichment(res$prom_entrez,
               "N_N vs N_P | Hypo DMC (Promoter)",
               "Supp_NNvsNP_DMC_Hypo_Promoter", include_kegg = TRUE)
run_enrichment(res$all_entrez,
               "N_N vs N_P | Hypo DMC (Global)",
               "Supp_NNvsNP_DMC_Hypo_Global", include_kegg = TRUE)

# Hypo DMC
dmc_file <- file.path(my_dmc_dir, DMC_PARENTAL, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hypo")
run_enrichment(res$all_entrez,
               "FN vs FP | Hypo DMC (Global)",
               "Supp_FNvsFP_DMC_Hypo_Global", include_kegg = TRUE)

# DMR Hyper
# DMR: N_N_vs_N_P, P_N_vs_P_P, N_B_vs_P_B
for (comp in c(DMR_ACQUIRE, "P_N_vs_P_P", DMR_PARENTAL)) {
  dmr_file <- file.path(my_dmr_dir, comp, "DMR_Significant.csv")
  res <- read_and_annotate(dmr_file, "hyper")
  run_enrichment(res$all_entrez,
                 paste0(comp, " | Hyper DMR (Global)"),
                 paste0("Supp_", comp, "_DMR_Hyper_Global"), include_kegg = TRUE)
}

# N_N vs N_P Hyper DMC
dmc_file <- file.path(my_dmc_dir, DMC_ACQUIRE, "Significant_DMCs_all.csv")
res <- read_and_annotate(dmc_file, "hyper")
run_enrichment(res$all_entrez,
               "N_N vs N_P | Hyper DMC (Global)",
               "Supp_NNvsNP_DMC_Hyper_Global", include_kegg = TRUE)

# --- Supplementary: Entropy full pipeline (scRNA entropy + ChromHMM + TCGA mediation) ---

#!/usr/bin/env Rscript
# (Entropy Analysis Full Pipeline)
#
#   Rscript entropy_full_pipeline.R [part]
#
# part = 1            (download, Python/scanpy)
# part = 2  EPIC×ChromHMM   (root2)
# part = 3  TCGA +     (download)
# part = 4  ×             (download, Python/scanpy)
# part = all  ()
#

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
RUN_PART <- if (length(args) > 0) args[1] else "all"

FIG_DIR <- "/home/download/csc_article/fig6/entropy/"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Entropy Analysis Pipeline\n")
cat(sprintf("  Running: Part %s\n", RUN_PART))
cat(sprintf("  Output:  %s\n", FIG_DIR))
cat("============================================================\n\n")

# ║  PART 1:  (download)                       ║
# ║  : Python + scanpy                                        ║
if (RUN_PART %in% c("1", "all")) {

cat("\n", strrep("=", 60), "\n")
cat("  PART 1: Single-cell Transcriptomic Entropy\n")
cat(strrep("=", 60), "\n")

py_script <- file.path(FIG_DIR, "run_part1_sc_entropy.py")
writeLines(con = py_script, text = '#!/usr/bin/env python3
"""Part 1: Single-cell transcriptomic entropy"""
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os, warnings; warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family":"Arial","font.size":10,
    "pdf.fonttype":42,"savefig.dpi":300,"savefig.bbox":"tight"})

fig_dir = "/home/download/csc_article/fig6/entropy/"
os.makedirs(fig_dir, exist_ok=True)

print("Loading data...")
adata = sc.read_h5ad("/data2/adata/adata_subset_0214.h5ad")
adata = adata[adata.obs["malignant"] == "tumor"].copy()

if "raw_counts" in adata.layers:
    adata.X = adata.layers["raw_counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)

#  CSC -> CSC(classic)
adata.obs["final_type"] = adata.obs["final_type"].replace({"CSC": "CSC(classic)"})
print(f"Tumor cells: {adata.shape[0]}")
print(adata.obs["final_type"].value_counts().to_string())

# --- Shannon ---
print("\\nComputing transcriptomic entropy...")
X = adata.X
if scipy.sparse.issparse(X): X = X.toarray()

entropies = np.zeros(X.shape[0])
for i in range(X.shape[0]):
    x = X[i, :]
    p = x[x > 0]; p = p / p.sum()
    entropies[i] = -np.sum(p * np.log2(p))
    if (i+1) % 10000 == 0: print(f"  {i+1}/{X.shape[0]}")

adata.obs["transcriptomic_entropy"] = entropies
print(f"  Range: {entropies.min():.2f} - {entropies.max():.2f}")

type_order = [
    "CSC(classic)", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)",
    "Hepatocytes (Stem)", "Hepatocytes (Metabolic-H)", "Hepatocytes (Metabolic-L)",
    "Hepatocytes (cycling)", "Hepatocytes", "Hepatocytes (Quiescent)",
    "Cholangiocytes",
]
type_order = [t for t in type_order if t in adata.obs["final_type"].values]

stats = adata.obs.groupby("final_type")["transcriptomic_entropy"].agg(["mean","median","std","count"])
stats = stats.loc[[t for t in type_order if t in stats.index]]
print("\\n" + stats.to_string())

adata.obs["CSC_group"] = adata.obs["final_type"].apply(
    lambda x: "CSC" if "CSC" in str(x) else "non-CSC tumor")
csc_ent = adata.obs.loc[adata.obs["CSC_group"]=="CSC", "transcriptomic_entropy"]
non_ent = adata.obs.loc[adata.obs["CSC_group"]=="non-CSC tumor", "transcriptomic_entropy"]
U, p = mannwhitneyu(csc_ent, non_ent, alternative="greater")
print(f"\\nCSC vs non-CSC: CSC_mean={csc_ent.mean():.4f}, nonCSC_mean={non_ent.mean():.4f}, p={p:.2e}")

# --- : Boxplot ---
color_map = {
    "CSC(classic)":"#D62728", "CSC(progenitor_cycling)":"#FF7F0E",
    "CSC(progenitor_quiet)":"#E377C2", "Hepatocytes (Stem)":"#8C564B",
    "Hepatocytes (Metabolic-H)":"#1F77B4", "Hepatocytes (Metabolic-L)":"#AEC7E8",
    "Hepatocytes (cycling)":"#17BECF", "Hepatocytes":"#9467BD",
    "Hepatocytes (Quiescent)":"#C5B0D5", "Cholangiocytes":"#2CA02C",
}
plot_df = adata.obs[adata.obs["final_type"].isin(type_order)].copy()
plot_df["final_type"] = pd.Categorical(plot_df["final_type"], categories=type_order, ordered=True)

fig, ax = plt.subplots(figsize=(12, 5))
bp_data = [plot_df.loc[plot_df["final_type"]==t, "transcriptomic_entropy"].values for t in type_order]
parts = ax.boxplot(bp_data, positions=range(len(type_order)), widths=0.6,
                   patch_artist=True, showfliers=False,
                   medianprops=dict(color="black", linewidth=1.5))
for patch, t in zip(parts["boxes"], type_order):
    patch.set_facecolor(color_map.get(t, "#999")); patch.set_alpha(0.7)

# strip (500)
for i, t in enumerate(type_order):
    vals = plot_df.loc[plot_df["final_type"]==t, "transcriptomic_entropy"].values
    jitter = np.random.normal(0, 0.08, size=len(vals))
    if len(vals) > 500:
        idx = np.random.choice(len(vals), 500, replace=False)
        vals, jitter = vals[idx], jitter[idx]
    ax.scatter(i + jitter, vals, s=3, alpha=0.3, c=color_map.get(t, "#999"), zorder=2)

ax.set_xticks(range(len(type_order)))
short = [t.replace("Hepatocytes","Hep") for t in type_order]
ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Transcriptomic Entropy (bits)", fontsize=11)
ax.set_title("Transcriptomic Entropy by Cell Type (Tumor Cells)", fontsize=12, fontweight="bold")
n_csc = sum(1 for t in type_order if "CSC" in t)
ax.axvline(n_csc - 0.5, color="gray", ls=":", lw=0.8, alpha=0.5)
ax.text(n_csc/2 - 0.5, ax.get_ylim()[1]*0.98, "CSC", ha="center", fontsize=10,
        color="#D62728", fontweight="bold")
ax.text((n_csc + len(type_order))/2 - 0.5, ax.get_ylim()[1]*0.98, "non-CSC tumor",
        ha="center", fontsize=10, color="#1F77B4", fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
for fmt in ["pdf","png"]:
    plt.savefig(f"{fig_dir}Boxplot_TranscriptomicEntropy_CellType.{fmt}",
                dpi=300 if fmt=="pdf" else 200, bbox_inches="tight")
plt.close()
print("  -> Boxplot saved")

# --- Violin: CSC vs non-CSC ---
fig, ax = plt.subplots(figsize=(4, 5))
data = [csc_ent.values, non_ent.values]
vp = ax.violinplot(data, positions=[0,1], showmeans=True, showmedians=True)
vp["bodies"][0].set_facecolor("#D62728"); vp["bodies"][1].set_facecolor("#1F77B4")
for b in vp["bodies"]: b.set_alpha(0.6)
ax.set_xticks([0,1]); ax.set_xticklabels(["CSC","non-CSC\\ntumor"], fontsize=11)
ax.set_ylabel("Transcriptomic Entropy (bits)")
ax.set_title(f"p = {p:.2e} (Mann-Whitney)", fontsize=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
for fmt in ["pdf","png"]:
    plt.savefig(f"{fig_dir}Violin_Entropy_CSCvsNonCSC.{fmt}", dpi=300 if fmt=="pdf" else 200)
plt.close()
print("  -> Violin saved")

# --- Pairwise comparisons ---
print("\\n--- Pairwise comparisons ---")
csc_types = [t for t in type_order if "CSC" in t]
non_types = [t for t in type_order if "CSC" not in t]
pairs = []
for ct in csc_types:
    for nt in non_types:
        v1 = plot_df.loc[plot_df["final_type"]==ct, "transcriptomic_entropy"].values
        v2 = plot_df.loc[plot_df["final_type"]==nt, "transcriptomic_entropy"].values
        if len(v1)>5 and len(v2)>5:
            _, pp = mannwhitneyu(v1, v2, alternative="greater")
            pairs.append({"CSC_type":ct, "nonCSC_type":nt,
                          "CSC_mean":v1.mean(), "nonCSC_mean":v2.mean(), "pval":pp})
df_p = pd.DataFrame(pairs)
if len(df_p)>0:
    _, df_p["FDR"], _, _ = multipletests(df_p["pval"], method="fdr_bh")
    df_p["sig"] = df_p["FDR"].apply(lambda q: "***" if q<0.001 else ("**" if q<0.01 else ("*" if q<0.05 else "ns")))
    print(df_p.to_string(index=False))
    df_p.to_csv(f"{fig_dir}Pairwise_Entropy_Comparisons.csv", index=False)

adata.obs[["final_type","CSC_group","transcriptomic_entropy"]].to_csv(
    f"{fig_dir}Cell_TranscriptomicEntropy.csv")
print(f"\\nPart 1 done. Outputs in {fig_dir}")
')

cat("  Running Python script...\n")
ret <- system(sprintf("python3 %s", py_script))
if (ret != 0) cat("  WARNING: Python script returned non-zero exit code!\n")

} # end Part 1

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

# comparisonboxplot
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

# ║  PART 3: TCGA +  (download)                ║
if (RUN_PART %in% c("3", "all")) {

cat("\n", strrep("=", 60), "\n")
cat("  PART 3: TCGA Paired Entropy Analysis\n")
cat(strrep("=", 60), "\n")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TCGAbiolinks)
})

meth_rds <- "/home/download/JWD_LIHC/tcga_lihc_methylation_final.rds"
message(">>> Loading TCGA methylation: ", meth_rds)
meth_obj <- readRDS(meth_rds)
if (inherits(meth_obj, "SummarizedExperiment") || inherits(meth_obj, "RangedSummarizedExperiment")) {
  tcga_beta <- assay(meth_obj)
} else {
  tcga_beta <- meth_obj
}
rm(meth_obj); gc(verbose = FALSE)
message(sprintf("  450K: %d probes × %d samples", nrow(tcga_beta), ncol(tcga_beta)))

# GDCprepare
tcga_rdata <- "/home/download/csc_article/TCGA_LIHC_entropy.RData"

if (file.exists(tcga_rdata)) {
  message(">>> Loading cached expression: ", tcga_rdata)
  load(tcga_rdata)
} else {
  message(">>> Fetching TCGA-LIHC expression via GDCprepare...")
  tryCatch({
    query <- GDCquery(project = "TCGA-LIHC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
    data_obj <- GDCprepare(query)
    tcga_expr <- assay(data_obj, "fpkm_unstrand")
    rownames(tcga_expr) <- as.data.frame(rowData(data_obj))$gene_name
    save(tcga_expr, file = tcga_rdata)
    message(sprintf("  Expression: %d genes × %d samples", nrow(tcga_expr), ncol(tcga_expr)))
    message("  Cached to: ", tcga_rdata)
  }, error = function(e) {
    message("  ERROR loading expression: ", e$message)
    tryCatch({
      query <- GDCquery(project = "TCGA-LIHC",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")
      data_obj <- GDCprepare(query, directory = "/home/download/JWD_LIHC/GDCdata")
      tcga_expr <<- assay(data_obj, "fpkm_unstrand")
      rownames(tcga_expr) <- as.data.frame(rowData(data_obj))$gene_name
      save(tcga_expr, file = tcga_rdata)
    }, error = function(e2) {
      message("  FAILED: ", e2$message)
    })
  })
}

if (exists("tcga_expr")) {
  # --- 3c.  ---
  colnames(tcga_beta) <- substr(colnames(tcga_beta), 1, 15)
  colnames(tcga_expr) <- substr(colnames(tcga_expr), 1, 15)

  # tumor samples only (barcode position 14 == '0')
  tumor_beta <- colnames(tcga_beta)[substr(colnames(tcga_beta), 14, 14) == "0"]
  tumor_expr <- colnames(tcga_expr)[substr(colnames(tcga_expr), 14, 14) == "0"]
  common <- sort(intersect(tumor_beta, tumor_expr))
  common <- common[!duplicated(common)]
  tcga_expr <- tcga_expr[, !duplicated(colnames(tcga_expr))]
  tcga_beta <- tcga_beta[, !duplicated(colnames(tcga_beta))]
  common <- sort(intersect(colnames(tcga_beta), colnames(tcga_expr)))
  common <- common[substr(common, 14, 14) == "0"]
  message(sprintf("  Paired tumor samples: %d", length(common)))

  # --- 3d.  ---
  message(">>> Computing transcriptomic entropy...")
  expr_sub <- tcga_expr[, common]
  expr_sub <- expr_sub[!duplicated(rownames(expr_sub)), ]

  expr_entropy <- sapply(common, function(s) {
    x <- as.numeric(expr_sub[, s])
    x <- x[x > 0]
    if (sum(x) == 0) return(NA)
    p <- x / sum(x)
    -sum(p * log2(p))
  })

  # --- 3e.  ---
  message(">>> Computing methylation entropy...")
  beta_sub <- tcga_beta[, common]

  meth_entropy <- sapply(common, function(s) {
    b <- as.numeric(beta_sub[, s])
    b <- b[!is.na(b) & b > 0.01 & b < 0.99]
    if (length(b) == 0) return(NA)
    h <- -(b * log2(b) + (1 - b) * log2(1 - b))
    mean(h)
  })

  # --- 3f. DNMT ---
  dnmt_genes <- c("DNMT1", "DNMT3A", "DNMT3B")
  dnmt_found <- dnmt_genes[dnmt_genes %in% rownames(expr_sub)]
  message(sprintf("  DNMT genes found: %s", paste(dnmt_found, collapse = ", ")))

  dnmt_mat <- t(expr_sub[dnmt_found, common, drop = FALSE])
  dnmt_z <- scale(dnmt_mat)
  dnmt_score <- rowMeans(dnmt_z, na.rm = TRUE)

  q33 <- quantile(dnmt_score, 0.33, na.rm = TRUE)
  q66 <- quantile(dnmt_score, 0.66, na.rm = TRUE)
  dnmt_group <- ifelse(dnmt_score < q33, "DNMT_low",
                ifelse(dnmt_score > q66, "DNMT_high", "DNMT_mid"))

  # --- 3g.  ---
  df_tcga <- data.frame(
    sample = common,
    expr_entropy = expr_entropy,
    meth_entropy = meth_entropy,
    DNMT_score = dnmt_score,
    DNMT_group = dnmt_group,
    stringsAsFactors = FALSE
  )
  for (g in dnmt_found) {
    df_tcga[[g]] <- as.numeric(expr_sub[g, common])
  }
  df_tcga <- df_tcga[complete.cases(df_tcga[, c("expr_entropy","meth_entropy","DNMT_score")]), ]
  message(sprintf("  Complete cases: %d", nrow(df_tcga)))

  # --- 3h.  ---
  rho_all <- cor.test(df_tcga$meth_entropy, df_tcga$expr_entropy, method = "spearman")
  rho_dnmt_meth <- cor.test(df_tcga$DNMT_score, df_tcga$meth_entropy, method = "spearman")
  rho_dnmt_expr <- cor.test(df_tcga$DNMT_score, df_tcga$expr_entropy, method = "spearman")

  cat("\n--- Correlation Results ---\n")
  cat(sprintf("  meth_entropy ~ expr_entropy: rho=%.4f, p=%.2e\n",
              rho_all$estimate, rho_all$p.value))
  cat(sprintf("  DNMT_score ~ meth_entropy:   rho=%.4f, p=%.2e\n",
              rho_dnmt_meth$estimate, rho_dnmt_meth$p.value))
  cat(sprintf("  DNMT_score ~ expr_entropy:   rho=%.4f, p=%.2e\n",
              rho_dnmt_expr$estimate, rho_dnmt_expr$p.value))

  cat("\n--- Stratified by DNMT level ---\n")
  for (grp in c("DNMT_low","DNMT_mid","DNMT_high")) {
    sub <- df_tcga[df_tcga$DNMT_group == grp, ]
    if (nrow(sub) > 10) {
      rt <- cor.test(sub$meth_entropy, sub$expr_entropy, method = "spearman")
      cat(sprintf("  [%10s] n=%3d, rho=%.4f, p=%.3e\n", grp, nrow(sub), rt$estimate, rt$p.value))
    }
  }

  # --- 3i.  ---
  message("\n>>> Plotting...")
  colors_grp <- c("DNMT_low"="#2196F3","DNMT_mid"="#9E9E9E","DNMT_high"="#F44336")
  df_tcga$DNMT_group <- factor(df_tcga$DNMT_group, levels = c("DNMT_low","DNMT_mid","DNMT_high"))

  # Scatter: meth vs expr entropy
  p1 <- ggplot(df_tcga, aes(meth_entropy, expr_entropy, color = DNMT_group)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = colors_grp) +
    labs(title = sprintf("TCGA-LIHC (n=%d)\nSpearman ρ = %.3f, p = %.2e",
                         nrow(df_tcga), rho_all$estimate, rho_all$p.value),
         x = "Methylation Entropy (mean CpG H)", y = "Transcriptomic Entropy (bits)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(file.path(FIG_DIR, "Scatter_MethEnt_vs_ExprEnt.pdf"), p1, width = 7, height = 6, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, "Scatter_MethEnt_vs_ExprEnt.png"), p1, width = 7, height = 6, dpi = 300)
  message("  -> Scatter saved")

  # TriPanel
  p_tri1 <- ggplot(df_tcga, aes(DNMT_score, meth_entropy)) +
    geom_point(color = "#607D8B", alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(title = sprintf("DNMT → Meth Entropy\nρ=%.3f", rho_dnmt_meth$estimate),
         x = "DNMT Score", y = "Meth Entropy") +
    theme_bw(base_size = 10)

  p_tri2 <- ggplot(df_tcga, aes(meth_entropy, expr_entropy)) +
    geom_point(color = "#607D8B", alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(title = sprintf("Meth Entropy → Expr Entropy\nρ=%.3f", rho_all$estimate),
         x = "Meth Entropy", y = "Expr Entropy") +
    theme_bw(base_size = 10)

  p_tri3 <- ggplot(df_tcga, aes(DNMT_score, expr_entropy)) +
    geom_point(color = "#607D8B", alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    labs(title = sprintf("DNMT → Expr Entropy\nρ=%.3f", rho_dnmt_expr$estimate),
         x = "DNMT Score", y = "Expr Entropy") +
    theme_bw(base_size = 10)

  p_tri <- p_tri1 | p_tri2 | p_tri3
  ggsave(file.path(FIG_DIR, "TriPanel_DNMT_MethEnt_ExprEnt.pdf"), p_tri,
         width = 15, height = 5, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, "TriPanel_DNMT_MethEnt_ExprEnt.png"), p_tri, width = 15, height = 5, dpi = 300)
  message("  -> TriPanel saved")

  write.csv(df_tcga, file.path(FIG_DIR, "TCGA_Entropy_ForMediation.csv"), row.names = FALSE)

  # --- 3j.  ---
  message("\n>>> Running mediation analysis...")
  tryCatch({
    library(mediation)

    df_tcga$DNMT_z <- scale(df_tcga$DNMT_score)
    df_tcga$meth_z <- scale(df_tcga$meth_entropy)
    df_tcga$expr_z <- scale(df_tcga$expr_entropy)

    model_m <- lm(meth_z ~ DNMT_z, data = df_tcga)
    model_y <- lm(expr_z ~ DNMT_z + meth_z, data = df_tcga)

    set.seed(42)
    med_result <- mediate(model_m, model_y,
                          treat = "DNMT_z", mediator = "meth_z",
                          boot = TRUE, sims = 5000)

    cat("\n--- Mediation Results ---\n")
    print(summary(med_result))

    sink(file.path(FIG_DIR, "Mediation_Results.txt"))
    print(summary(med_result))
    sink()

    pdf(file.path(FIG_DIR, "Mediation_Plot.pdf"), width = 6, height = 4)
    plot(med_result, main = "Mediation: DNMT → Meth Entropy → Expr Entropy")
    dev.off()
    message("  -> Mediation results saved")

  }, error = function(e) {
    message("  Mediation analysis failed: ", e$message)
    message("  Install with: install.packages('mediation')")
  })

} else {
  message("  SKIPPED: TCGA expression data not available")
}

message("\n>>> Part 3 complete.")

} # end Part 3

# ║  PART 4:  ×  (download)                       ║
# ║  : Python + scanpy + Part 1                         ║
if (RUN_PART %in% c("4", "all")) {

cat("\n", strrep("=", 60), "\n")
cat("  PART 4: Metabolic Scoring × Entropy\n")
cat(strrep("=", 60), "\n")

py_script4 <- file.path(FIG_DIR, "run_part4_metabolic.py")
writeLines(con = py_script4, text = '#!/usr/bin/env python3
"""Part 4: Metabolic scoring × transcriptomic entropy"""
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import os, warnings; warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family":"Arial","font.size":10,
    "pdf.fonttype":42,"savefig.dpi":300,"savefig.bbox":"tight"})

fig_dir = "/home/download/csc_article/fig6/entropy/"

print("Loading data...")
adata = sc.read_h5ad("/data2/adata/adata_subset_0214.h5ad")
adata = adata[adata.obs["malignant"] == "tumor"].copy()
if "raw_counts" in adata.layers:
    adata.X = adata.layers["raw_counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# CSC -> CSC(classic)
adata.obs["final_type"] = adata.obs["final_type"].replace({"CSC": "CSC(classic)"})

ent_file = f"{fig_dir}Cell_TranscriptomicEntropy.csv"
if os.path.exists(ent_file):
    ent_df = pd.read_csv(ent_file, index_col=0)
    common = adata.obs.index.intersection(ent_df.index)
    adata = adata[common].copy()
    adata.obs["transcriptomic_entropy"] = ent_df.loc[common, "transcriptomic_entropy"].values
    print(f"  Loaded entropy for {len(common)} cells")
else:
    print(f"  ERROR: {ent_file} not found! Run Part 1 first.")
    exit(1)

adata.obs["CSC_group"] = adata.obs["final_type"].apply(
    lambda x: "CSC" if "CSC" in str(x) else "non-CSC tumor")

gene_sets = {
    "SAM_cycle": ["MAT1A","MAT2A","MAT2B","AHCY","AHCYL1","MTR","MTRR",
                  "BHMT","BHMT2","CBS","CTH","GNMT","NNMT","PEMT","GAMT"],
    "one_carbon": ["MTHFR","MTHFD1","MTHFD2","MTHFD1L","SHMT1","SHMT2",
                   "TYMS","DHFR","GART","ATIC","ALDH1L1","ALDH1L2"],
    "methionine_metabolism": ["MAT1A","MAT2A","AHCY","MTR","BHMT",
                              "CBS","CTH","DNMT1","DNMT3A","DNMT3B"],
}

print("\\nScoring gene sets...")
for name, genes in gene_sets.items():
    found = [g for g in genes if g in adata.var_names]
    print(f"  {name}: {len(found)}/{len(genes)} found")
    if len(found) >= 3:
        sc.tl.score_genes(adata, found, score_name=f"score_{name}")

print("\\n--- Correlations with entropy ---")
score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
results = []
for sc_col in score_cols:
    vals = adata.obs[sc_col].values
    ent = adata.obs["transcriptomic_entropy"].values
    mask = ~(np.isnan(vals) | np.isnan(ent))
    rho, p = spearmanr(vals[mask], ent[mask])
    name = sc_col.replace("score_","")
    results.append({"gene_set":name, "rho":rho, "pval":p})
    print(f"  {name:25s}  rho={rho:.4f}  p={p:.2e}")

    fig, ax = plt.subplots(figsize=(6, 5))
    for grp, color in [("CSC","#D62728"), ("non-CSC tumor","#1F77B4")]:
        sub = adata.obs[adata.obs["CSC_group"]==grp]
        if len(sub) > 2000: sub = sub.sample(2000, random_state=42)
        ax.scatter(sub[sc_col], sub["transcriptomic_entropy"],
                   c=color, s=5, alpha=0.2, label=grp)
    ax.set_xlabel(f"{name} score"); ax.set_ylabel("Transcriptomic Entropy")
    ax.set_title(f"{name} vs Entropy\\nSpearman ρ={rho:.3f}, p={p:.2e}")
    ax.legend(fontsize=9, markerscale=3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}Scatter_{name}_vs_Entropy.pdf", bbox_inches="tight")
    plt.savefig(f"{fig_dir}Scatter_{name}_vs_Entropy.png", dpi=200, bbox_inches="tight")
    plt.close()

pd.DataFrame(results).to_csv(f"{fig_dir}Metabolic_Entropy_Correlations.csv", index=False)
print(f"\\nPart 4 done. Outputs in {fig_dir}")
')

cat("  Running Python script...\n")
ret <- system(sprintf("python3 %s", py_script4))
if (ret != 0) cat("  WARNING: Python script returned non-zero exit code!\n")

} # end Part 4

cat("\n", strrep("=", 60), "\n")
cat("  ALL DONE!\n")
cat(strrep("=", 60), "\n")
cat(sprintf("
Output directory:
  Figures: %s
  
Part descriptions:
  Part 1: scRNA-seq entropy
    - Boxplot_TranscriptomicEntropy_CellType.pdf
    - Violin_Entropy_CSCvsNonCSC.pdf
    - Cell_TranscriptomicEntropy.csv
    
  Part 2: EPIC methylation entropy x ChromHMM  
    - Fig_MethylationEntropy_ChromHMM.pdf
    - EPIC_Entropy_by_ChromHMM.csv
    
  Part 3: TCGA paired entropy
    - Scatter_MethEnt_vs_ExprEnt.pdf
    - TriPanel_DNMT_MethEnt_ExprEnt.pdf
    - Mediation_Results.txt
    
  Part 4: Metabolic scoring
    - Scatter_*_vs_Entropy.pdf
    - Metabolic_Entropy_Correlations.csv

Usage:
  Rscript entropy_full_pipeline.R all
  
  Rscript entropy_full_pipeline.R 1
  Rscript entropy_full_pipeline.R 2   # EPIC×ChromHMM (root2)
  Rscript entropy_full_pipeline.R 3
  Rscript entropy_full_pipeline.R 4
", FIG_DIR))

# --- Supplementary: Compartment analysis (compartmap) ---

#!/usr/bin/env Rscript
# compartmap:  scRNA-seq  A/B Compartment
# 1.  /home/root2/methylation/output/Compartmap/

# ======================== 0.  ========================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("compartmap", "SummarizedExperiment", "GenomicRanges",
              "rtracklayer", "EnsDb.Hsapiens.v86")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
for (pkg in c("ggplot2", "patchwork", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(compartmap)
library(SummarizedExperiment)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

# ======================== 1.  ========================
output_dir <- "/home/root2/methylation/output/Compartmapsc"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cell_types <- c(
  "Bipotent Progenitor",
  "Hepatocytes",
  "CSC",
  "CSC(progenitor_cycling)",
  "CSC(progenitor_quiet)",
  "Hepatocytes (Metabolic-L)"
)

chromosomes <- paste0("chr", 1:22)
resolution  <- 1e6
max_cells   <- 2000

save_plot <- function(p, path_no_ext, w, h) {
  ggsave(paste0(path_no_ext, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(path_no_ext, ".png"), p, width = w, height = h, dpi = 300)
}

save_base_plot <- function(path_no_ext, w_in, h_in, expr) {
  # PDF
  cairo_pdf(paste0(path_no_ext, ".pdf"), width = w_in, height = h_in)
  tryCatch(expr, error = function(e) { plot.new(); text(0.5, 0.5, e$message) })
  dev.off()
  # PNG
  png(paste0(path_no_ext, ".png"), width = w_in * 150, height = h_in * 150, res = 150)
  tryCatch(expr, error = function(e) { plot.new(); text(0.5, 0.5, e$message) })
  dev.off()
}

# ======================== 2.  Seurat  ========================
# seu_obj <- readRDS("/path/to/seu_obj.rds")
print(table(seu_obj$final_type))

# ======================== 3.  rowRanges ========================
edb <- EnsDb.Hsapiens.v86
all_genes <- genes(edb, columns = c("gene_name", "gene_id"))
names(all_genes) <- all_genes$gene_name
seqlevelsStyle(all_genes) <- "UCSC"

gene_names    <- rownames(seu_obj)
matched_genes <- intersect(gene_names, names(all_genes))

gene_gr <- all_genes[matched_genes]
gene_gr <- gene_gr[!duplicated(names(gene_gr))]
matched_genes <- names(gene_gr)

# ======================== 4.  compartmap ========================

compartment_results <- list()

for (ct in cell_types) {
  message("\n", strrep("=", 60))
  message(strrep("=", 60))

  cells <- colnames(seu_obj)[seu_obj$final_type == ct]
  if (length(cells) == 0) { message("    [skip] no cells"); next }

  if (length(cells) > max_cells) {
    set.seed(42)
    cells <- sample(cells, max_cells)
  }

  counts_mat <- GetAssayData(seu_obj, layer = "counts")[matched_genes, cells]
  if (inherits(counts_mat, "sparseMatrix")) counts_mat <- as.matrix(counts_mat)

  rse <- SummarizedExperiment(
    assays    = list(counts = counts_mat),
    rowRanges = gene_gr[matched_genes],
    colData   = DataFrame(row.names = cells)
  )
  colnames(rse) <- cells
  message("    RSE: ", nrow(rse), " genes x ", ncol(rse), " cells")

  # TF-IDF
  tfidf_mat <- transformTFIDF(assay(rse, "counts"))
  assay(rse, "counts") <- t(tfidf_mat)

  ct_results <- list()
  for (chrom in chromosomes) {
    message("    ", chrom, "...", appendLF = FALSE)
    res <- tryCatch({
      comp <- scCompartments(
        rse, chr = chrom, res = resolution,
        group = TRUE, bootstrap = TRUE, num.bootstraps = 10,
        genome = "hg38", assay = "rna"
      )
      comp_fix <- fixCompartments(comp, min.conf = 0.8)
      message(" OK")
      comp_fix
    }, error = function(e) {
      message(" SKIP (", e$message, ")")
      NULL
    })

    if (!is.null(res)) {
      if (is(res, "GRangesList")) res <- res[[1]]
      ct_results[[chrom]] <- res
    }
  }

  if (length(ct_results) > 0) {
    all_gr <- do.call(c, unname(ct_results))
    compartment_results[[ct]] <- all_gr
    saveRDS(all_gr, file.path(output_dir, paste0("Comp_",
                                                 gsub("[^A-Za-z0-9]", "_", ct), ".rds")))
  }
}

# ======================== 5.  ========================

gr_to_dt <- function(gr, ct_name) {
  score_col <- if ("flip.score" %in% colnames(mcols(gr))) "flip.score" else "score"
  conf_col  <- if ("flip.conf.est" %in% colnames(mcols(gr))) "flip.conf.est" else
    if ("conf.est" %in% colnames(mcols(gr))) "conf.est" else NA

  dt <- data.table(
    chr         = as.character(seqnames(gr)),
    start       = start(gr),
    end         = end(gr),
    mid         = (start(gr) + end(gr)) / 2,
    eigenvector = mcols(gr)[[score_col]],
    compartment = gr$compartments,
    cell_type   = ct_name
  )
  if (!is.na(conf_col)) dt$confidence <- mcols(gr)[[conf_col]]
  dt
}

comp_dt_list <- lapply(names(compartment_results), function(ct) {
  gr_to_dt(compartment_results[[ct]], ct)
})
all_comp <- rbindlist(comp_dt_list, fill = TRUE)
fwrite(all_comp, file.path(output_dir, "All_Compartments.csv"))

# ======================== 6. plotAB  ========================

for (ct in names(compartment_results)) {
  gr <- compartment_results[[ct]]
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)

  for (chrom in c("chr1", "chr14", "chr19")) {
    chr_gr <- gr[seqnames(gr) == chrom]
    if (length(chr_gr) == 0) next

    save_base_plot(
      file.path(output_dir, paste0("plotAB_", ct_safe, "_", chrom)),
      8, 2.7,
      {
        plotAB(chr_gr, chr = chrom, with.ci = TRUE, median.conf = TRUE)
        title(main = paste0(ct, " - ", chrom))
      }
    )
  }
}

# ======================== 7.  ========================

heatmap_chroms <- c("chr1", "chr14", "chr19")

for (ct in cell_types) {
  if (!(ct %in% names(compartment_results))) next
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)

  cells <- colnames(seu_obj)[seu_obj$final_type == ct]
  if (length(cells) > max_cells) { set.seed(42); cells <- sample(cells, max_cells) }

  counts_mat <- GetAssayData(seu_obj, layer = "counts")[matched_genes, cells]
  if (inherits(counts_mat, "sparseMatrix")) counts_mat <- as.matrix(counts_mat)

  rse_tmp <- SummarizedExperiment(
    assays    = list(counts = counts_mat),
    rowRanges = gene_gr[matched_genes],
    colData   = DataFrame(row.names = cells)
  )
  colnames(rse_tmp) <- cells
  tfidf_tmp <- transformTFIDF(assay(rse_tmp, "counts"))
  assay(rse_tmp, "counts") <- t(tfidf_tmp)

  for (chrom in heatmap_chroms) {
    save_base_plot(
      file.path(output_dir, paste0("CorMatrix_", ct_safe, "_", chrom)),
      5.5, 5.5,
      {
        rmt <- getDenoisedCorMatrix(
          rse_tmp, res = resolution, chr = chrom,
          genome = "hg38", assay = "rna", iter = 2
        )
        plotCorMatrix(rmt, uppertri = TRUE)
        title(main = paste0(ct, " - ", chrom), cex.main = 1.1)
      }
    )
  }
}

# ======================== 8.  ========================

# ---- 8a.  ----
plot_track <- function(dt, ct_name, chrom) {
  sub <- as.data.frame(dt[cell_type == ct_name & chr == chrom])
  if (nrow(sub) == 0) return(NULL)
  ggplot(sub, aes(x = mid / 1e6, y = eigenvector)) +
    geom_area(data = sub[sub$eigenvector > 0, ], fill = "#E41A1C", alpha = 0.7) +
    geom_area(data = sub[sub$eigenvector <= 0, ], fill = "#377EB8", alpha = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(title = ct_name, x = paste(chrom, "(Mb)"), y = "Eigenvector") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

for (chrom in c("chr1", "chr14", "chr19")) {
  plots <- list()
  for (ct in names(compartment_results)) {
    p <- plot_track(all_comp, ct, chrom)
    if (!is.null(p)) plots[[ct]] <- p
  }
  if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = 1)
    save_plot(combined, file.path(output_dir, paste0("Track_", chrom)),
              14, 2.5 * length(plots))
  }
}

# ---- 8b. Compartment  ----
compute_switch <- function(dt, ct1, ct2) {
  d1 <- dt[cell_type == ct1, .(chr, start, comp1 = compartment, e1 = eigenvector)]
  d2 <- dt[cell_type == ct2, .(chr, start, comp2 = compartment, e2 = eigenvector)]
  m <- merge(d1, d2, by = c("chr", "start"))
  m[, switch := fifelse(comp1 == "open" & comp2 == "closed", "A->B",
                        fifelse(comp1 == "closed" & comp2 == "open", "B->A", "stable"))]
  m[, comparison := paste0(ct1, " vs ", ct2)]
  m
}

comparisons <- list(
  c("Bipotent Progenitor", "Hepatocytes"),
  c("Bipotent Progenitor", "CSC"),
  c("Hepatocytes",         "CSC"),
  c("CSC",                 "CSC(progenitor_cycling)"),
  c("CSC",                 "CSC(progenitor_quiet)"),
  c("CSC(progenitor_cycling)", "CSC(progenitor_quiet)"),
  c("CSC",                 "Hepatocytes (Metabolic-L)")
)

switch_all <- list()
for (pair in comparisons) {
  if (all(pair %in% names(compartment_results))) {
    switch_all[[paste(pair, collapse = " vs ")]] <- compute_switch(all_comp, pair[1], pair[2])
  }
}

if (length(switch_all) > 0) {
  switch_dt <- rbindlist(switch_all)
  fwrite(switch_dt, file.path(output_dir, "Compartment_Switches.csv"))

  sw_summ <- switch_dt[, .N, by = .(comparison, switch)]
  p_sw <- ggplot(sw_summ, aes(x = comparison, y = N, fill = switch)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = c("A->B" = "#E41A1C", "B->A" = "#377EB8", "stable" = "grey70")) +
    labs(title = "Compartment Switches Between Cell Types",
         x = "", y = "Number of 1Mb bins", fill = "") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_sw, file.path(output_dir, "Switch_Summary"), 12, 6)

  for (comp_name in names(switch_all)) {
    sw <- switch_all[[comp_name]]
    pair_names <- strsplit(comp_name, " vs ")[[1]]

    p <- ggplot(sw, aes(x = e1, y = e2, color = switch)) +
      geom_point(alpha = 0.4, size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, color = "grey80") +
      geom_vline(xintercept = 0, color = "grey80") +
      scale_color_manual(values = c("A->B" = "#E41A1C", "B->A" = "#377EB8", "stable" = "grey70")) +
      labs(title = comp_name,
           x = paste("Eigenvector:", pair_names[1]),
           y = paste("Eigenvector:", pair_names[2])) +
      theme_bw(base_size = 12)
    fname <- gsub("[^A-Za-z0-9_]", "_", comp_name)
    save_plot(p, file.path(output_dir, paste0("Scatter_", fname)), 7, 6)
  }
}

# ---- 8c. A/B  ----
pct_dt <- all_comp[, .(
  A_pct = 100 * sum(compartment == "open") / .N,
  B_pct = 100 * sum(compartment == "closed") / .N
), by = cell_type]

pct_long <- melt(pct_dt, id.vars = "cell_type", variable.name = "Comp", value.name = "Pct")
pct_long[, Comp := fifelse(Comp == "A_pct", "A (Open)", "B (Closed)")]

ct_order <- c("Bipotent Progenitor", "Hepatocytes", "Hepatocytes (Metabolic-L)",
              "CSC", "CSC(progenitor_cycling)", "CSC(progenitor_quiet)")
pct_long$cell_type <- factor(pct_long$cell_type,
                             levels = ct_order[ct_order %in% pct_long$cell_type])

p_pct <- ggplot(pct_long, aes(x = cell_type, y = Pct, fill = Comp)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("A (Open)" = "#377EB8", "B (Closed)" = "#E41A1C")) +
  labs(title = "A/B Compartment Proportion by Cell Type",
       x = "", y = "Proportion (%)", fill = "") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_pct, file.path(output_dir, "Proportion"), 9, 5)

key_pairs <- list(
  c("Bipotent Progenitor", "CSC"),
  c("Hepatocytes", "CSC"),
  c("Bipotent Progenitor", "Hepatocytes")
)

for (pair in key_pairs) {
  if (!all(pair %in% names(compartment_results))) next

  for (chrom in c("chr1", "chr14")) {
    d1 <- all_comp[cell_type == pair[1] & chr == chrom]
    d2 <- all_comp[cell_type == pair[2] & chr == chrom]
    mg <- merge(d1[, .(chr, start, mid, e1 = eigenvector)],
                d2[, .(chr, start, e2 = eigenvector)], by = c("chr", "start"))
    mg[, delta := e2 - e1]
    mg_df <- as.data.frame(mg)

    p1 <- plot_track(all_comp, pair[1], chrom)
    p2 <- plot_track(all_comp, pair[2], chrom)
    p3 <- ggplot(mg_df, aes(x = mid / 1e6, y = delta)) +
      geom_area(data = mg_df[mg_df$delta > 0, ], fill = "#E41A1C", alpha = 0.6) +
      geom_area(data = mg_df[mg_df$delta <= 0, ], fill = "#377EB8", alpha = 0.6) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      labs(title = paste0("\u0394 (", pair[2], " \u2212 ", pair[1], ")"),
           x = paste(chrom, "(Mb)"), y = "\u0394 Eigenvector") +
      theme_classic(base_size = 10)

    fname <- paste0("Diff_", gsub("[^A-Za-z0-9]", "_", pair[1]), "_vs_",
                    gsub("[^A-Za-z0-9]", "_", pair[2]), "_", chrom)
    save_plot(p1 / p2 / p3, file.path(output_dir, fname), 14, 8)
  }
}

# ======================== 9. Domain inflections ========================

for (ct in names(compartment_results)) {
  gr <- compartment_results[[ct]]
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)

  inflections <- tryCatch({
    getDomainInflections(gr, what = "flip.score", res = resolution,
                        chrs = chromosomes, genome = "hg38")
  }, error = function(e) {
    tryCatch(
      getDomainInflections(gr, what = "score", res = resolution,
                           chrs = chromosomes, genome = "hg38"),
      error = function(e2) NULL
    )
  })

  if (!is.null(inflections)) {
    saveRDS(inflections, file.path(output_dir, paste0("Inflections_", ct_safe, ".rds")))
    message("    ", ct, ": ", length(inflections), " inflection points")
  }
}

# ======================== 10.  ========================
cat("\n", strrep("=", 60), "\n")
cat(strrep("=", 60), "\n\n")

print(as.data.frame(pct_dt))

if (length(switch_all) > 0) {
  sw_summary <- switch_dt[, .N, by = .(comparison, switch)]
  sw_wide <- dcast(sw_summary, comparison ~ switch, value.var = "N", fill = 0)
  print(as.data.frame(sw_wide))
}

ev_summary <- all_comp[, .(
  mean_eigenvector = round(mean(eigenvector, na.rm = TRUE), 5),
  n_bins = .N
), by = cell_type]
print(as.data.frame(ev_summary))

# --- Supplementary: Chr14 compartment interaction heatmap ---

# ======================== 7.  ( getDenoisedCorMatrix) ========================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (ct in cell_types) {
  if (!(ct %in% names(compartment_results))) next
  
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  
  cells <- colnames(seu_obj)[seu_obj$final_type == ct]
  
  if (length(cells) < 10) {
    next
  }
  
  if (length(cells) > max_cells) { 
    set.seed(42)
    cells <- sample(cells, max_cells) 
  }
  
  counts_mat <- GetAssayData(seu_obj, layer = "counts")[matched_genes, cells]
  if (inherits(counts_mat, "sparseMatrix")) counts_mat <- as.matrix(counts_mat)
  
  rse_tmp <- SummarizedExperiment(
    assays    = list(counts = counts_mat),
    rowRanges = gene_gr[matched_genes],
    colData   = DataFrame(row.names = cells)
  )
  colnames(rse_tmp) <- cells
  tfidf_tmp <- transformTFIDF(assay(rse_tmp, "counts"))
  assay(rse_tmp, "counts") <- t(tfidf_tmp)
  
  rmt <- NULL
  tryCatch({
    rmt <- getDenoisedCorMatrix(
      rse_tmp, 
      res    = resolution, 
      chr    = "chr14",
      genome = "hg38", 
      assay  = "rna", 
      iter   = 2
    )
    
    if (is.null(rmt) || all(is.na(rmt))) {
      rmt <- NULL
    } else if (var(as.vector(rmt), na.rm = TRUE) == 0) {
      rmt <- NULL
    }
  }, error = function(e) {
  })
  
  pdf(file.path(output_dir, paste0("CorMatrix_", ct_safe, "_chr14.pdf")),
      width = 6, height = 6)
  
  if (!is.null(rmt)) {
    tryCatch({
      p <- plotCorMatrix(rmt, uppertri = TRUE)
      
      if (!is.null(p)) {
        print(p)
      }
      
      grid::grid.text(paste0(ct, " - chr14"), 
                      x = 0.5, y = 0.98, 
                      gp = grid::gpar(fontsize = 12, fontface = "bold"))
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste0("Plot Failed: ", e$message), cex = 0.8)
    })
  } else {
    plot.new()
    text(0.5, 0.5, paste0("No valid correlation matrix\nfor ", ct), cex = 0.9)
  }
  
  dev.off()
}

# --- Supplementary: CSC epigenetic full analysis (TriOmics + DMC + plasticity) ---

#!/usr/bin/env Rscript
#
#  TriOmics (RNA + ATAC + Methylation) + DMC Volcano + Genomic Annotation
#  + Plasticity + PCA + Imprinting + Repeat Elements + Motif
#
#
# PART A —  (scRNA x scATAC x Methylation)
#      Fig_Venn_CSC_UP/DOWN_TriOmics.pdf
#      Fig_Scatter_RNA_vs_ATAC.pdf
#      Fig_Scatter_RNA_vs_Meth.pdf
#      Fig_Quadrant_Epigenetic.pdf
#      Fig_Heatmap_TriOmics.pdf
#      Fig_Bar_Concordance.pdf
#
# PART B — DMC  (FP vs FN)
#      Fig4C3_DMC_Volcano.pdf
#      Fig4C4_DMC_GenomicAnnotation.pdf
#
#      Fig4C5a_Plasticity_Maintain.pdf
#      Fig4C5b_Plasticity_Acquire.pdf
#      Fig4C6_Meth_Memory_Heatmap.pdf
#      Fig4C7_PCA_6Groups.pdf
#      Fig_Meth_Consistency_*.pdf
#
#      Fig_Imprinted_*.pdf
#      Fig_Repeat_*.pdf
#      Fig_Motif_*.pdf
#

# ----  AnnotationDbi::select  dplyr::select  ----
library(data.table)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(VennDiagram)
library(grid)
library(gridExtra)
library(Seurat)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# select/filter/rename  dplyr
select  <- dplyr::select
filter  <- dplyr::filter
rename  <- dplyr::rename
mutate  <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by

cat("\n"); message(strrep("#", 70))
message("#  CSC Epigenetic Integration Pipeline (Consolidated)")
message(strrep("#", 70))

# 0. PATH CONFIG
SEURAT_PATH    <- "/data2/adata/atlas.RData"
ATAC_UP_PATH   <- "/data2/adata/sctour_result_hep/export_for_r_stem/CSC_Up_DA_peaks_Annotated.csv"
ATAC_DOWN_PATH <- "/data2/adata/sctour_result_hep/export_for_r_stem/CSC_Down_DA_peaks_Annotated.csv"

DMC_DIR   <- "/home/download/scTE/dma/"
DMR_DIR   <- "/home/download/scTE/Batch_DMR_Final_v5/"
INPUT_DIR <- "/home/download/csc_article/methylation/input/"

WGCNA_DIR <- "/home/download/csc_article/WGCNA_V2/"
OUT_DIR   <- "/home/download/csc_article/fig5/Epigenetic_Integration/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCATTER_W <- 5
SCATTER_H <- 5

METH_COMPARISONS <- list(
  list(name = "FN_vs_FP",   dmc = "FN_vs_FP",   dmr = "N_B_vs_P_B",
       desc = "CSC vs non-CSC (baseline)"),
  list(name = "N_N_vs_N_P", dmc = "N_N_vs_N_P", dmr = "N_N_vs_N_P",
       desc = "non-CSC: maintain vs acquire CSC phenotype (2W)"),
  list(name = "P_N_vs_P_P", dmc = "P_N_vs_P_P", dmr = "P_N_vs_P_P",
       desc = "CSC: lose vs maintain stemness (2W)")
)

DNMT_GENES <- c("DNMT1", "DNMT3A", "DNMT3B", "DNMT3L",
                "TET1", "TET2", "TET3", "UHRF1", "UHRF2")

# ---- Helper: gene annotation () ----
get_gene_symbols <- function(entrez_ids) {
  if (length(entrez_ids) == 0) return(data.frame(ENTREZID = character(), SYMBOL = character()))
  tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = unique(entrez_ids),
                          columns = "SYMBOL", keytype = "ENTREZID"),
    error = function(e) data.frame(ENTREZID = character(), SYMBOL = character())
  )
}

# ---- Helper: annotate methylation sites to genes ----
annotate_meth_to_genes <- function(meth_df) {
  if (!"chr" %in% colnames(meth_df) || nrow(meth_df) == 0) {
    return(data.frame(SYMBOL = character(), meth_diff = numeric()))
  }
  gr <- makeGRangesFromDataFrame(meth_df, keep.extra.columns = TRUE)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- trim(gr)

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  prom_gr <- suppressMessages(
    promoters(genes(txdb, single.strand.genes.only = FALSE), upstream = 3000, downstream = 3000)
  )
  # unlist if GRangesList
  if (is(prom_gr, "GRangesList")) { prom_gr <- unlist(prom_gr); prom_gr <- prom_gr[!duplicated(names(prom_gr))] }

  hits <- findOverlaps(gr, prom_gr)
  if (length(hits) == 0) return(data.frame(SYMBOL = character(), meth_diff = numeric()))

  result <- data.frame(idx = queryHits(hits), entrez = names(prom_gr)[subjectHits(hits)],
                       stringsAsFactors = FALSE)
  sym_map <- get_gene_symbols(unique(result$entrez))
  result <- merge(result, sym_map, by.x = "entrez", by.y = "ENTREZID", all.x = TRUE)

  if ("meth.diff" %in% colnames(meth_df)) {
    # Map idx back to original meth_df rows (after keepStandardChromosomes filtering)
    # Need to track which rows survived filtering
    gr_full <- makeGRangesFromDataFrame(meth_df, keep.extra.columns = TRUE)
    std_idx <- which(as.character(seqnames(gr_full)) %in% seqlevels(keepStandardChromosomes(gr_full, pruning.mode = "coarse")))
    result$meth_diff <- meth_df$meth.diff[std_idx[result$idx]]
  }

  result <- result[!is.na(result$SYMBOL), ]
  result %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(meth_diff = mean(meth_diff, na.rm = TRUE), .groups = "drop")
}

# ---- Helper: sample group ----
get_group <- function(id) {
  if (grepl("^P_\\d+_B$", id)) return("FP")
  if (grepl("^N_\\d+_B$", id)) return("FN")
  if (grepl("^P_P_\\d+$", id)) return("PP")
  if (grepl("^P_N_\\d+$", id)) return("PN")
  if (grepl("^N_\\d+_P$", id)) return("NP")
  if (grepl("^N_\\d+_N$", id)) return("NN")
  return("Other")
}

# ---- Helper: scatter plot ----
plot_epi_scatter <- function(df, xvar, yvar, xlabel, ylabel,
                              color = "#D62728", title = "", color_var = NULL) {
  ok <- is.finite(df[[xvar]]) & is.finite(df[[yvar]])
  ct <- cor.test(df[[xvar]][ok], df[[yvar]][ok], method = "spearman")

  p <- ggplot(df[ok, ], aes(.data[[xvar]], .data[[yvar]])) +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey70", linetype = "dashed")

  if (!is.null(color_var) && color_var %in% colnames(df)) {
    p <- p + geom_point(aes(color = .data[[color_var]]), alpha = 0.4, size = 1.5) +
      scale_color_manual(values = c("UP" = "#D62728", "DOWN" = "#1F77B4"))
  } else {
    p <- p + geom_point(alpha = 0.3, size = 1.5, color = color)
  }

  p + geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
    labs(title = title,
         subtitle = sprintf("R=%.3f  P=%.2e  (n=%d)", ct$estimate, ct$p.value, sum(ok)),
         x = xlabel, y = ylabel) +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1)
}

###############################################################################
#                                                                             #
# PART A: scRNA DEG + scATAC + Methylation                         #
#                                                                             #
###############################################################################

message("\n", strrep("=", 60))
message("  PART A: TriOmics Integration")
message(strrep("=", 60))

# ---- A1: CSC DEG ----
message("\n--- A1: CSC DEG ---")
deg_path <- file.path(OUT_DIR, "CSC_vs_OtherTumor_DEG.csv")

if (file.exists(deg_path)) {
  message("  Loading cached DEG")
  deg_df <- read.csv(deg_path, row.names = 1)
} else {
  message("  Computing DEG from Seurat...")
  load(SEURAT_PATH)
  seu_rna <- NULL
  for (vn in ls()) {
    obj <- tryCatch(get(vn), error = function(e) NULL)
    if (inherits(obj, "Seurat")) { seu_rna <- obj; message("  Found: ", vn); break }
  }
  if (is.null(seu_rna)) stop("No Seurat object found")

  ct_col <- ifelse("final_type" %in% colnames(seu_rna@meta.data), "final_type", "cell_type")
  if ("malignant" %in% colnames(seu_rna@meta.data)) {
    seu_rna <- subset(seu_rna, malignant == "tumor")
    message("  Filtered to tumor: ", ncol(seu_rna))
  }
  all_types   <- unique(seu_rna@meta.data[[ct_col]])
  csc_types   <- all_types[grepl("^CSC", all_types)]
  tumor_types <- all_types[grepl("Hepatocyte|Hep", all_types)]

  keep <- colnames(seu_rna)[seu_rna@meta.data[[ct_col]] %in% c(csc_types, tumor_types)]
  seu_sub <- subset(seu_rna, cells = keep)
  DefaultAssay(seu_sub) <- "RNA"
  if ("raw_counts" %in% Layers(seu_sub)) seu_sub[["RNA"]]$counts <- seu_sub[["RNA"]]$raw_counts
  seu_sub <- NormalizeData(seu_sub)
  seu_sub$is_CSC <- ifelse(seu_sub@meta.data[[ct_col]] %in% csc_types, "CSC", "Other_Tumor")
  Idents(seu_sub) <- "is_CSC"

  deg_df <- FindMarkers(seu_sub, ident.1 = "CSC", ident.2 = "Other_Tumor",
                        min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
  deg_df$gene <- rownames(deg_df)
  write.csv(deg_df, deg_path)
  rm(seu_rna, seu_sub); gc()
}
if (!"gene" %in% colnames(deg_df)) deg_df$gene <- rownames(deg_df)
deg_up   <- deg_df %>% dplyr::filter(avg_log2FC > 0.5, p_val_adj < 0.05)
deg_down <- deg_df %>% dplyr::filter(avg_log2FC < -0.5, p_val_adj < 0.05)
message("  UP: ", nrow(deg_up), " | DOWN: ", nrow(deg_down))

# ---- A2: scATAC DA ----
message("\n--- A2: scATAC DA ---")
atac_genes_up <- atac_genes_down <- character(0)
atac_fc <- data.frame()

if (file.exists(ATAC_UP_PATH)) {
  atac_up <- read.csv(ATAC_UP_PATH)
  atac_genes_up <- unique(na.omit(atac_up$gene_name))
  if ("avg_log2FC" %in% colnames(atac_up)) {
    atac_fc_up <- atac_up %>% dplyr::filter(!is.na(gene_name)) %>%
      dplyr::group_by(gene_name) %>% dplyr::summarise(atac_log2FC = max(avg_log2FC), .groups = "drop")
  }
  message("  ATAC UP genes: ", length(atac_genes_up))
}
if (file.exists(ATAC_DOWN_PATH)) {
  atac_down <- read.csv(ATAC_DOWN_PATH)
  atac_genes_down <- unique(na.omit(atac_down$gene_name))
  if ("avg_log2FC" %in% colnames(atac_down)) {
    atac_fc_down <- atac_down %>% dplyr::filter(!is.na(gene_name)) %>%
      dplyr::group_by(gene_name) %>% dplyr::summarise(atac_log2FC = min(avg_log2FC), .groups = "drop")
  }
  message("  ATAC DOWN genes: ", length(atac_genes_down))
}
if (exists("atac_fc_up") && exists("atac_fc_down")) {
  atac_fc <- dplyr::bind_rows(atac_fc_up, atac_fc_down) %>%
    dplyr::group_by(gene_name) %>% dplyr::summarise(atac_log2FC = mean(atac_log2FC), .groups = "drop")
}

# ---- A3: Methylation (3 comparisons) ----
message("\n--- A3: Methylation ---")
meth_genes_hyper <- meth_genes_hypo <- character(0)
meth_fc <- data.frame()
all_meth_results <- list()

for (comp in METH_COMPARISONS) {
  message("  ", comp$name, ": ", comp$desc)
  dmc_path <- file.path(DMC_DIR, comp$dmc, "Significant_DMCs_all.csv")
  dmr_path <- file.path(DMR_DIR, comp$dmr, "DMR_Significant.csv")
  mpath <- if (file.exists(dmc_path)) dmc_path else if (file.exists(dmr_path)) dmr_path else NULL
  if (is.null(mpath)) { message("    Not found, skip"); next }

  mdf <- fread(mpath)
  meth_annot <- annotate_meth_to_genes(as.data.frame(mdf))
  message("    ", nrow(meth_annot), " genes")
  if (nrow(meth_annot) > 0)
    all_meth_results[[comp$name]] <- meth_annot %>% dplyr::mutate(comparison = comp$name)

  if (comp$name == "FN_vs_FP" && nrow(meth_annot) > 0) {
    meth_genes_hyper <- meth_annot$SYMBOL[meth_annot$meth_diff > 0]
    meth_genes_hypo  <- meth_annot$SYMBOL[meth_annot$meth_diff < 0]
    meth_fc <- meth_annot %>% dplyr::rename(gene_name = SYMBOL)
    message("    [Primary] Hyper: ", length(meth_genes_hyper), " | Hypo: ", length(meth_genes_hypo))
  }
}

# ---- A4: Integrate ----
message("\n--- A4: Build integrated table ---")
integrated <- deg_df %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  dplyr::select(gene, rna_log2FC = avg_log2FC, rna_padj = p_val_adj)

if (nrow(atac_fc) > 0) integrated <- dplyr::left_join(integrated, atac_fc, by = c("gene" = "gene_name"))
if (nrow(meth_fc) > 0) integrated <- dplyr::left_join(integrated, meth_fc, by = c("gene" = "gene_name"))

integrated <- integrated %>% dplyr::mutate(
  rna_dir  = ifelse(rna_log2FC > 0, "UP", "DOWN"),
  atac_dir = dplyr::case_when(is.na(atac_log2FC) ~ "NS", atac_log2FC > 0 ~ "Open", atac_log2FC < 0 ~ "Closed", TRUE ~ "NS"),
  meth_dir = dplyr::case_when(is.na(meth_diff) ~ "NS", meth_diff > 0 ~ "Hyper", meth_diff < 0 ~ "Hypo", TRUE ~ "NS"),
  concordant = dplyr::case_when(
    rna_dir == "UP" & atac_dir == "Open" & meth_dir == "Hypo" ~ "Full_Concordant_UP",
    rna_dir == "DOWN" & atac_dir == "Closed" & meth_dir == "Hyper" ~ "Full_Concordant_DOWN",
    rna_dir == "UP" & atac_dir == "Open" ~ "RNA+ATAC_Concordant",
    rna_dir == "DOWN" & atac_dir == "Closed" ~ "RNA+ATAC_Concordant",
    rna_dir == "UP" & meth_dir == "Hypo" ~ "RNA+Meth_Concordant",
    rna_dir == "DOWN" & meth_dir == "Hyper" ~ "RNA+Meth_Concordant",
    TRUE ~ "Other"
  )
)

hub_genes <- character(0)
hub_path <- file.path(WGCNA_DIR, "hub_genes_per_module.csv")
if (file.exists(hub_path)) {
  hub_df <- fread(hub_path)
  hub_genes <- unique(hub_df$gene)
  integrated$is_hub <- integrated$gene %in% hub_genes
  integrated$hub_module <- hub_df$module[match(integrated$gene, hub_df$gene)]
}

write.csv(integrated, file.path(OUT_DIR, "Table_TriOmics_Integrated.csv"), row.names = FALSE)
message("  Integrated: ", nrow(integrated), " genes")
print(table(integrated$concordant))

# ---- A5: Venn ----
message("\n--- A5: Venn ---")
venn_up <- list(`RNA UP` = deg_up$gene, `ATAC Open` = atac_genes_up, `Meth Hypo` = meth_genes_hypo)
venn_up <- venn_up[sapply(venn_up, length) > 0]
if (length(venn_up) >= 2) {
  pdf(file.path(OUT_DIR, "Fig_Venn_CSC_UP_TriOmics.pdf"), width = 6, height = 6)
  grid.draw(venn.diagram(venn_up, filename = NULL,
    fill = c("#E74C3C", "#3498DB", "#2ECC71")[1:length(venn_up)], alpha = 0.5,
    cat.fontsize = 12, main = "CSC Upregulated Genes\n(RNA UP & ATAC Open & Meth Hypo)", main.fontsize = 13))
  dev.off(); message("  -> Fig_Venn_CSC_UP")
}
venn_dn <- list(`RNA DOWN` = deg_down$gene, `ATAC Closed` = atac_genes_down, `Meth Hyper` = meth_genes_hyper)
venn_dn <- venn_dn[sapply(venn_dn, length) > 0]
if (length(venn_dn) >= 2) {
  pdf(file.path(OUT_DIR, "Fig_Venn_CSC_DOWN_TriOmics.pdf"), width = 6, height = 6)
  grid.draw(venn.diagram(venn_dn, filename = NULL,
    fill = c("#9B59B6", "#E67E22", "#1ABC9C")[1:length(venn_dn)], alpha = 0.5,
    cat.fontsize = 12, main = "CSC Downregulated Genes\n(RNA DOWN & ATAC Closed & Meth Hyper)", main.fontsize = 13))
  dev.off(); message("  -> Fig_Venn_CSC_DOWN")
}

# ---- A6: Scatter RNA vs ATAC ----
if ("atac_log2FC" %in% colnames(integrated)) {
  df_p <- integrated %>% dplyr::filter(!is.na(atac_log2FC))
  if (nrow(df_p) > 5) {
    p <- plot_epi_scatter(df_p, "rna_log2FC", "atac_log2FC",
      "RNA log2FC (CSC vs Other Tumor)", "ATAC log2FC (CSC vs Other)",
      title = "CSC: Transcription vs Chromatin Accessibility", color_var = "rna_dir")
    ggsave(file.path(OUT_DIR, "Fig_Scatter_RNA_vs_ATAC.pdf"), p, width = SCATTER_W, height = SCATTER_H)
    message("  -> Fig_Scatter_RNA_vs_ATAC")
  }
}

# ---- A7: Scatter RNA vs Meth ----
if ("meth_diff" %in% colnames(integrated)) {
  df_p <- integrated %>% dplyr::filter(!is.na(meth_diff))
  if (nrow(df_p) > 5) {
    p <- plot_epi_scatter(df_p, "rna_log2FC", "meth_diff",
      "RNA log2FC (CSC vs Other Tumor)", "Methylation Difference (%)",
      title = "CSC: Transcription vs Promoter Methylation", color_var = "rna_dir")
    ggsave(file.path(OUT_DIR, "Fig_Scatter_RNA_vs_Meth.pdf"), p, width = SCATTER_W, height = SCATTER_H)
    message("  -> Fig_Scatter_RNA_vs_Meth")
  }
}

# ---- A8: Quadrant ATAC x Meth ----
has_both <- integrated %>% dplyr::filter(!is.na(atac_log2FC), !is.na(meth_diff))
if (nrow(has_both) > 10) {
  p <- ggplot(has_both, aes(atac_log2FC, meth_diff)) +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey70", linetype = "dashed") +
    geom_point(aes(color = rna_dir, size = abs(rna_log2FC)), alpha = 0.5) +
    scale_color_manual(values = c("UP" = "#D62728", "DOWN" = "#1F77B4")) +
    scale_size_continuous(range = c(1, 4), name = "|RNA log2FC|") +
    annotate("text", x = max(has_both$atac_log2FC)*0.7, y = min(has_both$meth_diff)*0.7,
             label = "Active\n(Open+Hypo)", color = "#D62728", fontface = "bold", size = 3.5) +
    annotate("text", x = min(has_both$atac_log2FC)*0.7, y = max(has_both$meth_diff)*0.7,
             label = "Silenced\n(Closed+Hyper)", color = "#1F77B4", fontface = "bold", size = 3.5) +
    labs(title = "Epigenetic Regulation Landscape (CSC)",
         x = "ATAC log2FC", y = "Methylation Difference (%)", color = "RNA") +
    theme_bw(base_size = 12) + theme(aspect.ratio = 1)
  ggsave(file.path(OUT_DIR, "Fig_Quadrant_Epigenetic.pdf"), p, width = 6, height = 6)
  message("  -> Fig_Quadrant_Epigenetic")
}

# ---- A9: Heatmap concordant genes ----
conc_genes <- integrated %>%
  dplyr::filter(concordant != "Other") %>%
  dplyr::arrange(desc(abs(rna_log2FC))) %>% head(50)
if (nrow(conc_genes) >= 5) {
  mat_cols <- intersect(c("rna_log2FC", "atac_log2FC", "meth_diff"), colnames(conc_genes))
  mat <- conc_genes %>% dplyr::select(gene, dplyr::all_of(mat_cols)) %>%
    tibble::column_to_rownames("gene") %>% as.matrix()
  mat_sc <- apply(mat, 2, function(x) { s <- sd(x, na.rm=T); if(is.na(s)||s==0) x else (x-mean(x,na.rm=T))/s })
  ann <- data.frame(Direction = conc_genes$rna_dir, row.names = conc_genes$gene)
  if ("hub_module" %in% colnames(conc_genes)) ann$Hub <- ifelse(is.na(conc_genes$hub_module), "", conc_genes$hub_module)
  col_lbl <- c(rna_log2FC="RNA\nlog2FC", atac_log2FC="ATAC\nlog2FC", meth_diff="Meth\nDiff(%)")
  pdf(file.path(OUT_DIR, "Fig_Heatmap_TriOmics.pdf"), width = 5, height = max(6, nrow(mat_sc)*0.25))
  pheatmap(mat_sc, cluster_cols=FALSE, labels_col=col_lbl[colnames(mat_sc)], annotation_row=ann,
           color=colorRampPalette(c("#2166AC","white","#B2182B"))(100), border_color=NA, fontsize_row=7,
           main="Epigenetically Concordant CSC Genes")
  dev.off(); message("  -> Fig_Heatmap_TriOmics")
}

# ---- A10: Bar concordance ----
conc_summary <- integrated %>% dplyr::count(concordant) %>% dplyr::arrange(desc(n)) %>%
  dplyr::mutate(concordant = factor(concordant, levels = concordant))
p <- ggplot(conc_summary, aes(concordant, n, fill = concordant)) +
  geom_col(alpha = 0.85) + geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Epigenetic Concordance with CSC DEGs", x = "", y = "Number of Genes") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(OUT_DIR, "Fig_Bar_Concordance.pdf"), p, width = 7, height = 5)

# ---- A11: Methylation cross-comparison scatter ----
if (length(all_meth_results) >= 2) {
  pairs <- list(
    list("FN_vs_FP", "P_N_vs_P_P", "Baseline CSC vs 2W Maintained CSC"),
    list("FN_vs_FP", "N_N_vs_N_P", "Baseline CSC vs 2W Acquired CSC")
  )
  for (pr in pairs) {
    if (!pr[[1]] %in% names(all_meth_results) || !pr[[2]] %in% names(all_meth_results)) next
    shared <- dplyr::inner_join(all_meth_results[[pr[[1]]]], all_meth_results[[pr[[2]]]], by = "SYMBOL", suffix = c("_a", "_b"))
    if (nrow(shared) < 10) next
    p <- plot_epi_scatter(shared, "meth_diff_a", "meth_diff_b",
      paste0("Meth Diff (%) - ", pr[[1]]), paste0("Meth Diff (%) - ", pr[[2]]),
      color = "#8E44AD", title = pr[[3]])
    ggsave(file.path(OUT_DIR, sprintf("Fig_Meth_Consistency_%s_vs_%s.pdf", pr[[1]], pr[[2]])),
           p, width = SCATTER_W, height = SCATTER_H)
    message("  -> Meth consistency: ", pr[[3]])
  }
}

# ---- A12: 3-comparison heatmap ----
if (length(all_meth_results) >= 2) {
  wide <- all_meth_results[[1]] %>% dplyr::select(SYMBOL, meth_diff)
  colnames(wide)[2] <- names(all_meth_results)[1]
  for (i in 2:length(all_meth_results)) {
    tmp <- all_meth_results[[i]] %>% dplyr::select(SYMBOL, meth_diff)
    colnames(tmp)[2] <- names(all_meth_results)[i]
    wide <- dplyr::full_join(wide, tmp, by = "SYMBOL")
  }
  wide$n_present <- rowSums(!is.na(wide[, -1]))
  ws <- wide %>% dplyr::filter(n_present >= 2) %>% dplyr::select(-n_present)
  if (nrow(ws) > 20) {
    ws_top <- ws %>% dplyr::arrange(desc(abs(.data[[names(all_meth_results)[1]]]))) %>% head(50)
    mat <- ws_top %>% tibble::column_to_rownames("SYMBOL") %>% as.matrix()
    pdf(file.path(OUT_DIR, "Fig_Heatmap_Meth_3Comparisons.pdf"), width = 5, height = max(6, nrow(mat)*0.25))
    pheatmap(mat, cluster_cols = FALSE, color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
             border_color = NA, fontsize_row = 6, main = "Methylation Consistency across Comparisons")
    dev.off(); message("  -> Fig_Heatmap_Meth_3Comparisons")
  }
}

###############################################################################
#                                                                             #
#  PART B: DMC Volcano + Genomic Annotation                                   #
#                                                                             #
###############################################################################

message("\n", strrep("=", 60))
message("  PART B: DMC Volcano & Genomic Annotation")
message(strrep("=", 60))

# Load all tested sites
all_tested_path <- file.path(DMC_DIR, "FN_vs_FP", "All_tested_sites.csv")
sig_dmc_path    <- file.path(DMC_DIR, "FN_vs_FP", "Significant_DMCs_all.csv")
all_dmc <- if (file.exists(all_tested_path)) fread(all_tested_path) else fread(sig_dmc_path)

if (!"qvalue" %in% colnames(all_dmc)) all_dmc$qvalue <- all_dmc$p_val_adj
if (!"meth.diff" %in% colnames(all_dmc)) all_dmc$meth.diff <- all_dmc$meth_diff

all_dmc$neg_log10q <- -log10(pmax(all_dmc$qvalue, 1e-300))
all_dmc$sig_class <- "NS"
all_dmc$sig_class[all_dmc$meth.diff > 10 & all_dmc$qvalue < 0.05] <- "Hyper"
all_dmc$sig_class[all_dmc$meth.diff < -10 & all_dmc$qvalue < 0.05] <- "Hypo"
n_hyper <- sum(all_dmc$sig_class == "Hyper"); n_hypo <- sum(all_dmc$sig_class == "Hypo")
message("  Hyper: ", n_hyper, " | Hypo: ", n_hypo, " | Total: ", nrow(all_dmc))

message("  Annotating significant DMCs only...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- suppressMessages(genes(txdb))  # single.strand.genes.only=TRUE (default, faster)
promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 500)

sig_dmc_sub <- all_dmc[all_dmc$sig_class != "NS", ]
sig_gr <- GRanges(seqnames = sig_dmc_sub$chr,
                  ranges = IRanges(start = sig_dmc_sub$start, end = sig_dmc_sub$end))
sig_gr <- keepStandardChromosomes(sig_gr, pruning.mode = "coarse")
sig_gr <- trim(sig_gr)

# sig_dmc_sub  sig_gr
std_mask <- paste0(sig_dmc_sub$chr, ":", sig_dmc_sub$start) %in%
            paste0(seqnames(sig_gr), ":", start(sig_gr))
sig_dmc_sub <- sig_dmc_sub[std_mask, ]
sig_gr <- GRanges(seqnames = sig_dmc_sub$chr,
                  ranges = IRanges(start = sig_dmc_sub$start, end = sig_dmc_sub$end))

prom_hits <- findOverlaps(sig_gr, promoters_gr)
body_hits <- findOverlaps(sig_gr, genes_gr)

sig_dmc_sub$region <- "Intergenic"
sig_dmc_sub$region[unique(queryHits(body_hits))] <- "Gene Body"
sig_dmc_sub$region[unique(queryHits(prom_hits))] <- "Promoter"
message("  Promoter: ", sum(sig_dmc_sub$region == "Promoter"),
        " | Gene Body: ", sum(sig_dmc_sub$region == "Gene Body"),
        " | Intergenic: ", sum(sig_dmc_sub$region == "Intergenic"))

# Gene symbol for significant promoter DMCs (vectorized)
sig_dmc_sub$gene_symbol <- NA_character_
if (length(prom_hits) > 0) {
  pe <- names(promoters_gr)[subjectHits(prom_hits)]
  sm <- get_gene_symbols(unique(pe))
  e2s <- setNames(sm$SYMBOL, sm$ENTREZID)
  qi <- queryHits(prom_hits)
  sv <- e2s[pe]
  ok <- !is.na(sv) & !duplicated(qi)
  sig_dmc_sub$gene_symbol[qi[ok]] <- sv[ok]
  message("  Gene-annotated DMCs: ", sum(!is.na(sig_dmc_sub$gene_symbol)))
}

# ----  (DNMT + Hub + DEG top)  ----
message("  Targeted annotation for key genes on ALL CpG sites...")

key_genes <- unique(c(DNMT_GENES, hub_genes))
key_entrez <- tryCatch(
  AnnotationDbi::select(org.Hs.eg.db, keys = key_genes, columns = "ENTREZID", keytype = "SYMBOL"),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
)
key_entrez <- key_entrez[!is.na(key_entrez$ENTREZID), ]

if (nrow(key_entrez) > 0) {
  key_prom <- promoters_gr[names(promoters_gr) %in% key_entrez$ENTREZID]
  message("  Key gene promoters: ", length(key_prom))

  all_gr <- GRanges(seqnames = all_dmc$chr,
                    ranges = IRanges(start = all_dmc$start, end = all_dmc$end))
  key_hits <- findOverlaps(all_gr, key_prom)
  message("  Key gene hits in all CpGs: ", length(key_hits))

  all_dmc$gene_symbol <- NA_character_
  all_dmc$region <- NA_character_
  if (length(key_hits) > 0) {
    kqi <- queryHits(key_hits)
    kei <- names(key_prom)[subjectHits(key_hits)]
    ks2 <- setNames(key_entrez$SYMBOL, key_entrez$ENTREZID)
    ksv <- ks2[kei]
    kok <- !is.na(ksv) & !duplicated(kqi)
    all_dmc$gene_symbol[kqi[kok]] <- ksv[kok]
    all_dmc$region[kqi[kok]] <- "Promoter"
  }
}

# sig_dmc_sub  all_dmc
sig_key <- paste0(sig_dmc_sub$chr, ":", sig_dmc_sub$start)
all_key <- paste0(all_dmc$chr, ":", all_dmc$start)
match_idx <- match(sig_key, all_key)
valid_match <- !is.na(match_idx)
all_dmc$region[match_idx[valid_match]] <- sig_dmc_sub$region[valid_match]
need_sym <- valid_match & !is.na(sig_dmc_sub$gene_symbol) & is.na(all_dmc$gene_symbol[match_idx[valid_match]])
if (any(need_sym)) {
  all_dmc$gene_symbol[match_idx[need_sym]] <- sig_dmc_sub$gene_symbol[need_sym]
}

# Labels for volcano
all_dmc$label <- NA_character_
label_pool <- unique(c(DNMT_GENES, hub_genes))
# + promoter +
li <- which(all_dmc$gene_symbol %in% label_pool & all_dmc$sig_class != "NS")
if (length(li) > 0) {
  li <- li[order(all_dmc$qvalue[li])][1:min(30, length(li))]
  all_dmc$label[li] <- all_dmc$gene_symbol[li]
}
di <- which(all_dmc$gene_symbol %in% DNMT_GENES)
if (length(di) > 0) {
  di <- di[order(-abs(all_dmc$meth.diff[di]))][1:min(8, length(di))]
  all_dmc$label[di] <- all_dmc$gene_symbol[di]
}
message("  Labels: ", sum(!is.na(all_dmc$label)), " CpGs")

# ---- B1: Volcano ----
p_vol <- ggplot(all_dmc, aes(meth.diff, neg_log10q)) +
  geom_point(aes(color = sig_class), alpha = 0.4, size = 0.8) +
  geom_point(data = all_dmc[!is.na(all_dmc$label), ], aes(color = sig_class), size = 2.5, alpha = 0.9) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25, segment.size = 0.3,
                  fontface = "italic", na.rm = TRUE) +
  scale_color_manual(values = c(Hyper="#B2182B", Hypo="#2166AC", NS="grey70"),
    labels = c(sprintf("Hyper (n=%d)", n_hyper), sprintf("Hypo (n=%d)", n_hypo), "NS")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  labs(title = "DMC: CSC (FP) vs non-CSC (FN)", subtitle = "DNMT/Hub genes highlighted",
       x = "Methylation Difference (%)", y = expression(-log[10](FDR)), color = "") +
  theme_bw(base_size = 12) + theme(legend.position = c(0.85, 0.85))
ggsave(file.path(OUT_DIR, "Fig4C3_DMC_Volcano.pdf"), p_vol, width = 7, height = 6)
message("  -> Fig4C3_DMC_Volcano")

# ---- B2: Genomic Annotation ----
sig_dmc_sub <- all_dmc[all_dmc$sig_class != "NS", ]

# Stacked bar
annot_s <- sig_dmc_sub %>%
  dplyr::group_by(sig_class, region) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sig_class) %>% dplyr::mutate(pct = n / sum(n) * 100)

p_ann <- ggplot(annot_s, aes(sig_class, pct, fill = region)) +
  geom_col(position = "stack", alpha = 0.85, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)), position = position_stack(vjust = 0.5),
            size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = c(Promoter="#D62728", `Gene Body`="#FF7F0E", Intergenic="#1F77B4")) +
  labs(title = "Genomic Distribution of Significant DMCs", x = "", y = "Percentage (%)", fill = "Region") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "Fig4C4_DMC_GenomicAnnotation.pdf"), p_ann, width = 5.5, height = 5)
message("  -> Fig4C4_DMC_GenomicAnnotation")

# vs Background
bg_s <- all_dmc %>% dplyr::group_by(region) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(pct = n/sum(n)*100, group = "All CpGs")
h_s <- sig_dmc_sub %>% dplyr::filter(sig_class=="Hyper") %>% dplyr::group_by(region) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% dplyr::mutate(pct = n/sum(n)*100, group = "CSC Hyper-DMC")
l_s <- sig_dmc_sub %>% dplyr::filter(sig_class=="Hypo") %>% dplyr::group_by(region) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% dplyr::mutate(pct = n/sum(n)*100, group = "CSC Hypo-DMC")
comb <- dplyr::bind_rows(bg_s, h_s, l_s)
comb$group <- factor(comb$group, levels = c("All CpGs", "CSC Hyper-DMC", "CSC Hypo-DMC"))

p_ann2 <- ggplot(comb, aes(group, pct, fill = region)) +
  geom_col(position = "fill", alpha = 0.85, color = "white") +
  scale_fill_manual(values = c(Promoter="#D62728", `Gene Body`="#FF7F0E", Intergenic="#1F77B4")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "DMC Genomic Region vs Background", x = "", y = "Proportion", fill = "Region") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(OUT_DIR, "Fig4C4_DMC_Annotation_vsBackground.pdf"), p_ann2, width = 6, height = 5)
message("  -> Fig4C4_DMC_Annotation_vsBackground")

###############################################################################
#                                                                             #
#  PART C: Plasticity + Memory Heatmap + PCA                                  #
#                                                                             #
###############################################################################

message("\n", strrep("=", 60))
message("  PART C: CSC Methylation Plasticity")
message(strrep("=", 60))

# Load raw beta values
all_files  <- list.files(INPUT_DIR, pattern = "\\.cgMerged\\.935K\\.txt$", full.names = TRUE)
sample_ids <- gsub("\\.cgMerged\\.935K\\.txt$", "", basename(all_files))

if (length(all_files) == 0) {
  message("  WARNING: No raw methylation files in ", INPUT_DIR, " - skipping PART C")
} else {

  sample_grps <- sapply(sample_ids, get_group)
  message("  Samples: ", paste(paste0(sample_ids, "(", sample_grps, ")"), collapse = ", "))

  beta_list <- list()
  for (i in seq_along(all_files)) {
    dt <- fread(all_files[i], select = c(1, 2, 7, 8))
    colnames(dt) <- c("chr", "pos", "numC", "numT")
    dt[, coverage := numC + numT]
    dt <- dt[coverage >= 5]
    dt[, beta := numC / coverage]
    dt[, cpg_id := paste0(chr, ":", pos)]
    beta_list[[sample_ids[i]]] <- dt[, .(cpg_id, beta)]
    message("  ", sample_ids[i], ": ", nrow(dt))
  }

  # Merge
  beta_wide <- beta_list[[1]]; setnames(beta_wide, "beta", sample_ids[1])
  for (i in 2:length(beta_list)) {
    tmp <- beta_list[[i]]; setnames(tmp, "beta", sample_ids[i])
    beta_wide <- merge(beta_wide, tmp, by = "cpg_id", all = FALSE)
  }
  message("  Shared CpGs: ", nrow(beta_wide))

  scols <- sample_ids
  fp_s <- scols[sample_grps[scols] == "FP"]
  fn_s <- scols[sample_grps[scols] == "FN"]
  pp_s <- scols[sample_grps[scols] == "PP"]
  pn_s <- scols[sample_grps[scols] == "PN"]
  np_s <- scols[sample_grps[scols] == "NP"]
  nn_s <- scols[sample_grps[scols] == "NN"]

  # Significant DMC CpG IDs
  sig_cpg_ids <- paste0(sig_dmc_sub$chr, ":", sig_dmc_sub$start)
  sig_beta <- beta_wide[beta_wide$cpg_id %in% sig_cpg_ids, ]
  message("  Significant DMC sites with beta: ", nrow(sig_beta))

  compute_diff <- function(dt, g1, g2) {
    m1 <- rowMeans(as.matrix(dt[, ..g1, drop = FALSE]), na.rm = TRUE)
    m2 <- rowMeans(as.matrix(dt[, ..g2, drop = FALSE]), na.rm = TRUE)
    (m2 - m1) * 100
  }

  if (nrow(sig_beta) > 10) {
    baseline_diff <- compute_diff(sig_beta, fn_s, fp_s)

    # ---- C1: Plasticity - Maintain ----
    if (length(pp_s) > 0 && length(pn_s) > 0) {
      maintain_diff <- compute_diff(sig_beta, pn_s, pp_s)
      df_a <- data.frame(baseline = baseline_diff, maintain = maintain_diff)
      df_a <- df_a[is.finite(df_a$baseline) & is.finite(df_a$maintain), ]
      p <- plot_epi_scatter(df_a, "baseline", "maintain",
        "Baseline DMC: FP - FN (%)", "2W Maintained: PP - PN (%)",
        color = "#D62728", title = "CSC Methylation Memory: Maintain Stemness") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50")
      ggsave(file.path(OUT_DIR, "Fig4C5a_Plasticity_Maintain.pdf"), p, width = SCATTER_W, height = SCATTER_H)
      message("  -> Fig4C5a_Plasticity_Maintain")
    }

    # ---- C2: Plasticity - Acquire ----
    if (length(np_s) > 0 && length(nn_s) > 0) {
      acquire_diff <- compute_diff(sig_beta, nn_s, np_s)
      df_b <- data.frame(baseline = baseline_diff, acquire = acquire_diff)
      df_b <- df_b[is.finite(df_b$baseline) & is.finite(df_b$acquire), ]
      p <- plot_epi_scatter(df_b, "baseline", "acquire",
        "Baseline DMC: FP - FN (%)", "2W Acquired: NP - NN (%)",
        color = "#1F77B4", title = "CSC Methylation Acquisition: non-CSC gaining stemness") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50")
      ggsave(file.path(OUT_DIR, "Fig4C5b_Plasticity_Acquire.pdf"), p, width = SCATTER_W, height = SCATTER_H)
      message("  -> Fig4C5b_Plasticity_Acquire")
    }

    # ---- C3: Memory Heatmap ----
    mat_heat <- data.frame(cpg_id = sig_beta$cpg_id, Baseline = baseline_diff)
    if (length(pp_s) > 0 && length(pn_s) > 0) mat_heat$Maintain <- compute_diff(sig_beta, pn_s, pp_s)
    if (length(np_s) > 0 && length(nn_s) > 0) mat_heat$Acquire  <- compute_diff(sig_beta, nn_s, np_s)

    gene_map <- all_dmc[all_dmc$sig_class != "NS" & !is.na(all_dmc$gene_symbol), ]
    gene_map$cpg_id <- paste0(gene_map$chr, ":", gene_map$start)
    mat_heat <- merge(mat_heat, gene_map[, c("cpg_id", "gene_symbol")], by = "cpg_id", all.x = TRUE)
    mha <- mat_heat[!is.na(mat_heat$gene_symbol), ]
    mha <- mha[order(-abs(mha$Baseline)), ]
    mha <- mha[!duplicated(mha$gene_symbol), ]
    mha_top <- head(mha, 50)

    if (nrow(mha_top) >= 5) {
      mat_p <- mha_top[, grep("Baseline|Maintain|Acquire", colnames(mha_top))]
      rownames(mat_p) <- mha_top$gene_symbol
      pdf(file.path(OUT_DIR, "Fig4C6_Meth_Memory_Heatmap.pdf"), width = 5, height = max(6, nrow(mat_p)*0.28))
      pheatmap(as.matrix(mat_p), cluster_cols = FALSE, breaks = seq(-50, 50, length.out = 101),
               color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
               border_color = NA, fontsize_row = 7,
               labels_col = c("Baseline\nFP vs FN", "2W Maintain\nPP vs PN", "2W Acquire\nNP vs NN")[1:ncol(mat_p)],
               main = "CSC Methylation Memory (Top 50 Promoter DMCs)")
      dev.off(); message("  -> Fig4C6_Meth_Memory_Heatmap")
    }
  }

  # ---- C4: PCA ----
  message("\n--- PCA ---")
  pca_mat <- as.matrix(sig_beta[, ..scols])
  rownames(pca_mat) <- sig_beta$cpg_id
  pca_mat <- pca_mat[complete.cases(pca_mat), ]
  message("  PCA input: ", nrow(pca_mat), " x ", ncol(pca_mat))

  if (nrow(pca_mat) > 10 && ncol(pca_mat) >= 4) {
    pca_res <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
    pca_df <- data.frame(
      PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
      sample = rownames(pca_res$x), group = sample_grps[rownames(pca_res$x)])
    ve <- summary(pca_res)$importance[2, 1:2] * 100
    pca_df$group <- factor(pca_df$group, levels = c("FP","FN","PP","PN","NP","NN"))

    gcol <- c(FP="#D62728", FN="#1F77B4", PP="#FF9896", PN="#AEC7E8", NP="#FF7F0E", NN="#9EDAE5")
    gshp <- c(FP=16, FN=16, PP=17, PN=17, NP=15, NN=15)

    p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = group)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE) +
      scale_color_manual(values = gcol) + scale_shape_manual(values = gshp) +
      stat_ellipse(aes(group = group), type = "norm", level = 0.68, linetype = "dashed",
                   linewidth = 0.5, show.legend = FALSE) +
      labs(title = "PCA of CSC-associated DMCs",
           x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
      theme_bw(base_size = 12) + theme(aspect.ratio = 1)
    ggsave(file.path(OUT_DIR, "Fig4C7_PCA_6Groups.pdf"), p_pca, width = 7, height = 6)
    message("  -> Fig4C7_PCA_6Groups")
  }

  # ---- C5: Global beta + Entropy ----
  message("\n--- Global Beta & Entropy ---")

  global_stats <- data.frame()
  for (sid in scols) {
    bvals <- as.numeric(beta_wide[[sid]])
    bvals <- bvals[!is.na(bvals)]
    # Entropy: H = -[b*log2(b) + (1-b)*log2(1-b)]
    b_clip <- pmax(pmin(bvals, 0.999), 0.001)
    entropy <- -( b_clip * log2(b_clip) + (1 - b_clip) * log2(1 - b_clip) )
    global_stats <- rbind(global_stats, data.frame(
      sample = sid, group = sample_grps[sid],
      mean_beta = mean(bvals), median_beta = median(bvals),
      mean_entropy = mean(entropy), sd_beta = sd(bvals)
    ))
  }

  global_stats$group <- factor(global_stats$group, levels = c("FP","FN","PP","PN","NP","NN"))

  p_beta <- ggplot(global_stats, aes(group, mean_beta, fill = group)) +
    geom_col(alpha = 0.8, color = "grey30", linewidth = 0.3) +
    geom_point(size = 2) +
    scale_fill_manual(values = gcol) +
    labs(title = "Global Methylation Level", x = "", y = "Mean Beta") +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  p_ent <- ggplot(global_stats, aes(group, mean_entropy, fill = group)) +
    geom_col(alpha = 0.8, color = "grey30", linewidth = 0.3) +
    geom_point(size = 2) +
    scale_fill_manual(values = gcol) +
    labs(title = "Epigenetic Entropy", subtitle = "Lower = more ordered methylation state",
         x = "", y = "Mean Shannon Entropy") +
    theme_bw(base_size = 12) + theme(legend.position = "none")

  pdf(file.path(OUT_DIR, "Fig4C1_Global_Beta_Entropy.pdf"), width = 10, height = 4.5)
  gridExtra::grid.arrange(p_beta, p_ent, ncol = 2)
  dev.off()
  message("  -> Fig4C1_Global_Beta_Entropy")
  write.csv(global_stats, file.path(OUT_DIR, "Table_Global_Stats.csv"), row.names = FALSE)

} # end if raw files exist

###############################################################################
#                                                                             #
#  PART D: Imprinted Genes / Repeat Elements                                  #
#                                                                             #
###############################################################################

if (exists("beta_wide") && exists("scols")) {

message("\n", strrep("=", 60))
message("  PART D: Imprinted Genes & Repeat Elements")
message(strrep("=", 60))

# ---- D1: Imprinted Genes ----
imprinted_genes <- c(
  "IGF2","PEG3","PEG10","MEST","NDN","SNRPN","MAGEL2","MKRN3","DLK1","RTL1",
  "DIRAS3","PLAGL1","SGCE","KCNQ1OT1","NAP1L5","INPP5F","H19","CDKN1C","MEG3",
  "GRB10","UBE3A","PHLDA2","SLC22A18","TP73","GNAS","KCNQ1","IGF2R","WT1")

imp_entrez <- get_gene_symbols(imprinted_genes)  # reversed: SYMBOL -> ENTREZID won't work
imp_entrez <- tryCatch(
  AnnotationDbi::select(org.Hs.eg.db, keys = imprinted_genes, columns = "ENTREZID", keytype = "SYMBOL"),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
)
imp_entrez <- imp_entrez[!is.na(imp_entrez$ENTREZID), ]

if (nrow(imp_entrez) > 0) {
  imp_prom <- promoters(genes_gr[names(genes_gr) %in% imp_entrez$ENTREZID], upstream = 3000, downstream = 3000)

  # Ensure chr and pos columns exist in beta_wide
  if (!"chr" %in% colnames(beta_wide)) {
    parts <- tstrsplit(beta_wide$cpg_id, ":", fixed = TRUE)
    beta_wide[, chr := parts[[1]]]
    beta_wide[, pos := as.integer(parts[[2]])]
  }
  cpg_gr_all <- GRanges(seqnames = beta_wide$chr, ranges = IRanges(start = beta_wide$pos, width = 1))
  hits <- findOverlaps(cpg_gr_all, imp_prom)

  if (length(hits) > 0) {
    ent2sym <- setNames(imp_entrez$SYMBOL, imp_entrez$ENTREZID)
    imp_data <- data.frame()
    for (j in seq_along(queryHits(hits))) {
      qi <- queryHits(hits)[j]
      eid <- names(imp_prom)[subjectHits(hits)[j]]
      sym <- ent2sym[eid]
      if (is.na(sym)) next
      row_vals <- as.numeric(beta_wide[qi, ..scols])
      imp_data <- rbind(imp_data, data.frame(gene = sym, t(row_vals)))
    }
    if (nrow(imp_data) > 0) {
      colnames(imp_data) <- c("gene", scols)
      imp_mean <- imp_data %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(as.numeric(.), na.rm = TRUE)), .groups = "drop")

      imp_mean$FP_mean <- rowMeans(imp_mean[, fp_s, drop = FALSE], na.rm = TRUE)
      imp_mean$FN_mean <- rowMeans(imp_mean[, fn_s, drop = FALSE], na.rm = TRUE)
      imp_mean$diff <- imp_mean$FP_mean - imp_mean$FN_mean

      # Barplot
      imp_plot <- as.data.frame(imp_mean) %>%
        dplyr::filter(!is.na(diff)) %>% dplyr::arrange(diff) %>%
        dplyr::mutate(gene = factor(gene, levels = gene),
                      direction = ifelse(diff > 0, "CSC Hyper", "CSC Hypo"))

      p_imp <- ggplot(imp_plot, aes(gene, diff, fill = direction)) +
        geom_col(alpha = 0.85, color = "grey30", linewidth = 0.3) +
        scale_fill_manual(values = c("CSC Hyper" = "#B2182B", "CSC Hypo" = "#2166AC")) +
        coord_flip() +
        labs(title = "Imprinted Gene Methylation: CSC vs non-CSC", x = "", y = "Delta-beta (FP - FN)", fill = "") +
        theme_bw(base_size = 12) + theme(legend.position = "bottom")
      ggsave(file.path(OUT_DIR, "Fig_Imprinted_Barplot.pdf"), p_imp,
             width = 6, height = max(4, nrow(imp_plot)*0.25))
      message("  -> Fig_Imprinted_Barplot")
      write.csv(imp_mean, file.path(OUT_DIR, "Table_Imprinted_Methylation.csv"), row.names = FALSE)
    }
  }
}

# ---- D2: Repeat Elements ----
message("\n--- Repeat Elements ---")
rmsk_path <- file.path(OUT_DIR, "hg38_rmsk.txt.gz")
if (!file.exists(rmsk_path)) {
  message("  Downloading RepeatMasker...")
  tryCatch(
    download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz", rmsk_path, quiet = TRUE),
    error = function(e) message("  Download failed, skip repeat analysis"))
}

if (file.exists(rmsk_path)) {
  rmsk <- fread(rmsk_path, select = c(6, 7, 8, 12, 13))
  colnames(rmsk) <- c("chr", "start", "end", "repClass", "repFamily")

  repeat_classes <- list(LINE1 = rmsk[repFamily == "L1"], Alu = rmsk[repFamily == "Alu"],
                         LTR = rmsk[repClass == "LTR"], Satellite = rmsk[repClass == "Satellite"])

  if (!"chr" %in% colnames(beta_wide)) {
    parts <- tstrsplit(beta_wide$cpg_id, ":", fixed = TRUE)
    beta_wide[, chr := parts[[1]]]; beta_wide[, pos := as.integer(parts[[2]])]
  }
  cpg_gr_all <- GRanges(seqnames = beta_wide$chr, ranges = IRanges(start = beta_wide$pos, width = 1))

  rep_results <- list()
  for (rn in names(repeat_classes)) {
    rdt <- repeat_classes[[rn]]
    if (nrow(rdt) == 0) next
    r_gr <- GRanges(seqnames = rdt$chr, ranges = IRanges(start = rdt$start, end = rdt$end))
    hits <- findOverlaps(cpg_gr_all, r_gr)
    if (length(hits) == 0) next
    ridx <- unique(queryHits(hits))
    rb <- beta_wide[ridx, ..scols]
    per_s <- data.frame(sample = scols, beta = colMeans(as.matrix(rb), na.rm = TRUE),
                        group = sample_grps[scols], repeat_class = rn)
    rep_results[[rn]] <- per_s
    fp_v <- per_s$beta[per_s$group == "FP"]; fn_v <- per_s$beta[per_s$group == "FN"]
    wt <- tryCatch(wilcox.test(fp_v, fn_v), error = function(e) list(p.value = NA))
    message(sprintf("  %s: %d CpGs, FP=%.4f FN=%.4f P=%.3e", rn, length(ridx),
                    mean(fp_v), mean(fn_v), wt$p.value))
  }

  if (length(rep_results) > 0) {
    rep_df <- dplyr::bind_rows(rep_results)
    rep_df$group <- factor(rep_df$group, levels = c("FP","FN","PP","PN","NP","NN"))

    p_rep <- ggplot(rep_df, aes(group, beta, fill = group)) +
      geom_boxplot(alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = c(FP="#D62728", FN="#1F77B4", PP="#FF9896", PN="#AEC7E8", NP="#FF7F0E", NN="#9EDAE5")) +
      facet_wrap(~ repeat_class, scales = "free_y", nrow = 1) +
      labs(title = "Repeat Element Methylation across Groups", x = "", y = "Mean Beta") +
      theme_bw(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "Fig_Repeat_AllGroups.pdf"), p_rep, width = 14, height = 4.5)
    message("  -> Fig_Repeat_AllGroups")
    write.csv(rep_df, file.path(OUT_DIR, "Table_Repeat_Methylation.csv"), row.names = FALSE)
  }
}

} # end if beta_wide exists

###############################################################################
message("\n", strrep("#", 70))
message("#  ALL DONE!")
message(strrep("#", 70))
message(sprintf("\n  Output directory: %s\n", OUT_DIR))
message("  PART A - TriOmics:
    Fig_Venn_CSC_UP/DOWN_TriOmics.pdf
    Fig_Scatter_RNA_vs_ATAC.pdf / RNA_vs_Meth.pdf
    Fig_Quadrant_Epigenetic.pdf
    Fig_Heatmap_TriOmics.pdf / Fig_Bar_Concordance.pdf
    Fig_Meth_Consistency_*.pdf / Fig_Heatmap_Meth_3Comparisons.pdf

  PART B - DMC Baseline:
    Fig4C3_DMC_Volcano.pdf
    Fig4C4_DMC_GenomicAnnotation.pdf / _vsBackground.pdf

  PART C - Plasticity:
    Fig4C1_Global_Beta_Entropy.pdf
    Fig4C5a_Plasticity_Maintain.pdf / Fig4C5b_Plasticity_Acquire.pdf
    Fig4C6_Meth_Memory_Heatmap.pdf
    Fig4C7_PCA_6Groups.pdf

  PART D - Imprinting & Repeats:
    Fig_Imprinted_Barplot.pdf
    Fig_Repeat_AllGroups.pdf
")
