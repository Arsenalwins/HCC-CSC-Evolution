#!/usr/bin/env python
# Figure 3a,b,c — CNV score paired comparisons

import scanpy as sc
import pandas as pd
import numpy as np
import glob
import os
import itertools
import re
from scipy import stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs): return iterable

GLOBAL_PALETTE = {
    'CSC':                       '#D62728',
    'CSC(progenitor_cycling)':   '#FF9896',
    'CSC(progenitor_quiet)':     '#FF7F0E',
    'Hepatocytes':               '#1F77B4',
    'Hepatocytes (Metabolic-H)': '#AEC7E8',
    'Hepatocytes (Metabolic-L)': '#98DF8A',
    'Hepatocytes (cycling)':     '#2CA02C',
    'Hepatocytes (Stem)':        '#17BECF',
    'Hepatocytes (Quiescent)':   '#BCBD22',
    'Cholangiocytes':            '#9467BD',
    'Biliary-like':              '#C5B0D5',
    'Bipotent Progenitor':       '#8C564B',
}
DEFAULT_COLOR = '#D3D3D3'

def get_color(cell_type):
    return GLOBAL_PALETTE.get(str(cell_type).strip(), DEFAULT_COLOR)

# 1. Barcode
def find_best_match(adata_ids, file_ids):
    funcs = {
        'exact': lambda x: str(x),
        'pure_alnum': lambda x: re.sub(r'[^a-zA-Z0-9]', '', str(x)),
        'core_seq': lambda x: re.search(r'([ACGT]{12,20})', str(x)).group(1) if re.search(r'([ACGT]{12,20})', str(x)) else str(x),
        'no_suffix_alnum': lambda x: re.sub(r'[^a-zA-Z0-9]', '', re.sub(r'[-_.][0-9]+$', '', str(x))),
        'dot_to_dash': lambda x: str(x).replace('.', '-'),
        'dash_to_dot': lambda x: str(x).replace('-', '.'),
    }
    best_score, best_mapping, best_strategy = 0, {}, "None"
    for name, func in funcs.items():
        a_mapped = {x: func(x) for x in adata_ids}
        f_mapped = {func(x): x for x in file_ids}
        common = set(a_mapped.values()) & set(f_mapped.keys())
        if len(common) > best_score:
            best_score = len(common)
            best_strategy = name
            best_mapping = {a: f_mapped[a_mapped[a]] for a in adata_ids if a_mapped[a] in common}
    return best_score, best_strategy, best_mapping

def find_patient_file(pid, copykat_dirs, patient_id_mapping=None):
    search_pid = str(patient_id_mapping.get(pid, pid)) if patient_id_mapping else str(pid)
    safe_pid = re.escape(search_pid)
    gse_match = re.match(r'(GSE\d+)', search_pid)
    gse_prefix = gse_match.group(1) if gse_match else None
    core_id = search_pid.split('_')[-1]
    safe_core = re.escape(core_id)
    for d in copykat_dirs:
        if not os.path.exists(d): continue
        files = [f for f in os.listdir(d) if f.endswith("_CNA_results.txt") and "raw" not in f]
        for f in files:
            if re.search(rf"(^|[_.-]){safe_pid}([_.-]|$)", f):
                return os.path.join(d, f)
        if len(core_id) >= 4:
            for f in files:
                if gse_prefix and gse_prefix not in f: continue
                if re.search(rf"(^|[_.-]){safe_core}([_.-]|$)", f):
                    return os.path.join(d, f)
    return None

def extract_dataset(adata_path, cell_type_col, patient_col, copykat_dirs,
                    type_mapping=None, exclude_types=None, patient_id_mapping=None):
    dataset_name = os.path.basename(adata_path)
    print(f"\n{'='*60}")
    adata = sc.read_h5ad(adata_path)
    if patient_col not in adata.obs.columns:
        for col in ['sample', 'Sample', 'patient', 'patients', 'PatientID']:
            if col in adata.obs.columns: patient_col = col; break
    if patient_id_mapping: print(f"   Patient ID mapping: {patient_id_mapping}")

    adata.obs['cnv_score'] = np.nan
    unique_patients = adata.obs[patient_col].dropna().unique()
    total_filled = 0
    for pid in tqdm(unique_patients, desc="matching"):
        patient_cells = adata.obs[adata.obs[patient_col] == pid].index.tolist()
        if not patient_cells: continue
        target_file = find_patient_file(pid, copykat_dirs, patient_id_mapping)
        if not target_file:
            continue
        try:
            df_ck = pd.read_csv(target_file, sep="\t", engine='c')
            num_df = df_ck.select_dtypes(include=[np.number])
            cols = [c for c in num_df.columns if c not in ['chrom', 'chrompos', 'abspos']]
            if not cols: continue
            scores = (num_df[cols] ** 2).mean(axis=0)
            file_ids = scores.index.tolist()
            match_count, strategy, id_mapping = find_best_match(patient_cells, file_ids)
            if match_count > 0:
                score_series = pd.Series(
                    {cid: scores[id_mapping[cid]] for cid in patient_cells if cid in id_mapping}, dtype=float)
                adata.obs.loc[score_series.index, 'cnv_score'] = score_series
                total_filled += len(score_series.dropna())
            else:
                pass
        except Exception as e:
            pass

    subset = adata.obs[[patient_col, cell_type_col, 'cnv_score']].dropna().copy()
    subset.index.name = 'cell_id'
    subset.reset_index(inplace=True)
    subset.rename(columns={patient_col: 'patient', cell_type_col: 'cell_type'}, inplace=True)
    subset['dataset_source'] = dataset_name
    subset['cell_type'] = subset['cell_type'].astype(str)
    if type_mapping: subset['cell_type'] = subset['cell_type'].replace(type_mapping)
    if exclude_types: subset = subset[~subset['cell_type'].isin(exclude_types)]
    return subset

def get_stars(pval):
    if pval < 0.001: return "***"
    elif pval < 0.01: return "**"
    elif pval < 0.05: return "*"
    else: return "ns"

def plot_cnv_pdf_style(df_wide, type_a, type_b, ax):
    pair_data = df_wide[[type_a, type_b]].dropna()
    n_patients = len(pair_data)
    if n_patients < 2: return False

    try:
        stat, pval = stats.wilcoxon(pair_data[type_a], pair_data[type_b], alternative='two-sided')
    except ValueError:
        pval = 1.0

    plot_data = pair_data.reset_index().melt(
        id_vars=pair_data.index.name, value_vars=[type_a, type_b],
        var_name='Cell_Type', value_name='CNV_Score')
    palette = {type_a: get_color(type_a), type_b: get_color(type_b)}

    sns.boxplot(data=plot_data, x='Cell_Type', y='CNV_Score', order=[type_a, type_b],
                hue='Cell_Type', palette=palette, legend=False,
                showfliers=False, width=0.45, ax=ax, zorder=1,
                boxprops=dict(alpha=0.6, edgecolor='black', linewidth=0.8),
                medianprops=dict(color='black', linewidth=1.8),
                whiskerprops=dict(color='black', linewidth=0.8),
                capprops=dict(color='black', linewidth=0.8))

    for idx in pair_data.index:
        y1, y2 = pair_data.loc[idx, type_a], pair_data.loc[idx, type_b]
        ax.plot([0, 1], [y1, y2], color='gray', alpha=0.35, linewidth=0.7, zorder=2)

    sns.stripplot(data=plot_data, x='Cell_Type', y='CNV_Score', order=[type_a, type_b],
                  hue='Cell_Type', palette=palette, legend=False,
                  size=5, alpha=0.85, jitter=False, ax=ax, zorder=3,
                  edgecolor='black', linewidth=0.4)

    y_max = plot_data['CNV_Score'].max()
    y_min = plot_data['CNV_Score'].min()
    y_range = max(y_max - y_min, 0.01)

    from matplotlib.ticker import MultipleLocator
    ax.set_ylim(y_min - y_range * 0.08, y_max + y_range * 0.45)
    
    if y_range < 0.05:
        major_step, minor_step = 0.01, 0.005
    elif y_range < 0.2:
        major_step, minor_step = 0.05, 0.01
    elif y_range < 1.0:
        major_step, minor_step = 0.2, 0.05
    else:
        major_step, minor_step = 1.0, 0.2
    
    ax.yaxis.set_major_locator(MultipleLocator(major_step))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_step))
    ax.tick_params(axis='y', which='major', labelsize=9, length=4)
    ax.tick_params(axis='y', which='minor', length=2)

    line_y = y_max + y_range * 0.1
    ax.plot([0, 0, 1, 1], [line_y, line_y + y_range * 0.03, line_y + y_range * 0.03, line_y], lw=1, c='black')

    stars = get_stars(pval)
    if pval < 0.0001:
        p_text = f"{stars}\n$\\it{{P}}$ < 0.0001"
    else:
        p_text = f"{stars}\n$\\it{{P}}$ = {pval:.4f}"
    ax.text(0.5, line_y + y_range * 0.05, p_text, ha='center', va='bottom', fontsize=10, color='black')

    # ----  &  ----
    ax.set_title(f"{type_a} vs {type_b}\n$n={n_patients}$ patients", fontsize=11, pad=12)
    if ax.get_subplotspec().is_first_col():
        ax.set_ylabel("Mean CNV Score per patient", fontsize=11)
    else:
        ax.set_ylabel("")
    ax.set_xlabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.8, 1.8)
    ax.set_xticklabels([type_a, type_b], fontsize=9, rotation=30, ha='right')
    return True

dirs_data1 = ["/data1/copykat3", "/home/download/pysc/input/copykat"]
dirs_data2 = ["/data2/MACS/downstream/copykat_0308"]

df_1 = extract_dataset(
    adata_path="/data2/adata/adata_subset_0214.h5ad",
    cell_type_col="final_type",
    patient_col="patients",
    copykat_dirs=dirs_data1
)

_tmp = sc.read_h5ad("/data2/MACS/downstream/adata_hep_0308.h5ad")
_actual_pids = _tmp.obs['patient'].dropna().unique()
patient_id_mapping_2 = {}
for p in _actual_pids:
    patient_id_mapping_2[p] = f"Patient{p}"
    patient_id_mapping_2[str(p)] = f"Patient{p}"
del _tmp

df_2 = extract_dataset(
    adata_path="/data2/MACS/downstream/adata_hep_0308.h5ad",
    cell_type_col="predicted_cell_type",
    patient_col="patient",
    copykat_dirs=dirs_data2,
    type_mapping={
        "Cholangiocyte": "Cholangiocytes",
        "Hepatocyte": "Hepatocytes",
        "Hepatocytes (Metabolic Low)": "Hepatocytes (Metabolic-L)",
        "Hepatocytes (Metabolic High)": "Hepatocytes (Metabolic-H)",
        "Hepatocytes (Cycling)": "Hepatocytes (cycling)",
        "CSC(progenitor_prolif)": "CSC(progenitor_cycling)",
    },
    exclude_types=["Progenitor-like"],
    patient_id_mapping=patient_id_mapping_2
)

valid_dfs = [df for df in [df_1, df_2] if not df.empty]
if not valid_dfs: raise ValueError("No valid data extracted, aborting.")
df_combined = pd.concat(valid_dfs, ignore_index=True)

output_dir = "/home/download/csc_article/fig3/cnv"
os.makedirs(output_dir, exist_ok=True)

csv_path = os.path.join(output_dir, "CNV_Scores_All_Cells_Export.csv")
df_combined.to_csv(csv_path, index=False)

df_agg = df_combined.groupby(['patient', 'cell_type'], observed=False)['cnv_score'].agg(['mean', 'count']).reset_index()
df_agg.loc[df_agg['count'] < 2, 'mean'] = np.nan
df_wide = df_agg.pivot(index='patient', columns='cell_type', values='mean')
df_wide = df_wide.dropna(axis=1, how='all').dropna(axis=0, how='all')
type_order = df_wide.median().sort_values(ascending=False).index.tolist()

for ct in type_order:
    n = df_wide[ct].notna().sum()

key_types = [t for t in type_order if any(k in t for k in ['Bipotent', 'Cholan', 'Hepatocytes', 'CSC'])]
for t1, t2 in itertools.combinations(key_types, 2):
    n_overlap = len(df_wide[[t1, t2]].dropna())
    if n_overlap > 0:
        pass

# ----  PDF ----
output_pdf_path = os.path.join(output_dir, "CNV_Combined_Pairwise_Comparisons.pdf")

all_pairs = list(itertools.combinations(type_order, 2))

ordered_pairs = []
for t1, t2 in all_pairs:
    t1_is_target = 'csc' in t1.lower() or 'tumor' in t1.lower()
    t2_is_target = 'csc' in t2.lower() or 'tumor' in t2.lower()
    if t1_is_target and not t2_is_target: ordered_pairs.append((t2, t1))
    else: ordered_pairs.append((t1, t2))

drawable_pairs = []
skipped_pairs = []
for t1, t2 in ordered_pairs:
    n_overlap = len(df_wide[[t1, t2]].dropna())
    if n_overlap >= 2:
        drawable_pairs.append((t1, t2))
    else:
        skipped_pairs.append((t1, t2, n_overlap))

if skipped_pairs:
    for t1, t2, n in skipped_pairs:
        pass

rows, cols = 3, 3
plots_per_page = rows * cols

with PdfPages(output_pdf_path) as pdf:
    plot_idx = 0
    fig = None
    current_page_count = 0
    
    for t1, t2 in drawable_pairs:
        if current_page_count == plots_per_page or fig is None:
            if fig is not None:
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
            fig, axes = plt.subplots(rows, cols, figsize=(15, 16))
            axes = axes.flatten()
            current_page_count = 0
        
        success = plot_cnv_pdf_style(df_wide, t1, t2, axes[current_page_count])
        if success:
            current_page_count += 1
            plot_idx += 1
        else:
            pass
    
    if fig is not None and current_page_count > 0:
        for i in range(current_page_count, plots_per_page):
            fig.delaxes(axes[i])
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

print(f"    CSV: {csv_path}")
print(f"    PDF: {output_pdf_path}")
