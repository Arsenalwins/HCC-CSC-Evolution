#!/usr/bin/env python
# Figure S1 — HCC single-cell atlas functional subgroup identification

# --- S1a: Atlas marker gene dotplot ---
import scanpy as sc
import matplotlib.pyplot as plt
import os

output_dir = '/home/download/csc_article/fig1'
os.makedirs(output_dir, exist_ok=True)

sc.settings.set_figure_params(dpi=300, format='pdf', transparent=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

markers_dict = {
    "Immune": ['PTPRC'],
    'B cells': ['MZB1'],
    'Endothelial': ['VWF', 'PECAM1', 'PLPP1'],
    'Fibroblast': ['ACTA2', 'MYL9', 'TAGLN'],
    'Hepatocytes': ['ALB', 'APOA2', 'APOC3'],
    'Mast cells': ['TPSB2', 'TPSAB1', 'HPGDS', 'CPA3'],
    'Myeloid': ['LYZ', 'CD74'],
    'NK cells': ['GNLY', 'NKG7'],
    'NK1': ['FGFBP2', 'MYOM2', 'GZMB', 'GZMH'],
    'NK2': ['KLRB1', 'KLRF1', 'CD7', 'NKG7', 'CMC1', 'XCL2'],
    'T cells': ['CD3D', 'CD3E'],
    'TAM': ['C1QA', 'C1QB', 'CD163'],
    'Bipotent': ['KRT19', 'TNFRSF12A', "SOX9", "HNF4A"],
    "Stemness": ['EPCAM', 'PROM1']
}

dp = sc.pl.dotplot(
    adata,
    var_names=markers_dict,
    groupby='leiden',
    standard_scale='var',       
    return_fig=True
)

dp.style(
    cmap='Reds',               
    dot_edge_color='black',    
    dot_edge_lw=0.5            
)

dp.legend(
    colorbar_title='Scaled\nExpression', 
    size_title='Fraction of\nCells'      
)

output_path = os.path.join(output_dir, 'fig1_marker_dotplot.pdf')

dp.savefig(output_path, dpi=300, bbox_inches='tight')

# --- S1c: MP functional module UMAP ---
adata_subset=sc.read("/data2/adata/adata_subset_0214.h5ad")
import scanpy as sc
import matplotlib.pyplot as plt
import os

plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

sc.settings.set_figure_params(
    dpi=300, 
    dpi_save=300, 
    format='pdf', 
    transparent=True, 
    figsize=(6, 6)
)

fig, ax = plt.subplots(figsize=(7, 6))

# basis='umap_mp'  basis='umap_scANVI_State'
sc.pl.umap(
    adata_subset,
    color='cell_type_mp_final',
    frameon=False,
    title='',
    legend_fontsize=10,
    legend_fontweight='normal',
    legend_loc='right margin',
    size=10,
    show=False,
    ax=ax
)

ax.set_xlabel("UMAP 1", fontsize=14, fontweight='bold')
ax.set_ylabel("UMAP 2", fontsize=14, fontweight='bold')
ax.set_xticks([])
ax.set_yticks([])

output_dir = '/home/download/csc_article/fig1/'
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'UMAP_cell_type_mp_final.pdf')

plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

plt.show() 

plt.close()

# --- S1d: MP score dotplot ---
import scanpy as sc
import matplotlib.pyplot as plt
import os

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

sc.settings.set_figure_params(dpi=300, format='pdf', transparent=True)

mp_mapping = {
    'Gavish_MP_1_norm': 'MP2',
    'Gavish_MP_5_norm': 'MP5',
    'Gavish_MP_2_norm': 'MP3',
    'Gavish_MP_10_norm': 'MP9',
    'Gavish_MP_3_norm': 'MP1',
    'Gavish_MP_7_norm': 'MP4',
    'Gavish_MP_11_norm': 'MP12',
    'Gavish_MP_12_norm': 'MP7',
    'Gavish_MP_4_norm': 'MP11',
    'Gavish_MP_9_norm': 'MP10',
    'Gavish_MP_8_norm': 'MP_8'
}

for old_col, new_col in mp_mapping.items():
    if old_col in adata_subset.obs.columns:
        adata_subset.obs[new_col] = adata_subset.obs[old_col]
    else:
        pass

var_groups = {
    'Cell Cycle': ['MP2', 'MP5'],
    'Biliary Lineage': ['MP3', 'MP9'],
    'Liver Function': ['MP1', 'MP4', 'MP12', 'MP7'],
    'Metabolic/Hypoxia': ['MP11', 'MP10'],
    'Immune/Stress': ['MP_8']
}

output_dir = '/home/download/csc_article/fig1/'
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'dotplot_mp_functions.pdf')

n_clusters = len(adata_subset.obs['leiden_mp_norm'].unique())
n_features = len(mp_mapping)
figsize = (n_clusters * 0.4 + 2, n_features * 0.35 + 2)

dp = sc.pl.dotplot(
    adata_subset,
    var_names=var_groups,         
    groupby='cell_type_mp_final', 
    standard_scale=None,          
    swap_axes=True,               
    figsize=figsize,
    show=False,
    return_fig=True
)

dp.style(cmap='Reds', dot_edge_color='black', dot_edge_lw=0.5)

axes_dict = dp.get_axes()
ax = axes_dict['mainplot_ax']

for label in ax.get_xticklabels():
    label.set_rotation(-90) 

fig = ax.figure
fig.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

plt.show()

# --- S1b: NMF clustering (Jaccard similarity heatmap) ---

import pandas as pd
import numpy as np
import pickle
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from tqdm import tqdm

STEP1_FILE = "/home/download/pysc/input/gavish_nmf_output/step1_raw_programs.pkl"
OUTPUT_DIR = "nmf_visualization"
MP_FILE="/home/download/pysc/input/gavish_nmf_output/final_gavish_mps.pkl"

FILTER_THRESHOLD = 0.1 

CUSTOM_COLORS = [
    "#2f2053", # Deep Purple (Low)
    "#712580", # Purple
    "#be4071", # Pink/Red
    "#e47c62", # Orange
    "#ffffff", # Light Orange
    "#ffffff"  # Pale Yellow (High)
]

CUSTOM_COLORS = CUSTOM_COLORS[::-1] 
MY_CMAP = LinearSegmentedColormap.from_list("MyCustomGradient", CUSTOM_COLORS)
VMAX = 0.5 

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def plot_sorted_heatmap():
    if not os.path.exists(STEP1_FILE) or not os.path.exists(MP_FILE):
        return

    with open(STEP1_FILE, 'rb') as f:
        step1_data = pickle.load(f)
    with open(MP_FILE, 'rb') as f:
        mp_data = pickle.load(f)

    mapping_list = []
    
    mp_names = list(mp_data.keys())
    mp_sets = {mp: set(genes) for mp, genes in mp_data.items()}
    
    for prog_id, genes in tqdm(step1_data.items(), desc="Mapping"):
        prog_set = set(genes)
        best_mp = None
        max_score = -1
        
        for mp in mp_names:
            score = len(prog_set & mp_sets[mp]) / len(prog_set | mp_sets[mp])
            if score > max_score:
                max_score = score
                best_mp = mp
        
        if max_score >= FILTER_THRESHOLD:
            mapping_list.append({
                "Program_ID": prog_id,
                "Assigned_MP": best_mp,
                "Genes": genes,
                "Score": max_score
            })
            
    df_map = pd.DataFrame(mapping_list)
    if df_map.empty: return

    mp_stats = df_map.groupby('Assigned_MP')['Score'].agg(['count', 'mean'])
    
    # ascending=[False, False]
    mp_stats = mp_stats.sort_values(by=['count', 'mean'], ascending=[False, False])
    
    sorted_mps = mp_stats.index.tolist()
    mp_rank_map = {mp: i for i, mp in enumerate(sorted_mps)}
    
    df_map['MP_Rank'] = df_map['Assigned_MP'].map(mp_rank_map)
    
    df_sorted = df_map.sort_values(by=['MP_Rank', 'Score'], ascending=[True, False])
    
    sorted_progs = df_sorted['Program_ID'].tolist()
    sorted_genes = df_sorted['Genes'].tolist()
    assigned_mps = df_sorted['Assigned_MP'].tolist()
    
    n = len(sorted_progs)
    mat = np.zeros((n, n))
    gene_sets = [set(g) for g in sorted_genes]
    
    for i in tqdm(range(n), desc="Matrix"):
        for j in range(i, n):
            s1, s2 = gene_sets[i], gene_sets[j]
            val = len(s1 & s2) / len(s1 | s2) if len(s1 | s2) > 0 else 0
            mat[i, j] = val
            mat[j, i] = val

    unique_mps = df_sorted['Assigned_MP'].unique()
    palette = sns.color_palette("husl", len(unique_mps))
    mp_color_map = dict(zip(unique_mps, palette))
    row_colors = pd.Series(assigned_mps).map(mp_color_map)

    plt.figure(figsize=(16, 15))
    
    g = sns.clustermap(
        pd.DataFrame(mat, index=sorted_progs, columns=sorted_progs),
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors.values,
        col_colors=row_colors.values,
        cmap=MY_CMAP,
        vmin=0, vmax=VMAX,
        xticklabels=False, yticklabels=False,
        figsize=(14, 14),
        dendrogram_ratio=0.01,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        cbar_kws={'label': 'Jaccard Similarity'}
    )
    
    ncol = 2 if len(unique_mps) > 25 else 1
    legend_handles = [Patch(facecolor=mp_color_map[mp], label=mp) for mp in unique_mps]
    plt.legend(handles=legend_handles, title="Consensus MP (Ordered by Size)", 
               bbox_to_anchor=(1.05, 1), loc='upper left', ncol=ncol, frameon=False)
    
    save_path = os.path.join(OUTPUT_DIR, "sorted_best_heatmap.pdf")
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    plot_sorted_heatmap()
