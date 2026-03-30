#!/usr/bin/env python
# Figure S2 — BP identification, RNA velocity and pseudotime

# --- S2a: Bipotent progenitor identification (multi-parameter gating) ---
import os
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

adata_normal = sc.read_h5ad('/data2/adata/normal_annoted.h5ad')

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)

manual_thresholds = {
    'stem_min': 0.30, 
    'hep_max':  0.30, 
    'chol_min': 0.40, 
    'trop2_min': 0.20,
    'trop2_max': 3.00 
}

required_cols = ['Score_Fig3e_Hepa', 'Score_Fig3e_Prog', 'Score_Fig3e_Chol', 'Expr_TACSTD2']
df_scores = adata_normal.obs[required_cols].copy()

mask_final = (
    (df_scores['Score_Fig3e_Prog'] > manual_thresholds['stem_min']) & 
    (df_scores['Score_Fig3e_Hepa'] < manual_thresholds['hep_max']) & 
    (df_scores['Score_Fig3e_Chol'] > manual_thresholds['chol_min']) &
    (df_scores['Expr_TACSTD2']     < manual_thresholds['trop2_max'])
)

fig1, axs = plt.subplots(4, 1, figsize=(8, 12))

# Hep, Chol, Prog
configs = [
    ('Score_Fig3e_Hepa', '#E64B35', 'hep_max', 'Max <'),
    ('Score_Fig3e_Chol', '#4DBBD5', 'chol_min', 'Min >'),
    ('Score_Fig3e_Prog', '#00A087', 'stem_min', 'Min >')
]

for i, (col, color, thresh_key, label_pre) in enumerate(configs):
    sns.kdeplot(df_scores[col], fill=True, color=color, ax=axs[i])
    val = manual_thresholds[thresh_key]
    axs[i].axvline(val, color='red', linestyle='--', linewidth=1.5, label=f"{label_pre} {val}")
    axs[i].set_title(f'{i+1}. {col.split("_")[-1]} Score')
    axs[i].legend(loc='upper right', frameon=False)

# TROP2
sns.kdeplot(df_scores[df_scores['Expr_TACSTD2']>0]['Expr_TACSTD2'], fill=True, color='#F39C12', ax=axs[3])
axs[3].axvline(manual_thresholds['trop2_max'], color='red', linestyle='--', linewidth=1.5, label=f"Max < {manual_thresholds['trop2_max']}")
axs[3].set_title('4. TACSTD2 Expression')
axs[3].legend(loc='upper right', frameon=False)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Gating_KDE_HighRes.pdf'), dpi=300, bbox_inches='tight')
plt.gcf().set_dpi(72)
#plt.show()

fig2, axes = plt.subplots(1, 2, figsize=(14, 6))

# View 1: Chol vs Prog
sc1 = axes[0].scatter(df_scores['Score_Fig3e_Chol'], df_scores['Score_Fig3e_Prog'], 
                      c=df_scores['Score_Fig3e_Hepa'], cmap='viridis_r', s=3, alpha=0.5)
plt.colorbar(sc1, ax=axes[0], label='Hepa Score')
axes[0].axvline(manual_thresholds['chol_min'], color='black', linestyle=':', lw=1)
axes[0].axhline(manual_thresholds['stem_min'], color='black', linestyle=':', lw=1)
axes[0].set_title('Lineage View')

# View 2: TROP2 vs Prog
sc2 = axes[1].scatter(df_scores['Expr_TACSTD2'], df_scores['Score_Fig3e_Prog'], 
                      c=df_scores['Score_Fig3e_Chol'], cmap='plasma', s=3, alpha=0.5)
plt.colorbar(sc2, ax=axes[1], label='Chol Score')
rect = patches.Rectangle((0, manual_thresholds['stem_min']), manual_thresholds['trop2_max'], 1.5, 
                         linewidth=1.5, edgecolor='red', facecolor='none', label='Target Zone')
axes[1].add_patch(rect)
axes[1].set_title('TROP2 Gating View')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'Gating_Scatter_HighRes.pdf'), dpi=300, bbox_inches='tight')
plt.gcf().set_dpi(72)
#plt.show()

if sum(mask_final) > 0:
    adata_normal.obs['refined_type'] = adata_normal.obs['cell_type_mp'].astype(str)
    adata_normal.obs.loc[mask_final, 'refined_type'] = "Bipotent Progenitor"
    
    palette = {cat: "#e0e0e0" for cat in adata_normal.obs['refined_type'].unique()}
    palette["Bipotent Progenitor"] = "#E64B35"

    sc.pl.umap(adata_normal, color='refined_type', palette=palette, frameon=False, 
               title="Final Bipotent Progenitor Localization", show=False)
    
    plt.savefig(os.path.join(output_dir, 'Final_Bipotent_UMAP.pdf'), dpi=300, bbox_inches='tight')
    plt.gcf().set_dpi(80)
    #plt.show()

# --- S2c: Atlas parenchymal cell RNA velocity ---
# min_shared_counts:  spliced/unspliced
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# scVelo  PCA  Neighbors
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# 2.  (Streamplot)
# basis='umap'  adata X_umap
import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc

custom_colors = {
    'Biliary-like': '#F29E9A',
    'Bipotent Progenitor': '#DFB741',
    'CSC': '#DE8A3A',
    'CSC(progenitor_cycling)': '#D35E5B',  
    'CSC(progenitor_quiet)': '#C9AACB',    
    'Cholangiocytes': '#A1779C',           
    'Hepatocytes': '#5F80A5',
    'Hepatocytes (Metabolic-H)': '#5B9B53', 
    'Hepatocytes (Metabolic-L)': '#97CB7E', 
    'Hepatocytes (Quiescent)': '#A8CAE3',
    'Hepatocytes (Stem)': '#579B9A',
    'Hepatocytes (cycling)': '#94BABB'      
}

scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    color='final_type',
    palette=custom_colors,
    legend_loc='right margin',
    title='RNA Velocity (Stream)',
    show=False
)

# 3.  PDF
output_path = '/home/download/csc_article/fig2/all_result_corrected_colors.svg'
plt.savefig(output_path, bbox_inches='tight', dpi=200)

# plt.show()

# --- S2g,h: CytoTRACE2 pseudotime UMAP ---
import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)

# 2.  &  CytoTRACE2
def load_cytotrace2(base_path, adata):
    df_list = []
    for gse_id in os.listdir(base_path):
        file_path = os.path.join(base_path, gse_id, 'cytotrace2_results', 'cytotrace2_results.txt')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep=r'\s+', index_col=0)[
                ['CytoTRACE2_Score', 'CytoTRACE2_Potency', 'CytoTRACE2_Relative']]
            df_list.append(df)
    if not df_list:
        return
    master_df = pd.concat(df_list)
    adata.obs[['CytoTRACE2_Score', 'CytoTRACE2_Potency', 'CytoTRACE2_Relative']] = master_df.reindex(adata.obs_names)
    scored = adata.obs['CytoTRACE2_Score'].notna().sum()

load_cytotrace2('/home/download/pysc/input/cytotrace2_input_tumor', adata_tumor)
load_cytotrace2('/home/download/pysc/input/cytotrace2_input_normal', adata_normal)

cytotrace_cmap = mcolors.LinearSegmentedColormap.from_list(
    'cytotrace2', ['#440154', '#31688E', '#35B779', '#FDE725'])

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata_normal,
           color='CytoTRACE2_Score',
           cmap=cytotrace_cmap,
           vmin=0, vmax=0.2,
           frameon=False,
           title='CytoTRACE2 Score (Normal)',
           colorbar_loc='right',
           ax=axes[0],
           show=False)

sc.pl.umap(adata_tumor,
           color='CytoTRACE2_Score',
           cmap=cytotrace_cmap,
           vmin=0, vmax=0.6,
           frameon=False,
           title='CytoTRACE2 Score (Tumor)',
           colorbar_loc='right',
           ax=axes[1],
           show=False)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'CytoTRACE2_normal_tumor.svg'),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# --- Malignant status UMAP ---
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

sc.pl.umap(
    adata_subset,
    color='malignant',
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
output_path = os.path.join(output_dir, 'UMAP_malignant.pdf')

plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

plt.show() 

plt.close()

# --- S2b: CSC subpopulation RNA velocity (scVI/scANVI + scVelo) ---

target_types = ['CSC(progenitor_cycling)', 'CSC(progenitor_quiet)', 'CSC']

# 2.  adata_subset  adata_csc
adata_csc = adata[adata.obs['final_type'].isin(target_types)].copy()

import scanpy as sc

cell_list = [
    "CSC(progenitor_cycling)",
    "CSC(progenitor_quiet)",
    "CSC"
]

adata_filtered = adata[adata.obs["merged_cell_type"].isin(cell_list)].copy()

print("-" * 30)
print(adata_filtered.obs["merged_cell_type"].value_counts())

import scanpy as sc
import scvi
import numpy as np

if "raw_counts" in adata_filtered.layers:
    adata_filtered.X = adata_filtered.layers["raw_counts"].copy()
elif adata_filtered.raw is not None:
    adata_filtered.X = adata_filtered.raw.X.toarray() if hasattr(adata_filtered.raw.X, "toarray") else adata_filtered.raw.X
else:
    pass

adata_filtered.layers["counts"] = adata_filtered.X.copy()

sc.pp.normalize_total(adata_filtered, target_sum=1e4)
sc.pp.log1p(adata_filtered)

sc.pp.highly_variable_genes(
    adata_filtered,
    n_top_genes=2000,
    subset=True,
    flavor="seurat",
)

# 3.  scVI  (Setup AnnData)
scvi.model.SCVI.setup_anndata(
    adata_filtered, 
    layer="counts",
    batch_key="batch",
    labels_key="merged_cell_type"
)

# 3.  scVI  (Setup AnnData)
# batch_key="batchs"
scvi.model.SCVI.setup_anndata(
    adata_filtered, 
    layer="counts", 
    batch_key="batch", 
    labels_key="merged_cell_type" 
)

model_scvi = scvi.model.SCVI(adata_filtered, n_layers=2, n_latent=30, gene_likelihood="nb")
model_scvi.train(max_epochs=200, early_stopping=True)

model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi, 
    unlabeled_category="Unknown"
)

# 6.  SCANVI
model_scanvi.train(max_epochs=100, early_stopping=True)

# 7.  (Latent Representation)
adata_filtered.obsm["X_scANVI"] = model_scanvi.get_latent_representation(adata_filtered)

sc.pp.neighbors(adata_filtered, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata_filtered)

# adata_filtered.write("scanvi_integrated.h5ad", compression="gzip")

sc.pl.umap(
    adata_filtered, 
    color=["merged_cell_type", "batch"], 
    title=["Cell Type (SCANVI integrated)", "Batch Effect Check"],
    wspace=0.3,
    frameon=False
)

# 1.  Moments ( neighbors )
scv.pp.moments(adata_filtered, n_pcs=None, n_neighbors=None)

scv.tl.velocity(adata_filtered, mode='stochastic')
scv.tl.velocity_graph(adata_filtered)

scv.pl.velocity_embedding_stream(
    adata_filtered, 
    basis='umap', 
    color='merged_cell_type', 
    title='RNA Velocity in MP-space (Normal Cells)',
    legend_loc='right margin'
)

# --- S2e,f: scTour pseudotime (tumor-adjacent & tumor) ---

import sctour as sct
import scanpy as sc
import numpy as np
import scipy.sparse as sparse
adata=sc.read("/home/download/csc_article/monocle4bile/ref_out/MACS_import/bile.h5ad")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# (Highly Variable Genes, HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
adata.X = adata.layers['raw_counts'].copy()
if sparse.issparse(adata.X):
    adata.X = adata.X.toarray()

tnr = sct.train.Trainer(
    adata, 
    loss_mode='nb', 
    n_latent=20, 
    n_ode_hidden=25, 
    n_vae_hidden=25,
    nepoch=400, 
    lr=0.001
)

tnr.train()

import os
if not os.path.exists('model_dir'):
    os.makedirs('model_dir')
tnr.save_model(save_dir='model_dir', save_prefix='my_sctour_model')

pseudotime = tnr.get_time()
adata.obs['sctour_pseudotime'] = pseudotime

# get_latentsp  tuple: (mix_zs, zs, pred_zs)
mix_zs, _, _ = tnr.get_latentsp() 
adata.obsm['X_sctour'] = mix_zs

vf = tnr.get_vector_field(pseudotime, mix_zs)
adata.obsm['X_sctour_vf'] = vf

import pandas as pd
import matplotlib.pyplot as plt
import os

scanvi_umap_path = "/home/download/csc_article/monocle4bile/ref_out/MACS_import/umap_coords.csv"

if os.path.exists(scanvi_umap_path):
    umap_df = pd.read_csv(scanvi_umap_path, index_col=0)
    
    adata.obsm['X_umap'] = umap_df.loc[adata.obs_names, :].values
else:
    raise FileNotFoundError(f" scANVI UMAP file not found: {scanvi_umap_path}")

# 7.  ( scANVI UMAP)

adata = adata[np.argsort(adata.obs['sctour_pseudotime'].values), :]

color_key = "final_annotation" if "final_annotation" in adata.obs else "leiden"

fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(14, 6))

# --- : scANVI UMAP +  ---
sc.pl.umap(
    adata, 
    color=color_key, 
    ax=axs[0], 
    show=False, 
    frameon=False, 
    title="scANVI UMAP (Cell Types)",
    legend_loc="on data",
    legend_fontsize=8
)

# --- : scANVI UMAP + sctour  ---

sct.vf.plot_vector_field(
    adata, 
    zs_key='X_sctour',
    vf_key='X_sctour_vf',    # sctour
    use_rep_neigh='X_sctour',
    E_key='X_umap',
    color=color_key, 
    ax=axs[1], 
    show=False, 
    frameon=False, 
    title="scANVI UMAP + scTour Trajectory",
    
    size=100,
    alpha=0.1,
    stream_density=1.5,
    stream_linewidth=1.2,
    stream_arrowsize=1.2,
    stream_color='k'
)

plt.tight_layout()
plt.show()

# plt.savefig("sctour_trajectory_on_scanvi_umap.pdf")

import sctour as sct

fig, ax = plt.subplots(figsize=(7, 6))

sct.vf.plot_vector_field(
    adata, 
    zs_key='X_sctour',
    vf_key='X_sctour_vf',
    use_rep_neigh='X_sctour', 
    E_key='X_umap',          # scANVI/scVI  UMAP
    color='sctour_pseudotime',
    ax=ax,
    show=False, 
    frameon=False, 
    title="Pseudotime + Vector Field",
    
    cmap='viridis',
    size=100, alpha=0.8,
    stream_color='white',
    stream_linewidth=1.0,
    stream_density=1.5
)

plt.show()

sc.pl.violin(
    adata, 
    keys='sctour_pseudotime', 
    groupby='final_annotation', 
    rotation=90,
    order=None,              # ['Stem', 'Progenitor', 'Mature']
    ylabel='Inferred Pseudotime'
)

import sctour as sct
import scanpy as sc
import numpy as np
import scipy.sparse as sparse
adata=sc.read("/home/download/csc_article/monocle4bile/adata_normal_MACS.h5ad")

sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_genes(adata, min_cells=20)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# (Highly Variable Genes, HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
adata.X = adata.layers['raw_counts'].copy()
if sparse.issparse(adata.X):
    adata.X = adata.X.toarray()

tnr = sct.train.Trainer(
    adata, 
    loss_mode='nb', 
    
    n_latent=15,       
    
    alpha_recon_lec=0.3,
    alpha_recon_lode=0.7,
    
    n_ode_hidden=25, 
    n_vae_hidden=25,
    
    nepoch=1200,       
    lr=0.001
)

tnr.train()

adata.obs['sctour_pseudotime'] = tnr.get_time()
mix_zs, _, _ = tnr.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
adata.obsm['X_sctour'] = mix_zs
adata.obsm['X_sctour_vf'] = tnr.get_vector_field(adata.obs['sctour_pseudotime'].values, mix_zs)

import sctour as sct
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# 2.  (Vector Field on Original UMAP)

if 'X_sctour' not in adata.obsm:
    mix_zs, _, _ = tnr.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
    adata.obsm['X_sctour'] = mix_zs
    
if 'X_sctour_vf' not in adata.obsm:
    adata.obsm['X_sctour_vf'] = tnr.get_vector_field(adata.obs['sctour_pseudotime'].values, adata.obsm['X_sctour'])

fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(16, 7))

color_key = "final_annotation" if "final_annotation" in adata.obs else "leiden"

# ---  UMAP () ---
sc.pl.umap(
    adata, 
    color=color_key, 
    ax=axs[0], 
    show=False, 
    title='Original UMAP (Cell Types)',
    legend_loc='on data',
    legend_fontsize=8,
    frameon=False
)

# ---  UMAP + sctour  ---
sct.vf.plot_vector_field(
    adata, 
    zs_key='X_sctour',
    vf_key='X_sctour_vf',
    use_rep_neigh='X_sctour',
    
    E_key='X_umap',
    
    color=color_key, 
    ax=axs[1], 
    show=False, 
    title='Trajectory on Original UMAP',
    frameon=False,
    
    stream_density=1.5,
    stream_linewidth=1.0,
    stream_arrowsize=1.2,
    stream_color='k',
    size=100,
    alpha=0.1
)

plt.tight_layout()
plt.show()

# plt.savefig("sctour_trajectory_on_original_umap.pdf")

import numpy as np

# 1.  X_umap  NaN
nan_mask = np.isnan(adata.obsm['X_umap']).any(axis=1)
n_nans = nan_mask.sum()

if n_nans > 0:
    
    adata_clean = adata[~nan_mask].copy()
    
else:
    adata_clean = adata.copy()

if 'X_sctour_vf' in adata_clean.obsm:
    vf_nan_mask = np.isnan(adata_clean.obsm['X_sctour_vf']).any(axis=1)
    n_vf_nans = vf_nan_mask.sum()
    if n_vf_nans > 0:
        adata_clean = adata_clean[~vf_nan_mask].copy()

import matplotlib.pyplot as plt
import sctour as sct

fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(15, 6))

sc.pl.umap(
    adata_clean, 
    color=color_key, 
    ax=axs[0], 
    show=False, 
    frameon=False, 
    title="Cleaned UMAP"
)

sct.vf.plot_vector_field(
    adata_clean,
    zs_key='X_sctour', 
    vf_key='X_sctour_vf', 
    use_rep_neigh='X_sctour',
    E_key='X_umap',
    color=color_key, 
    ax=axs[1], 
    show=False, 
    title='Trajectory on Original UMAP',
    frameon=False,
    
    stream_density=1.5,
    stream_linewidth=1.0,
    stream_arrowsize=1.2,
    stream_color='k',
    size=100,
    alpha=0.1
)

plt.tight_layout()
plt.show()

import scanpy as sc
import sctour as sct
import pandas as pd
import numpy as np
import torch
from scipy import sparse
from scipy.sparse import csr_matrix

file_path = '/home/download/csc_article/monocle4bile/adata_normal_MACS.h5ad'
adata_sctour = sc.read_h5ad(file_path)

adata_sctour.X = adata_sctour.layers['raw_counts'].copy()

if sparse.issparse(adata_sctour.X):
    adata_sctour.X.data = np.round(adata_sctour.X.data)
else:
    adata_sctour.X = np.round(adata_sctour.X)

min_val = adata_sctour.X.min()

# 3.  HVG ( Seurat )

adata_temp = adata_sctour.copy()
sc.pp.normalize_total(adata_temp, target_sum=1e4)
sc.pp.log1p(adata_temp)

# 'seurat' flavor
sc.pp.highly_variable_genes(
    adata_temp, 
    n_top_genes=2000, 
    flavor='seurat', 
    subset=False
)

adata_sctour.var['highly_variable'] = adata_temp.var['highly_variable']
del adata_temp

adata_sctour = adata_sctour[:, adata_sctour.var['highly_variable']].copy()

if sparse.issparse(adata_sctour.X):
    adata_sctour.X = adata_sctour.X.toarray()

adata_sctour.X = adata_sctour.X.astype(np.float32)
adata_sctour.X[adata_sctour.X < 0] = 0

# print  <class 'numpy.ndarray'>

# 5.  sctour
device = 'cuda' if torch.cuda.is_available() else 'cpu'

tnode = sct.train.Trainer(
    adata_sctour, 
    loss_mode='nb',
    alpha_recon_lec=0.5, 
    alpha_recon_lode=0.5,
    use_gpu=(device == 'cuda')
)

tnode.train()

pseudotime = tnode.get_time()
adata_sctour.obs['sctour_pseudotime'] = pseudotime

csv_path = '/data2/adata/liver_pseudotime_scores_MACS.csv'
adata_sctour.obs[['sctour_pseudotime']].to_csv(csv_path)

# 7.  adata_normal
if 'adata_normal' in locals():
    ptime_series = adata_sctour.obs['sctour_pseudotime']
    adata_normal.obs['sctour_pseudotime'] = adata_normal.obs_names.map(ptime_series)
    
    na_count = adata_normal.obs['sctour_pseudotime'].isna().sum()
    if na_count > 0:
        pass
        
else:
    pass

import sctour as sct
import torch
import os
import pandas as pd

tnr = sct.train.Trainer(adata, loss_mode='nb', use_gpu=torch.cuda.is_available())
tnr.load_model(save_dir='model_dir', save_prefix='my_sctour_model')

adata.obs['sctour_pseudotime'] = tnr.get_time()

df_export = pd.DataFrame(adata.obs['sctour_pseudotime'])

output_path = '/data2/adata/liver_pseudotime_scores.csv'

# index=True  Barcode index_label
df_export.to_csv(output_path, index=True, index_label='barcode')

print(df_export.head())
