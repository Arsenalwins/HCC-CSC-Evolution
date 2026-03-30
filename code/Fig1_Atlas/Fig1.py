#!/usr/bin/env python
# Figure 1 — Single-cell atlas construction and CSC/BP identification

# --- 1b: scVI/scANVI integration + UMAP ---
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import os
adata = sc.read_h5ad('/data2/adata/annoted.h5ad')
adata_filtered = adata[adata.obs['leiden'] != 'unassign'].copy()
adata_filtered.obs['leiden'] = adata_filtered.obs['leiden'].cat.remove_unused_categories()
adata_filtered.X=adata_filtered.layers["raw_counts"].copy()
adata_filtered.layers["counts"]=adata_filtered.layers["raw_counts"].copy()
adata_hvg = adata_filtered[:, adata_filtered.var['highly_variable']].copy()

scvi.model.SCVI.setup_anndata(
    adata_hvg, 
    layer="counts", 
    batch_key="batch"  
)

vae = scvi.model.SCVI(adata_hvg)
vae.train() 

# 4.  scVI  scANVI ( HVG )
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata_hvg,               
    labels_key="leiden",           
    unlabeled_category="Unknown"   
)

scanvi_model.train(max_epochs=20, n_samples_per_label=100)

adata_filtered.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata_hvg)

sc.pp.neighbors(adata_filtered, use_rep="X_scANVI")
sc.tl.umap(adata_filtered, min_dist=0.3)

output_dir = '/data2/adata/'
os.makedirs(output_dir, exist_ok=True)

sc.settings.set_figure_params(dpi=300, format='pdf', transparent=True)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

fig, ax = plt.subplots(figsize=(6, 5))

sc.pl.umap(
    adata_filtered,
    color='leiden',
    frameon=False,
    title='scANVI Integrated UMAP',
    legend_loc='right margin',
    show=False,
    ax=ax
)

ax.set_xlabel("UMAP 1", fontsize=12, fontweight='bold')
ax.set_ylabel("UMAP 2", fontsize=12, fontweight='bold')
ax.set_xticks([])
ax.set_yticks([])

output_path = os.path.join(output_dir, 'umap_scanvi_integrated.pdf')
plt.savefig(output_path, dpi=300, bbox_inches='tight')

plt.show()

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.set_figure_params(dpi=300, dpi_save=300, format='pdf', figsize=(5, 5), transparent=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

adata = sc.read_h5ad('/data2/adata/annoted.h5ad')

cluster_col = 'leiden' 

adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')

n_clusters = len(adata.obs[cluster_col].cat.categories)

if n_clusters <= 20:
    palette = sns.color_palette("tab20", n_clusters)
else:
    palette = sns.color_palette("husl", n_clusters)

hex_colors = palette.as_hex()

# 4.  UMAP
fig, ax = plt.subplots(figsize=(6, 5))

sc.pl.umap(
    adata,
    color=cluster_col,
    palette=hex_colors,
    legend_loc='right margin',
    legend_fontsize=10,
    legend_fontweight='normal',
    frameon=False,
    title='UMAP of Cell Types',
    show=False,
    ax=ax
)

ax.set_xlabel("UMAP 1", fontsize=12, fontweight='bold')
ax.set_ylabel("UMAP 2", fontsize=12, fontweight='bold')

ax.set_xticks([])
ax.set_yticks([])

output_path = '/data2/adata/umap_publication.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight')

plt.show()

# --- 1c: CytoTRACE2 stemness score KDE ---
adata_subset=sc.read("/data2/adata/adata_subset_0214.h5ad")
adata_tumor=adata_subset[adata_subset.obs["malignant"]=="tumor"]
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_cytotrace2(base_path, adata):
    
    df_list = []
    for gse_id in os.listdir(base_path):
        file_path = os.path.join(base_path, gse_id, 'cytotrace2_results', 'cytotrace2_results.txt')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep=r'\s+', index_col=0)[['CytoTRACE2_Score', 'CytoTRACE2_Potency', 'CytoTRACE2_Relative']]
            df_list.append(df)
            
    if not df_list:
        return
        
    master_df = pd.concat(df_list)
    adata.obs[['CytoTRACE2_Score', 'CytoTRACE2_Potency', 'CytoTRACE2_Relative']] = master_df.reindex(adata.obs_names)
    
    scored_cells = adata.obs['CytoTRACE2_Score'].notna().sum()

path_tumor = '/home/download/pysc/input/cytotrace2_input_tumor'
load_cytotrace2(path_tumor, adata_tumor)

# 2.  adata_tumor  CytoTRACE2

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial']

sns.set_theme(style="ticks")
fig, ax = plt.subplots(figsize=(7, 5))

sns.kdeplot(
    data=adata_tumor.obs, 
    x='CytoTRACE2_Score',
    fill=True, 
    common_norm=False,
    palette='viridis', 
    alpha=0.5,
    linewidth=1.5,
    ax=ax
)

plt.title('Density Distribution of CytoTRACE2 Score in Tumor', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('CytoTRACE2 Score (Stem-like <-> Differentiated)', fontsize=12, fontweight='bold')
plt.ylabel('Density', fontsize=12, fontweight='bold')

sns.despine()

leg = ax.get_legend()
if leg is not None:
    leg.set_title('Potency')
    leg.set_bbox_to_anchor((1.05, 1))
    leg.set_frame_on(False)

output_path = '/home/download/csc_article/fig2/CytoTRACE2_density_tumor.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

plt.show()

# --- 1d: CSC subpopulation marker heatmap (AUCell) ---
target_types = ['CSC(progenitor_cycling)', 'CSC(progenitor_quiet)', 'CSC']

# 2.  adata_subset  adata_csc
adata_csc = adata_subset[adata_subset.obs['final_type'].isin(target_types)].copy()
import os
import scanpy as sc
import gseapy as gp
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

cluster_key = "final_type" 

if adata_csc.obs[cluster_key].dtype.name == 'category':
    adata_csc.obs[cluster_key] = adata_csc.obs[cluster_key].cat.remove_unused_categories()

custom_signatures = {
    'Hepatocytes': ['ALB', 'ASGR1', 'APOA2', 'APOC3', 'CYP3A4', 'GLUL'],
    'BEC': ['KRT7', 'CD24', 'MUC5B', 'MUC3A', 'SPP1'],
    'Bipotent': ['KRT19', 'TNFRSF12A', 'SOX9', 'HNF4A', 'SERPINE1', 'FKBP5', 'FN1', 
                 'GPC6', 'IL18', 'BICC1', 'CDH6', 'CREB5', 'HDAC9', 'FGFR2', 'SOX4', 'SOX6', 'PTCHD4'],
    'JHEP_CSC_Like': ['CD24', 'ICAM1', 'ACSL4', 'GOLGA8B', 'BAG3', 'RBM26', 'C17orf67'],
    'JHEP_CSC_Con': ['EPCAM', 'KRT19', 'CD44', 'CD24', 'PROM1', 'TACSTD2'],
    'Review_HCC_CSC': ['PROM1', 'CD44', 'EPCAM', 'CD24', 'KRT19', 'THY1', 'ANPEP'],
    'CURATED_Stemness_markers': ['ABCG2', 'BMI1', 'CD34', 'CD44', 'CTNNB1', 'EPAS1', 'EZH2', 
                                 'HIF1A', 'KDM5B', 'KLF4', 'LGR5', 'MYC', 'NANOG', 'NES', 
                                 'NOTCH1', 'POU5F1', 'PROM1', 'SOX2', 'TWIST1', 'ZFP42', 'ZSCAN4'],
    'NPJ_CSC': ['ENO1', 'HSPB1', 'GLO1', 'TKT', 'PON1', 'PTP4A1', 'KRT8', 'POLR2J', 'PERP', 'CD63', 
                'PRDX4', 'KDELR2', 'G6PD', 'GPC3', 'MALSU1', 'GSTM3', 'DEFB1', 'CYCS', 'S100A13', 
                'LAPTM4B', 'GGCT', 'CCT3', 'POLR2K', 'MPLKIP', 'PFDN2', 'ENY2', 'PPIA', 'ELF3', 
                'MRPL13', 'SEC61G', 'SNRPE', 'SQLE', 'CCT6A', 'COA6', 'TSTA3', 'CHCHD2', 'COX20', 
                'PUF60', 'ANXA4', 'STOML2', 'NDUFB3', 'CLTA', 'CMTM8', 'LDHA', 'GLB1', 'CCDC34', 
                'SPCS1', 'MDK', 'ALB', 'COA4', 'CCT5', 'ATP5C1', 'SEPP1', 'NUDT5', 'HMGCS1', 
                'PRDX3', 'VDAC1', 'BNIP3', 'NHP2', 'GAPDH', 'TUBB2A', 'TPI1', 'MRPS18B', 'TMEM106C', 
                'TUBB', 'CDK4', 'PFDN6', 'TSPAN8', 'CUTA', 'PSMB5', 'HMGA1', 'GPX2', 'SNRPC', 'CKB']
}

msigdb_targets = {
    'IL6_JAK_STAT3': 'IL-6/JAK/STAT3 Signaling', 'WNT_Signaling_Online': 'Wnt-beta Catenin Signaling',
    'NOTCH_Signaling': 'Notch Signaling', 'HEDGEHOG_Signaling': 'Hedgehog Signaling',
    'TGF_BETA_Signaling': 'TGF-beta Signaling', 'EMT_Signaling': 'Epithelial Mesenchymal Transition',
    'CellCycle_S_Phase': 'E2F Targets', 'CellCycle_G2M_Phase': 'G2-M Checkpoint',
    'Proliferation_MYC': 'Myc Targets V1', 'JNK_MAPK_Signaling': 'BIOCARTA_MAPK_PATHWAY',
    'ERK_Pathway': 'BIOCARTA_ERK_PATHWAY', 'MET_Pathway': 'PID_MET_PATHWAY',
    'TLR4_Signaling': 'REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE', 'BenPorath_Prolif': 'BENPORATH_PROLIFERATION'
}

hallmark_lib = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
c2_lib = gp.read_gmt("/home/download/pysc/input/c2.all.v2025.1.Hs.symbols.gmt")

final_gene_sets = {}
for name, genes in custom_signatures.items():
    final_gene_sets[name] = [g for g in genes if g in adata_csc.var_names]

for custom_name, db_name in msigdb_targets.items():
    genes = hallmark_lib.get(db_name, c2_lib.get(db_name, []))
    final_gene_sets[custom_name] = [g for g in genes if g in adata_csc.var_names]

# 2.  AUCell
counts = adata_csc.layers["raw_counts"] if "raw_counts" in adata_csc.layers else adata_csc.X
if hasattr(counts, "toarray"):
    df_counts = pd.DataFrame(counts.toarray(), index=adata_csc.obs_names, columns=adata_csc.var_names)
else:
    df_counts = pd.DataFrame(counts, index=adata_csc.obs_names, columns=adata_csc.var_names)

signatures = [GeneSignature(name, genes) for name, genes in final_gene_sets.items() if genes]
auc_mtx = aucell(df_counts, signatures, num_workers=16)

adata_csc.obs.drop(columns=[c for c in auc_mtx.columns if c in adata_csc.obs.columns], inplace=True, errors='ignore')
adata_csc.obs = pd.concat([adata_csc.obs, auc_mtx], axis=1)

# 3.  Z-score
pathway_structure = {
    "Cell Identity": ['Hepatocytes', 'BEC', 'Bipotent'],
    "CSC Features": ['JHEP_CSC_Like', 'JHEP_CSC_Con', 'Review_HCC_CSC', 'NPJ_CSC', 'CURATED_Stemness_markers'],
    "Proliferation": ['CellCycle_S_Phase', 'CellCycle_G2M_Phase', 'Proliferation_MYC', 'BenPorath_Prolif'],
    "Signaling Pathways": ['WNT_Signaling_Online', 'NOTCH_Signaling', 'HEDGEHOG_Signaling', 'TGF_BETA_Signaling',
                           'IL6_JAK_STAT3', 'JNK_MAPK_Signaling', 'ERK_Pathway', 'MET_Pathway', 'TLR4_Signaling', 'EMT_Signaling']
}

colors = {
    "Cell Identity": "# 00A087",      #
    "CSC Features": "# E64B35",       #
    "Proliferation": "# 3C5488",      #
    "Signaling Pathways": "# F39C12"  #
}

valid_keys = [k for sublist in pathway_structure.values() for k in sublist if k in adata_csc.obs.columns]

group_mean = adata_csc.obs[valid_keys + [cluster_key]].groupby(cluster_key).mean()
group_z = group_mean.apply(lambda x: (x - x.mean()) / (x.std() + 1e-9)).fillna(0).T

row_colors = pd.Series(index=valid_keys, dtype=object)
for cat, terms in pathway_structure.items():
    for term in terms:
        if term in row_colors.index: row_colors[term] = colors[cat]

custom_cmap = LinearSegmentedColormap.from_list("deepblue_red", ["#08306b", "#ffffff", "#e31a1c"])

figsize = (len(group_z.columns) * 0.8 + 3, len(valid_keys) * 0.4)

g = sns.clustermap(
    group_z, row_cluster=False, col_cluster=False,
    cmap=custom_cmap, center=0, vmin=-2.0, vmax=2.0,
    linewidths=0.5, linecolor='black', row_colors=row_colors, 
    figsize=figsize, cbar_pos=(1.05, 0.2, 0.03, 0.4)
)

g.ax_row_colors.set_xticks([])
g.ax_heatmap.yaxis.set_ticks_position('left')
g.ax_heatmap.yaxis.set_label_position('left')
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.set_xlabel(f"Tumor Clusters ({cluster_key})", fontsize=12, fontweight='bold', labelpad=10)

g.ax_heatmap.tick_params(axis='both', which='both', length=0)
g.cax.tick_params(axis='both', which='both', length=0)

plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=11)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=11)

patches = [mpatches.Patch(color=c, label=l) for l, c in colors.items()]
g.ax_heatmap.legend(handles=patches, title="Categories", bbox_to_anchor=(1.05, 1.05), loc='upper left', frameon=False)
g.cax.set_title("Z-score", fontsize=10, pad=10)

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'csc_Pathway_Heatmap.pdf')

plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

plt.show()

# --- 1e: Stacked violin plot ---
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np

markers_dict = {
    'Hepatocytes': ['ALB', "ASGR1", 'APOA2', 'APOC3', "CYP3A4", "GLUL"],
    "BEC": ["KRT7", "CD24", "MUC5B", "MUC3A", "SPP1"],
    'Bipotent': ['KRT19', 'TNFRSF12A', "SOX9", "HNF4A", "SERPINE1", "FKBP5",
                 "FN1", "GPC6", "IL18", "BICC1", "CDH6", "CREB5", "HDAC9",
                 "FGFR2", "SOX4", "SOX6", "PTCHD4"],
    "Malignant": ['APOH', 'GPC3', "AFP"],
    "Stemness": ['EPCAM', 'PROM1', "CD47", "ANPEP", "THY1"],
    "Immune": ['PTPRC'],
}

cell_order = [
    'Hepatocytes (Metabolic-H)',
    'Hepatocytes (Metabolic-L)',
    # ---  (CSC) ---
    'CSC(progenitor_cycling)',
    'CSC(progenitor_quiet)',
    'CSC',
    'Cholangiocytes',
    'Hepatocytes',
    'Bipotent Progenitor',
]

group_annotations = [
    (0, 1, 'Tumor', '#E64B35'),       # Hepatocytes Metabolic-H/L
    (2, 4, 'Tumor (CSC)', '# E64B35'),  # CSC
    (5, 7, 'Normal', '#4DBBD5'),       # Cholangiocytes / Hepatocytes / Bipotent
]

sns.set_theme(style='whitegrid', font_scale=1.1)
num_genes = sum(len(v) for v in markers_dict.values())
fig_width = num_genes * 0.45 + 4
fig_height = len(cell_order) * 0.9 + 2

# =================  category  =================
adata_subset.obs['final_type'] = adata_subset.obs['final_type'].astype('category')
adata_subset.obs['final_type'] = adata_subset.obs['final_type'].cat.reorder_categories(cell_order)

axes = sc.pl.stacked_violin(
    adata_subset,
    markers_dict,
    groupby='final_type',
    standard_scale='var',
    linewidth=0.8,
    colorbar_title='Standardized\nExpression',
    figsize=(fig_width, fig_height),
    return_fig=True,
)

ax_dict = axes.get_axes()
print("Available axes keys:", list(ax_dict.keys()))

ax_main = ax_dict.get('mainplot_ax') or ax_dict.get('main_plot_ax')

ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=90, ha='center')
ax_main.tick_params(axis='y', labelsize=11)
ax_main.tick_params(axis='x', labelsize=9)

# =================  /  =================
fig = ax_main.get_figure()
bbox = ax_main.get_position()

for start_idx, end_idx, label, color in group_annotations:
    n_total = len(cell_order)
    y_top = 1.0 - start_idx / n_total
    y_bot = 1.0 - (end_idx + 1) / n_total
    y_center = (y_top + y_bot) / 2.0

    fig_y_top = bbox.y0 + y_top * bbox.height
    fig_y_bot = bbox.y0 + y_bot * bbox.height
    fig_y_ctr = (fig_y_top + fig_y_bot) / 2.0
    fig_x = bbox.x1 + 0.01

    bracket_ax = fig.add_axes([0, 0, 1, 1], facecolor='none', zorder=10)
    bracket_ax.set_xlim(0, 1)
    bracket_ax.set_ylim(0, 1)
    bracket_ax.axis('off')

    bracket_ax.plot(
        [fig_x, fig_x], [fig_y_bot, fig_y_top],
        color=color, linewidth=2.5, solid_capstyle='round',
        transform=fig.transFigure, clip_on=False
    )
    tick_len = 0.008
    for yy in [fig_y_bot, fig_y_top]:
        bracket_ax.plot(
            [fig_x - tick_len, fig_x], [yy, yy],
            color=color, linewidth=2.5,
            transform=fig.transFigure, clip_on=False
        )
    bracket_ax.text(
        fig_x + 0.015, fig_y_ctr, label,
        rotation=270, va='center', ha='center',
        fontsize=11, fontweight='bold', color=color,
        transform=fig.transFigure
    )

plt.savefig("/home/download/csc_article/fig2/stacked_violin_markers.pdf",
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# --- 1f: RNA velocity streamlines (all / tumor / normal) ---
import scanpy as sc
import scvelo as scv
import pandas as pd
import anndata as ad
import os
import glob
import time
adata = sc.read_h5ad("/data1/velocity/liver_cancer_velocity_analyzed.h5ad")
import scanpy as sc

source_adata = sc.read_h5ad('/data2/adata/adata_subset_0214.h5ad')

# 2.  index (Barcode)
# adata.obs['final_type'] = source_adata.obs['final_type']
adata.obs['final_type'] = adata.obs.index.map(source_adata.obs['final_type'])

nan_count = adata.obs['final_type'].isna().sum()
if nan_count > 0:
    pass
else:
    pass

del source_adata
import scvelo as scv
import matplotlib.pyplot as plt
import numpy as np
import os

color_dict = {
    "Hepatocytes":              "# 4E79A7",   #
    "Hepatocytes (Metabolic-H)":"# 59A14F",   #
    "Hepatocytes (Metabolic-L)":"# 8CD17D",   #
    "Hepatocytes (Stem)":       "# 499894",   #
    "Hepatocytes (cycling)":    "# 86BCB6",   #
    "Hepatocytes (Quiescent)":  "# A0CBE8",   #

    # --- /CSC  (Tableau ) ---
    "CSC(progenitor_cycling)":  "# E15759",   #
    "CSC(progenitor_quiet)":    "# D4A6C8",   #
    "CSC":                      "# F28E2B",   #

    "Cholangiocytes":           "# B07AA1",   #
    "Biliary-like":             "# FF9D9A",   #
    "Bipotent Progenitor":      "# EDC948",   #
}
def set_colors(ad, col, cdict):
    """Set color palette in adata.uns"""
    cats = ad.obs[col].astype('category').cat.categories.tolist()
    ad.uns[f'{col}_colors'] = [cdict.get(c, '#999999') for c in cats]

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)

stream_params = dict(
    basis='umap',
    color='final_type',
    smooth=1.5,
    min_mass=3.0,
    density=1.2,
    cutoff_perc=10,
    linewidth=1.2,
    legend_loc='right margin',
    show=False,
)

# 1:  RNA Velocity Stream
mask = np.all(np.isfinite(adata.obsm['X_umap']), axis=1)
if 'velocity_umap' in adata.obsm:
    mask &= np.all(np.isfinite(adata.obsm['velocity_umap']), axis=1)
adata_clean = adata[mask].copy()
set_colors(adata_clean, 'final_type', color_dict)

scv.pl.velocity_embedding_stream(
    adata_clean,
    **stream_params,
    title='RNA Velocity (Stream) - All Cells',
)
plt.savefig(os.path.join(output_dir, 'RNA_Velocity_Stream_All.svg'),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# 2:  RNA Velocity Stream
set_colors(adata_tumor, 'final_type', color_dict)

scv.pp.moments(adata_tumor, n_pcs=None, n_neighbors=None)
scv.tl.velocity(adata_tumor, mode='stochastic')
scv.tl.velocity_graph(adata_tumor)

scv.pl.velocity_embedding_stream(
    adata_tumor,
    **stream_params,
    title='RNA Velocity (Stream) - Tumor Cells',
)
plt.savefig(os.path.join(output_dir, 'RNA_Velocity_Stream_Tumor.svg'),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# 3:  RNA Velocity Stream
set_colors(adata_final, 'final_type', color_dict)

scv.pl.velocity_embedding_stream(
    adata_final,
    **stream_params,
    title='RNA Velocity (Stream) - Normal Cells',
)
plt.savefig(os.path.join(output_dir, 'RNA_Velocity_Stream_Normal.svg'),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# --- 1f: Normal cell velocity (coordinate porting + scVelo) ---
import scanpy as sc
import matplotlib.pyplot as plt
import os

# 1.  (Porting Coordinates & Manifold)
ref_adata = sc.read_h5ad('/data2/adata/normal_annoted.h5ad')

common_cells = adata_normal.obs_names.intersection(ref_adata.obs_names)

adata_normal = adata_normal[common_cells].copy()

idx = ref_adata.obs_names.get_indexer(common_cells)
for key in ['X_umap', 'X_pca']:
    if key in ref_adata.obsm:
        adata_normal.obsm[key] = ref_adata.obsm[key][idx]

if 'final_type_colors' in ref_adata.uns:
    adata_normal.uns['final_type_colors'] = ref_adata.uns['final_type_colors']

if 'neighbors' in ref_adata.uns:
    adata_normal.uns['neighbors'] = ref_adata.uns['neighbors']
    adata_normal.obsp['distances'] = ref_adata.obsp['distances'][idx][:, idx]
    adata_normal.obsp['connectivities'] = ref_adata.obsp['connectivities'][idx][:, idx]

plt.gcf().set_dpi(72)

sc.pl.umap(
    adata_normal, 
    color=['final_type'] + [col for col in adata_normal.obs.columns if '_norm' in col][:1], 
    ncols=2, 
    frameon=False,
    wspace=0.3
)

output_h5ad = "/home/download/csc_article/fig2/normal_velocity_atlas.h5ad"
os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata_normal.write(output_h5ad)

del ref_adata
import scvelo as scv

sc.pp.neighbors(adata_normal, n_neighbors=20, n_pcs=30)

# 2.  Moments
scv.pp.moments(adata_normal, n_pcs=None, n_neighbors=None)

scv.tl.velocity(adata_normal, mode='stochastic')

scv.tl.velocity_graph(adata_normal)

scv.pl.velocity_embedding_stream(
    adata_normal, 
    basis='umap', 
    color='final_type', # 'leiden_mp_norm'
    palette='Set1',
    title='RNA Velocity in MP-space (Normal Cells)',
    legend_loc='right margin',
    dpi=100
)

output_path = "/home/download/csc_article/fig1/RNA_Velocity_Normal_Cells.pdf"
plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)
import os
import matplotlib.pyplot as plt
import scvelo as scv

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, 'ATLAS_RNA_Velocity_Normal_Cells.pdf')
svg_path = os.path.join(output_dir, 'ATLAS_RNA_Velocity_Normal_Cells.svg')

scv.pl.velocity_embedding_stream(
    adata_final, 
    basis='umap', 
    color='final_type',
    palette='tab20',
    
    smooth=1.5,
    min_mass=3.0,
    density=1.2,
    cutoff_perc=10,
    linewidth=1.2,
    
    legend_loc='right margin',
    title='RNA Velocity (Stream) - Smoothed',
    show=False
)

try:
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', transparent=True)
except ValueError as e:
    plt.savefig(svg_path, dpi=300, bbox_inches='tight', transparent=True)

plt.gcf().set_dpi(72)
plt.show()

# --- 1f: Tumor cell velocity (MP-space + scVelo) ---
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import pickle
from sklearn.preprocessing import MinMaxScaler
MP_PATH="/home/download/pysc/input/subset_nmf_output/final_gavish_mps.pkl"
target_clusters = ["tumor"]
# adata  loom  spliced/unspliced
adata_tumor = adata[adata.obs['malignant'].isin(target_clusters)].copy()

# --- B.  MP  (0-1 ) ---
# --- B.  MP  (0-1 ) ---
adata_tumor.X = adata_tumor.layers["raw_counts"].copy()

# 【】 tumorize_total  normalize_total
sc.pp.normalize_total(adata_tumor, target_sum=1e4) 
sc.pp.log1p(adata_tumor) 

with open(MP_PATH, 'rb') as f:
    mp_gene_sets = pickle.load(f)

norm_mp_cols = []
scaler = MinMaxScaler(feature_range=(0, 1))

for mp_name, genes in mp_gene_sets.items():
    if "Gavish_MP_" in mp_name and "_6" not in mp_name:
        valid_genes = [g for g in genes if g in adata_tumor.var_names]
        if len(valid_genes) > 0:
            sc.tl.score_genes(adata_tumor, gene_list=valid_genes, score_name=mp_name)
            raw_scores = adata_tumor.obs[mp_name].values.reshape(-1, 1)
            norm_col = f"{mp_name}_norm"
            adata_tumor.obs[norm_col] = scaler.fit_transform(raw_scores).flatten()
            norm_mp_cols.append(norm_col)
adata_tumor.obsm['X_gavish_mp_norm'] = adata_tumor.obs[norm_mp_cols].values

sc.pp.neighbors(adata_tumor, use_rep='X_gavish_mp_norm', n_neighbors=15, metric='euclidean')

sc.tl.leiden(adata_tumor, key_added='leiden_mp_norm', resolution=0.3)

sc.tl.umap(adata_tumor)

# 1.  Moments ( neighbors )
scv.pp.moments(adata_tumor, n_pcs=None, n_neighbors=None)

scv.tl.velocity(adata_tumor, mode='stochastic')
scv.tl.velocity_graph(adata_tumor)

scv.pl.velocity_embedding_stream(
    adata_tumor, 
    basis='umap', 
    color='final_type', 
    title='RNA Velocity in MP-space (tumor Cells)',
    legend_loc='right margin'
)
import os
import matplotlib.pyplot as plt
import scvelo as scv

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, 'RNA_Velocity_Tumor_Cells.pdf')
svg_path = os.path.join(output_dir, 'RNA_Velocity_Tumor_Cells.svg')

scv.pl.velocity_embedding_stream(
    adata_tumor, 
    basis='umap', 
    color='final_type', 
    palette='tab20',
    title='RNA Velocity in MP-space (tumor Cells)',
    legend_loc='right margin',
    
    # smooth=1.5, min_mass=3.0, density=1.2, cutoff_perc=10, linewidth=1.2,
    
    show=False
)

try:
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', transparent=True)
except ValueError as e:
    plt.savefig(svg_path, dpi=300, bbox_inches='tight', transparent=True)

plt.gcf().set_dpi(72)
plt.show()

# --- 1f: Export / import tumor UMAP coordinates ---
import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

output_dir = '/home/download/csc_article/fig2/'
os.makedirs(output_dir, exist_ok=True)

umap_tumor_df = pd.DataFrame(
    adata_tumor.obsm['X_umap'],
    index=adata_tumor.obs_names,
    columns=['UMAP_1', 'UMAP_2']
)
umap_tumor_df.to_csv("/home/download/csc_article/fig2/tumor_umap_coords.csv")
import pandas as pd

umap_ext = pd.read_csv("/home/download/csc_article/fig2/tumor_umap_coords.csv", index_col=0)

common_cells = adata_tumor.obs_names.intersection(umap_ext.index)

# adata_tumor
adata_tumor = adata_tumor[common_cells].copy()

adata_tumor.obsm['X_umap'] = umap_ext.loc[common_cells].values
