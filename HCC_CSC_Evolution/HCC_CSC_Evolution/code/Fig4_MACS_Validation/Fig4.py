#!/usr/bin/env python
# Figure 4 — EpCAM-enriched clinical samples validate the evolutionary model
# 4a: RNA velocity (all parenchymal / normal / tumor / CSC subtypes)

import scanpy as sc
import scvelo as scv
import scvi
import anndata as ad
import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from scipy.sparse import issparse

scv.settings.verbosity = 3
scv.set_figure_params('scvelo', dpi=150, color_map='viridis')

OUT_DIR = "/home/download/csc_article/fig3/macs/velocity/"
os.makedirs(OUT_DIR, exist_ok=True)
scv.settings.figdir = OUT_DIR
sc.settings.figdir = OUT_DIR

COLOR_DICT = {
    "Hepatocytes":              "#4E79A7",
    "Hepatocytes (Metabolic High)":"#59A14F",
    "Hepatocytes (Metabolic Low)":"#8CD17D",
    "Hepatocytes (Cycling)":    "#86BCB6",
    "Cholangiocyte":            "#B07AA1",
    "Cholangiocytes":           "#B07AA1",
    "Biliary-like":             "#FF9D9A",
    "Bipotent":                 "#EDC948",
    "Bipotent Progenitor":      "#EDC948",
    "Progenitor-like":          "#DBDB8D",
    "CSC":                      "#F28E2B",
    "CSC(progenitor prolif)":   "#E15759",
    "CSC(progenitor_prolif)":   "#E15759",
    "CSC(progenitor quiet)":    "#D4A6C8",
    "CSC(progenitor_quiet)":    "#D4A6C8",
    "Hepatocytes (Metabolic-H)":"#59A14F",
    "Hepatocytes (Metabolic-L)":"#8CD17D",
    "Hepatocytes (Stem)":       "#499894",
    "Hepatocytes (cycling)":    "#86BCB6",
    "Hepatocytes (Quiescent)":  "#A0CBE8",
    "CSC(progenitor_cycling)":  "#E15759",
}

def apply_colors(adata, key):
    """Set colors from COLOR_DICT for adata.uns['{key}_colors']"""
    cats = adata.obs[key].astype('category').cat.categories.tolist()
    colors = [COLOR_DICT.get(c, '#999999') for c in cats]
    adata.uns[f'{key}_colors'] = colors
    adata.obs[key] = adata.obs[key].astype('category')

adata_normal = sc.read("/home/download/csc_article/monocle4bile/adata_normal_MACS.h5ad")  # normal UMAP + manual_annotation
adata_csc_full = sc.read("/data2/MACS/downstream/CSC_annoted.h5ad")   # UMAP + cell_type
adata_csc_sub  = sc.read("/data2/MACS/downstream/CSC.h5ad")          # CSC  UMAP
adata_velo     = sc.read("/data2/MACS/downstream/normal_velocity_final.h5ad")  # velocity

# normal
mask = adata_normal.obs['manual_annotation'].notna()
adata_normal = adata_normal[mask].copy()

print(f"adata_normal:   {adata_normal.n_obs}")
print(f"adata_csc_full: {adata_csc_full.n_obs}")
print(f"adata_csc_sub:  {adata_csc_sub.n_obs}")
print(f"adata_velo:     {adata_velo.n_obs}")

# %% [2]  final_celltype
adata_all = adata_csc_full.copy()
adata_all.obs['final_celltype'] = adata_all.obs['cell_type'].astype(str)

# normal  manual_annotation
normal_in_all = adata_all.obs_names.intersection(adata_normal.obs_names)

if len(normal_in_all) > 100:
    adata_all.obs.loc[normal_in_all, 'final_celltype'] = (
        adata_normal.obs.loc[normal_in_all, 'manual_annotation'].astype(str).values
    )
else:
    def core_bc(bc): return re.sub(r'_\d+$', '', str(bc))
    normal_core = {core_bc(bc): bc for bc in adata_normal.obs_names}
    all_core = {core_bc(bc): bc for bc in adata_all.obs_names}
    for c in set(normal_core) & set(all_core):
        adata_all.obs.loc[all_core[c], 'final_celltype'] = (
            str(adata_normal.obs.loc[normal_core[c], 'manual_annotation'])
        )
    normal_in_all = pd.Index([all_core[c] for c in set(normal_core) & set(all_core)])

# tissue
adata_all.obs['tissue'] = 'Tumor'
adata_all.obs.loc[normal_in_all, 'tissue'] = 'Normal'

print("\nfinal_celltype:")
print(adata_all.obs['final_celltype'].value_counts())
print(f"\nNormal: {(adata_all.obs['tissue']=='Normal').sum()}, Tumor: {(adata_all.obs['tissue']=='Tumor').sum()}")

# %% [3]  loom  spliced/unspliced
DIR_LOOM = "/data1/velocity/MACS_velocity"
loom_info = {
    'T2': {'file': 'T2.loom', 'suffix': '_0'},
    'T1': {'file': 'T1.loom', 'suffix': '_1'},
    'P2': {'file': 'P2.loom', 'suffix': '_2'},
    'P1': {'file': 'P1.loom', 'suffix': '_3'},
}

loom_list = []
for smp, info in loom_info.items():
    fpath = os.path.join(DIR_LOOM, info['file'])
    if not os.path.exists(fpath):
        continue
    try:
        ldat = scv.read_loom(fpath, cache=True)
    except AttributeError:
        ldat = sc.read_loom(fpath)
    ldat.var_names_make_unique()
    new_bcs = []
    for bc in ldat.obs_names:
        match = re.search(r'[ACGT]{16}', str(bc))
        new_bcs.append(f"{match.group()}-1{info['suffix']}" if match else bc)
    ldat.obs_names = new_bcs
    loom_list.append(ldat)
    print(f"  {smp}: {ldat.n_obs} cells")

ldat_merged = ad.concat(loom_list, index_unique=None)
print(f"Loom total: {ldat_merged.n_obs}")

# %% [4] =====  velocity ( CSC_annoted  UMAP) =====
print("\n" + "=" * 60)
print("=" * 60)

adata_v_all = scv.utils.merge(adata_all, ldat_merged)
scv.pp.filter_and_normalize(adata_v_all, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_v_all, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_v_all, mode='stochastic')
scv.tl.velocity_graph(adata_v_all)

apply_colors(adata_v_all, 'final_celltype')
scv.pl.velocity_embedding_stream(
    adata_v_all, basis='umap', color='final_celltype',
    density=0.6, smooth=0.8, linewidth=1.0,
    legend_loc='right margin', frameon=False,
    title='All Parenchymal Cells',
    save="MACS_all_velocity.pdf"
)

# %% [5] ===== Normal velocity ( adata_normal_MACS  UMAP) =====
print("\n" + "=" * 60)
print("Normal velocity")
print("=" * 60)

normal_bcs = adata_all.obs_names[adata_all.obs['tissue'] == 'Normal']
adata_v_normal = adata_v_all[adata_v_all.obs_names.isin(normal_bcs)].copy()

# adata_normal_MACS  UMAP
common_n = adata_v_normal.obs_names.intersection(adata_normal.obs_names)

adata_v_normal = adata_v_normal[common_n].copy()
adata_v_normal.obsm['X_umap'] = adata_normal[common_n].obsm['X_umap'].copy()

sc.pp.neighbors(adata_v_normal, n_neighbors=30)
scv.pp.moments(adata_v_normal, n_pcs=30, n_neighbors=None)
scv.tl.velocity(adata_v_normal, mode='stochastic')
scv.tl.velocity_graph(adata_v_normal)

apply_colors(adata_v_normal, 'final_celltype')
scv.pl.velocity_embedding_stream(
    adata_v_normal, basis='umap', color='final_celltype',
    density=0.8, smooth=0.8, frameon=False,
    legend_loc='right margin',
    title='Normal Parenchymal (own UMAP)',
    save="MACS_normal_velocity.pdf"
)

# %% [6] ===== Tumor velocity (scVI → scANVI ) =====
print("\n" + "=" * 60)
print("=" * 60)

tumor_bcs = adata_all.obs_names[adata_all.obs['tissue'] == 'Tumor']
adata_v_tumor = adata_v_all[adata_v_all.obs_names.isin(tumor_bcs)].copy()

# counts
if 'counts' not in adata_v_tumor.layers:
    if 'raw_counts' in adata_v_tumor.layers:
        adata_v_tumor.layers['counts'] = adata_v_tumor.layers['raw_counts'].copy()
    elif 'spliced' in adata_v_tumor.layers:
        adata_v_tumor.layers['counts'] = adata_v_tumor.layers['spliced'].copy()
    else:
        adata_v_tumor.layers['counts'] = adata_v_tumor.X.copy()

ct = adata_v_tumor.layers['counts']
if issparse(ct):
    ct.data = np.rint(ct.data)
    adata_v_tumor.layers['counts'] = ct.astype(int)
else:
    adata_v_tumor.layers['counts'] = np.rint(ct).astype(int)

# HVG
adata_tumor_hvg = adata_v_tumor.copy()
sc.pp.highly_variable_genes(
    adata_tumor_hvg, layer='counts', n_top_genes=3000,
    flavor='seurat_v3', subset=True
)

# scVI
scvi.model.SCVI.setup_anndata(adata_tumor_hvg, layer='counts', batch_key='batch')
model_vi = scvi.model.SCVI(adata_tumor_hvg, n_layers=2, n_latent=20, gene_likelihood='nb')
model_vi.train(max_epochs=200, early_stopping=True)
adata_tumor_hvg.obsm['X_scVI'] = model_vi.get_latent_representation()

# scANVI
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model_vi, adata=adata_tumor_hvg,
    labels_key='final_celltype', unlabeled_category='Unknown'
)
scanvi_model.train(max_epochs=50, n_samples_per_label=100)
adata_tumor_hvg.obsm['X_scANVI'] = scanvi_model.get_latent_representation()

# UMAP
sc.pp.neighbors(adata_tumor_hvg, use_rep='X_scANVI', n_neighbors=30)
sc.tl.umap(adata_tumor_hvg)

adata_v_tumor.obsm['X_umap'] = np.zeros((adata_v_tumor.n_obs, 2))
adata_v_tumor.obsm['X_umap'][
    np.isin(adata_v_tumor.obs_names, adata_tumor_hvg.obs_names)
] = adata_tumor_hvg.obsm['X_umap']

keep = adata_v_tumor.obs_names.isin(adata_tumor_hvg.obs_names)
adata_v_tumor = adata_v_tumor[keep].copy()
adata_v_tumor.obsm['X_umap'] = adata_tumor_hvg[adata_v_tumor.obs_names].obsm['X_umap'].copy()

sc.pl.umap(adata_v_tumor, color='final_celltype', title='Tumor (scANVI UMAP)',
           frameon=False, legend_loc='right margin')

# velocity
scv.pp.moments(adata_v_tumor, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_v_tumor, mode='stochastic')
scv.tl.velocity_graph(adata_v_tumor)

apply_colors(adata_v_tumor, 'final_celltype')
scv.pl.velocity_embedding_stream(
    adata_v_tumor, basis='umap', color='final_celltype',
    density=0.8, smooth=0.8, frameon=False,
    legend_loc='right margin',
    title='Tumor Parenchymal (scANVI UMAP)',
    save="MACS_tumor_velocity.pdf"
)

# %% [7] ===== CSC  velocity ( CSC.h5ad  UMAP) =====
print("\n" + "=" * 60)
print("=" * 60)

csc_bcs = adata_csc_sub.obs_names
adata_v_csc = adata_v_all[adata_v_all.obs_names.isin(csc_bcs)].copy()

# CSC.h5ad  UMAP
common_csc = adata_v_csc.obs_names.intersection(adata_csc_sub.obs_names)

adata_v_csc = adata_v_csc[common_csc].copy()
adata_v_csc.obsm['X_umap'] = adata_csc_sub[common_csc].obsm['X_umap'].copy()

if 'csc_type' in adata_csc_sub.obs.columns:
    adata_v_csc.obs['csc_subtype'] = adata_csc_sub.obs.loc[common_csc, 'csc_type'].values
    color_csc = 'csc_subtype'
elif 'leiden_scvi_res0.3' in adata_csc_sub.obs.columns:
    adata_v_csc.obs['csc_cluster'] = adata_csc_sub.obs.loc[common_csc, 'leiden_scvi_res0.3'].values
    color_csc = 'csc_cluster'
else:
    color_csc = 'final_celltype'

sc.pp.neighbors(adata_v_csc, n_neighbors=min(30, len(common_csc) - 1))
scv.pp.moments(adata_v_csc, n_pcs=30, n_neighbors=None)
scv.tl.velocity(adata_v_csc, mode='stochastic')
scv.tl.velocity_graph(adata_v_csc)

apply_colors(adata_v_csc, color_csc)
scv.pl.velocity_embedding_stream(
    adata_v_csc, basis='umap', color=color_csc,
    density=1.0, smooth=0.8, frameon=False,
    legend_loc='right margin',
    title='CSC Subtypes (own UMAP)',
    save="MACS_CSC_velocity.pdf"
)

fig, axes = plt.subplots(2, 2, figsize=(20, 16))

scv.pl.velocity_embedding_stream(
    adata_v_all, basis='umap', color='final_celltype',
    ax=axes[0, 0], show=False, title='All Parenchymal',
    frameon=False, legend_loc='none', density=0.5, smooth=0.8
)

# Normal ( UMAP)
scv.pl.velocity_embedding_stream(
    adata_v_normal, basis='umap', color='final_celltype',
    ax=axes[0, 1], show=False, title='Normal (own UMAP)',
    frameon=False, legend_loc='none', density=0.8, smooth=0.8
)

# Tumor (scANVI UMAP)
scv.pl.velocity_embedding_stream(
    adata_v_tumor, basis='umap', color='final_celltype',
    ax=axes[1, 0], show=False, title='Tumor (scANVI UMAP)',
    frameon=False, legend_loc='none', density=0.8, smooth=0.8
)

# CSC ( UMAP)
scv.pl.velocity_embedding_stream(
    adata_v_csc, basis='umap', color=color_csc,
    ax=axes[1, 1], show=False, title='CSC (own UMAP)',
    frameon=False, legend_loc='right margin', density=1.0, smooth=0.8
)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "MACS_velocity_4panel.png"), dpi=300, bbox_inches='tight')
plt.show()

# %% [9]
adata_v_all.write(os.path.join(OUT_DIR, "MACS_all_velocity.h5ad"), compression='gzip')
adata_v_normal.write(os.path.join(OUT_DIR, "MACS_normal_velocity_final2.h5ad"), compression='gzip')
adata_v_tumor.write(os.path.join(OUT_DIR, "MACS_tumor_velocity.h5ad"), compression='gzip')
adata_v_csc.write(os.path.join(OUT_DIR, "MACS_CSC_velocity.h5ad"), compression='gzip')

# --- 4c: SNP-transcriptome pseudotime data preparation (Python) ---

#!/usr/bin/env python3
"""
MACS per-patient pipeline:
  1. Read from CSC_annoted.h5ad, use predicted_cell_type
  2. Compute MP scores
  3. MP-based leiden + UMAP
  4. Integrate SNP mutation info
  5. Export Matrix / Meta / UMAP for Monocle3
"""
import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import re
import os
import scipy.sparse as sp
from sklearn.preprocessing import MinMaxScaler

# ║                                                      ║

DATA_PATH = "/data2/MACS/downstream/CSC_annoted.h5ad"
MP_PATH   = "/home/download/pysc/input/subset_nmf_output/final_gavish_mps.pkl"
SNP_BASE  = "/home/download/tbsp"
OUT_DIR   = "/home/download/tbsp/MACS_MP_output"
os.makedirs(OUT_DIR, exist_ok=True)

PATIENTS = {
    "MACS_Patient1": ["P1", "T1"],
    "MACS_Patient2": ["P2", "T2"],
}

# ║                MP                               ║

with open(MP_PATH, 'rb') as f:
    mp_gene_sets = pickle.load(f)
target_mps = [k for k in mp_gene_sets if "Gavish_MP_" in k and k.split("_")[-1] != '6']
target_mps.sort(key=lambda x: int(x.split("_")[-1]))

# ║                                                 ║

adata = sc.read(DATA_PATH)

s_col = 'sample' if 'sample' in adata.obs.columns else 'batch'

# raw_counts
if 'raw_counts' not in adata.layers:
    if adata.X.max() > 50:
        adata.layers['raw_counts'] = adata.X.copy()
    elif adata.raw is not None:
        adata.layers['raw_counts'] = adata.raw.X.copy()
USE_LAYER = 'raw_counts' if 'raw_counts' in adata.layers else None

# ║             MP                         ║

for tag, samples in PATIENTS.items():
    print(f"\n{'='*60}")
    print(f"  {tag}: samples={samples}")
    print(f"{'='*60}")

    mask = adata.obs[s_col].isin(samples)
    ad_pat = adata[mask].copy()

    ad_pat.obs['CellType'] = ad_pat.obs['predicted_cell_type'].astype(str)

    adata_score = ad_pat.copy()
    if USE_LAYER and USE_LAYER in adata_score.layers:
        adata_score.X = adata_score.layers[USE_LAYER].copy()
        sc.pp.normalize_total(adata_score, target_sum=1e4)
        sc.pp.log1p(adata_score)
    elif adata_score.X.max() > 20:
        sc.pp.log1p(adata_score)

    scaler = MinMaxScaler()
    norm_cols = []
    for mp_name in target_mps:
        genes = mp_gene_sets[mp_name]
        valid = [g for g in genes if g in adata_score.var_names]
        if len(valid) > 3:
            sc.tl.score_genes(adata_score, gene_list=valid, score_name=mp_name,
                              ctrl_size=50, use_raw=False)
            col = f"{mp_name}_norm"
            ad_pat.obs[col] = scaler.fit_transform(
                adata_score.obs[mp_name].values.reshape(-1, 1)
            ).flatten()
            norm_cols.append(col)
    del adata_score

    print("  MP-based neighbors + leiden + UMAP...")
    mp_matrix = ad_pat.obs[norm_cols].values
    ad_pat.obsm['X_gavish_mp'] = mp_matrix

    sc.pp.neighbors(ad_pat, use_rep='X_gavish_mp', n_neighbors=50, key_added='nbr_mp')
    sc.tl.leiden(ad_pat, resolution=0.3, key_added='leiden_mp', neighbors_key='nbr_mp')
    sc.tl.umap(ad_pat, neighbors_key='nbr_mp', min_dist=0.5, spread=1.0)

    sc.pl.umap(ad_pat, color=['CellType', 'leiden_mp'], frameon=False,
               title=[f'{tag} CellType (MP UMAP)', f'{tag} Leiden (MP)'],
               save=f"_{tag}_MP_UMAP_check.png")

    # 5.  barcode
    barcodes_core = [re.sub(r'_\d+$', '', str(bc)) for bc in ad_pat.obs_names]

    total_muts = np.zeros(ad_pat.n_obs)
    tbsp_dir = f"tbsp_output_{tag}"
    for candidate in [
        os.path.join(SNP_BASE, "tbsp_out", tbsp_dir.replace("Patient1","P1").replace("Patient2","P2"), "SNP_matrix.tsv"),
        os.path.join(SNP_BASE, "tbsp_out", tbsp_dir, "SNP_matrix.tsv"),
        os.path.join(SNP_BASE, tbsp_dir, "SNP_matrix.tsv"),
        os.path.join(SNP_BASE, tbsp_dir.replace("Patient1","P1").replace("Patient2","P2"), "SNP_matrix.tsv"),
    ]:
        if os.path.exists(candidate):
            snp_df = pd.read_csv(candidate, sep='\t', index_col=0)
            bc_to_idx = {bc: i for i, bc in enumerate(barcodes_core)}
            for sbc in set(snp_df.columns) & set(barcodes_core):
                total_muts[bc_to_idx[sbc]] = snp_df[sbc].sum()
            break
    else:
        pass

    # 7.  MP-based UMAP
    umap_df = pd.DataFrame(ad_pat.obsm['X_umap'], index=barcodes_core,
                           columns=["UMAP_1", "UMAP_2"])
    umap_df = umap_df[~umap_df.index.duplicated(keep='first')]
    umap_df.to_csv(os.path.join(OUT_DIR, f"MP_UMAP_{tag}.csv"))

    # 8.  Meta
    meta_df = pd.DataFrame({
        'CellType': ad_pat.obs['CellType'].values,
        'Total_Mutations': total_muts.astype(int),
    }, index=barcodes_core)
    if s_col in ad_pat.obs.columns:
        meta_df['sample'] = ad_pat.obs[s_col].values
    meta_df = meta_df[~meta_df.index.duplicated(keep='first')]
    meta_df.to_csv(os.path.join(OUT_DIR, f"MP_Meta_{tag}.csv"))

    # 9.  Matrix (HVG)
    if USE_LAYER:
        ad_tmp = ad_pat.copy()
        ad_tmp.X = ad_tmp.layers[USE_LAYER].copy()
        sc.pp.normalize_total(ad_tmp, target_sum=1e4)
        sc.pp.log1p(ad_tmp)
        sc.pp.highly_variable_genes(ad_tmp, n_top_genes=2000, flavor='seurat')
        ad_hvg = ad_pat[:, ad_tmp.var.highly_variable].copy()
        rna_mat = ad_hvg.layers[USE_LAYER]
        rna_arr = rna_mat.toarray() if sp.issparse(rna_mat) else np.array(rna_mat)
        mat_df = pd.DataFrame(rna_arr.T, index=ad_hvg.var_names, columns=barcodes_core)
        mat_df = mat_df.loc[:, ~mat_df.columns.duplicated(keep='first')]
        mat_df.to_csv(os.path.join(OUT_DIR, f"MP_Matrix_{tag}.csv"))
        print(f"  Matrix: {mat_df.shape[0]} genes × {mat_df.shape[1]} cells")
        del ad_tmp, ad_hvg

print(f"\n{'='*60}")
print(f"{'='*60}")
