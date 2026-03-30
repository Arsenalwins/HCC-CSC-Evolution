#!/usr/bin/env python
# Figure 5c supplement + Figure S6 — WGCNA analysis plots

#!/usr/bin/env python3
"""
WGCNA full plotting script - generate all PDFs from saved data
Main figures:
  Fig1A  Transcriptional Similarity Clustermap (ME Pearson correlation)
  Fig1B  Pathway-Module DotPlot (Hallmark + KEGG, 4 focus modules)

Supplementary figures:
  FigS_W1  Module Sizes (focus modules highlighted)
  FigS_W2  ME z-score Heatmap (all significant modules)
  FigS_W3  Focus Module Violins (4 modules x cell types)
  FigS_W4  Module Score on UMAP
  FigS_W5  Hub Gene Grouped DotPlot (incl. Cholangiocytes)
  FigS_W6  Cohen's d Effect Size Heatmap

Required files (from celltype_specific_wgcna_v2.ipynb + similarity_enrichment_final.ipynb):
  - gene_modules.csv
  - Table_ME_per_Metacell.csv
  - Table_All_Gene_kME.csv / hub_genes_per_module.csv
  - Enrichment/Focus4_Enrichment_All.csv

"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
import scipy.sparse
import gseapy as gp
import os, warnings
warnings.filterwarnings('ignore')

# Config
ADATA_PATH = "/data2/adata/adata_subset_0214.h5ad"
OUTPUT_DIR = "/home/download/csc_article/WGCNA_V2/"
FIG_DIR    = "/home/download/csc_article/fig4/wgcna/"
ENRICH_DIR = f"{OUTPUT_DIR}Enrichment/"
os.makedirs(FIG_DIR, exist_ok=True)

CELLTYPE_COL = "final_type"
RAW_LAYER    = "raw_counts"
FIG_FORMAT   = "pdf"
DROP_CELLTYPE = "Biliary-like"

FOCUS = {
    'darkred':       {'label': 'CSC-M1',  'group': 'CSC',      'color': '#E74C3C'},
    'purple':        {'label': 'CSC-M2',  'group': 'CSC',      'color': '#E67E22'},
    'darkorange':    {'label': 'BP-M1',   'group': 'Bipotent', 'color': '#2980B9'},
    'skyblue':       {'label': 'BP-M2',   'group': 'Bipotent', 'color': '#27AE60'},
}

CT_COLORS = {
    'CSC':                      '#D62728',
    'CSC(progenitor_cycling)':  '#FF7F0E',
    'CSC(progenitor_quiet)':    '#E377C2',
    'Bipotent Progenitor':      '#9467BD',
    'Cholangiocytes':           '#8C564B',
    'Hepatocytes':              '#1F77B4',
    'Hepatocytes (Stem)':       '#17BECF',
    'Hepatocytes (cycling)':    '#BCBD22',
    'Hepatocytes (Metabolic-H)':'#2CA02C',
    'Hepatocytes (Metabolic-L)':'#7F7F7F',
    'Hepatocytes (Quiescent)':  '#AEC7E8',
}

CT_ORDER = [
    'CSC', 'CSC(progenitor_cycling)', 'CSC(progenitor_quiet)',
    'Bipotent Progenitor', 'Cholangiocytes',
    'Hepatocytes (Stem)', 'Hepatocytes (cycling)', 'Hepatocytes',
    'Hepatocytes (Metabolic-H)', 'Hepatocytes (Metabolic-L)',
    'Hepatocytes (Quiescent)',
]

GENE_SET_LIBS = [
    'GO_Biological_Process_2023',
    'GO_Molecular_Function_2023',
    'KEGG_2021_Human',
    'MSigDB_Hallmark_2020',
    'Reactome_2022',
]

# Load data
print("=" * 60)
print("Loading data...")
print("=" * 60)

gene_modules = pd.read_csv(f"{OUTPUT_DIR}gene_modules.csv", index_col=0)

# ME data
me_path = f"{OUTPUT_DIR}Table_ME_per_Metacell.csv"
if not os.path.exists(me_path):
    me_path = f"{OUTPUT_DIR}me_df.csv"
me_df = pd.read_csv(me_path, index_col=0)
me_df = me_df[me_df['celltype'] != DROP_CELLTYPE].copy()
mod_cols = [c for c in me_df.columns if c != 'celltype']
celltypes = sorted(me_df['celltype'].unique())

# Fix CT_ORDER
CT_ORDER = [c for c in CT_ORDER if c in celltypes]
for ct in celltypes:
    if ct not in CT_ORDER:
        CT_ORDER.append(ct)
extra_pal = sns.color_palette('husl', 20)
for i, ct in enumerate(celltypes):
    if ct not in CT_COLORS:
        CT_COLORS[ct] = matplotlib.colors.to_hex(extra_pal[i % 20])

print(f"  ME: {me_df.shape}, {len(celltypes)} cell types, {len(mod_cols)} modules")

# Mean ME & z-score
mean_me = me_df.groupby('celltype')[mod_cols].mean().loc[CT_ORDER]
corr_matrix = mean_me.T.corr(method='pearson')
zscore_me = mean_me.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0, axis=0)

# Enrichment results
enrich_path = f"{ENRICH_DIR}Focus4_Enrichment_All.csv"
enrich_results = {}
if os.path.exists(enrich_path):
    df_enrich = pd.read_csv(enrich_path)
    for mod in FOCUS:
        sub = df_enrich[df_enrich['module'] == mod]
        if len(sub) > 0:
            enrich_results[mod] = sub
    print(f"  Enrichment loaded: {len(df_enrich)} rows, {len(enrich_results)} modules")
else:
    print("  Enrichment CSV not found, will compute...")

# Hub genes
hub_path = f"{OUTPUT_DIR}hub_genes_per_module.csv"
if not os.path.exists(hub_path):
    hub_path = f"{OUTPUT_DIR}Table_All_Gene_kME.csv"
hub_df = pd.read_csv(hub_path)

# Kruskal-Wallis + Specificity (for heatmap & filtering)
print("\nComputing module specificity...")
kw_results = []
for mod in mod_cols:
    groups = [me_df.loc[me_df['celltype'] == ct, mod].dropna().values for ct in celltypes]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2: continue
    H, p = stats.kruskal(*groups)
    kw_results.append({'module': mod, 'H': H, 'p': p})
kw_df = pd.DataFrame(kw_results)
if len(kw_df) > 0:
    _, kw_df['p_adj'], _, _ = multipletests(kw_df['p'], method='fdr_bh')
    kw_df['sig'] = kw_df['p_adj'] < 0.05
    sig_mods = kw_df[kw_df['sig']].sort_values('p_adj')['module'].tolist()
    if len(sig_mods) < 5:
        sig_mods = kw_df.nsmallest(min(15, len(kw_df)), 'p_adj')['module'].tolist()
    print(f"  Significant modules: {kw_df['sig'].sum()} / {len(kw_df)}")

# Cohen's d specificity
spec_records = []
for mod in mod_cols:
    for ct in celltypes:
        ct_vals = me_df.loc[me_df['celltype'] == ct, mod].values
        other_vals = me_df.loc[me_df['celltype'] != ct, mod].values
        if len(ct_vals) < 2 or len(other_vals) < 2: continue
        U, p = stats.mannwhitneyu(ct_vals, other_vals, alternative='greater')
        pooled_std = np.sqrt((ct_vals.std()**2 + other_vals.std()**2) / 2)
        d = (ct_vals.mean() - other_vals.mean()) / pooled_std if pooled_std > 0 else 0
        z = zscore_me.loc[ct, mod] if ct in zscore_me.index else 0
        spec_records.append({
            'module': mod, 'celltype': ct, 'cohens_d': d,
            'p_value': p, 'z_score': z
        })
spec_df = pd.DataFrame(spec_records)
if len(spec_df) > 0:
    _, spec_df['fdr'], _, _ = multipletests(spec_df['p_value'], method='fdr_bh')
    spec_df['is_specific'] = (spec_df['fdr'] < 0.05) & (spec_df['cohens_d'] > 0.5)

# Fig1A: Transcriptional Similarity Clustermap
print("\n" + "=" * 60)
print("Fig1A: Transcriptional Similarity Clustermap")
print("=" * 60)

row_colors = pd.Series({ct: CT_COLORS[ct] for ct in CT_ORDER}, name='Cell Type')
cmap_sim = LinearSegmentedColormap.from_list('sim', ['#2166AC', '#F7F7F7', '#B2182B'], N=256)

g = sns.clustermap(
    corr_matrix.loc[CT_ORDER, CT_ORDER],
    method='ward',
    cmap=cmap_sim,
    vmin=-1, vmax=1, center=0,
    linewidths=0.4, linecolor='white',
    figsize=(7, 6.5),
    row_colors=row_colors, col_colors=row_colors,
    dendrogram_ratio=(0.15, 0.15),
    cbar_pos=(0.02, 0.82, 0.03, 0.12),
    cbar_kws={'label': 'Pearson r', 'ticks': [-1, -0.5, 0, 0.5, 1]},
    annot=True, fmt='.2f', annot_kws={'size': 7},
    tree_kws={'linewidths': 0.8, 'colors': '#333333'},
)
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.tick_params(axis='both', length=0)
g.fig.suptitle('Transcriptional Similarity\n(Module Eigengene Correlation)',
               fontsize=12, fontweight='bold', y=1.02)
g.savefig(f"{FIG_DIR}Fig1A_Transcriptional_Similarity.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
plt.close()
print("  -> Fig1A_Transcriptional_Similarity.pdf")

# Fig1B: Pathway-Module DotPlot (Hallmark + KEGG)
print("\n" + "=" * 60)
print("Fig1B: Pathway-Module DotPlot")
print("=" * 60)

# Run enrichment if not loaded
if not enrich_results:
    print("  Running Enrichr for 4 focus modules...")
    os.makedirs(ENRICH_DIR, exist_ok=True)
    for mod_name, info in FOCUS.items():
        print(f"    {info['label']}...")
        mod_genes = gene_modules[gene_modules['moduleColors'] == mod_name].index.tolist()
        try:
            enr = gp.enrichr(gene_list=mod_genes, gene_sets=GENE_SET_LIBS,
                             organism='human', outdir=None, cutoff=0.1, no_plot=True)
            res = enr.results.copy()
            res['module'] = mod_name
            res['label'] = info['label']
            enrich_results[mod_name] = res
        except Exception as e:
            print(f"    ERROR: {e}")
    if enrich_results:
        df_enrich = pd.concat(enrich_results.values(), ignore_index=True)
        df_enrich.to_csv(f"{ENRICH_DIR}Focus4_Enrichment_All.csv", index=False)
        print(f"  Saved {len(df_enrich)} enrichment results")

if enrich_results:
    summary_libs = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human']
    mods_ordered = list(FOCUS.keys())

    sel_terms = []
    for mod in mods_ordered:
        if mod not in enrich_results: continue
        sub = enrich_results[mod]
        sub = sub[sub['Gene_set'].isin(summary_libs)]
        sel_terms.extend(sub.sort_values('Adjusted P-value').head(5)['Term'].tolist())

    seen = set()
    unique_terms = [t for t in sel_terms if not (t in seen or seen.add(t))]

    dot_data = []
    for mod in mods_ordered:
        if mod not in enrich_results: continue
        sub = enrich_results[mod]
        for term in unique_terms:
            row = sub[sub['Term'] == term]
            if len(row) > 0:
                r = row.iloc[0]
                parts = r['Overlap'].split('/')
                gr = int(parts[0]) / int(parts[1]) if int(parts[1]) > 0 else 0
                dot_data.append({'module': mod, 'term': term,
                    'neg_log_p': -np.log10(max(r['Adjusted P-value'], 1e-50)),
                    'gene_ratio': gr, 'fdr': r['Adjusted P-value']})
            else:
                dot_data.append({'module': mod, 'term': term,
                    'neg_log_p': 0, 'gene_ratio': 0, 'fdr': 1.0})

    dot_df = pd.DataFrame(dot_data)
    dot_df['term_clean'] = dot_df['term'].apply(
        lambda x: x.split('(GO:')[0].split('(R-HSA')[0].split(' WP')[0].strip()[:48])

    term_order = (dot_df.groupby('term_clean')['neg_log_p'].max()
                  .sort_values(ascending=True).index.tolist())
    mod_labels = [FOCUS[m]['label'] for m in mods_ordered]
    mod_to_x = {m: i for i, m in enumerate(mods_ordered)}
    term_to_y = {t: i for i, t in enumerate(term_order)}

    fig, ax = plt.subplots(figsize=(5.5, max(5, len(term_order) * 0.3)))
    xs, ys, sizes, colors = [], [], [], []
    for _, row in dot_df.iterrows():
        if row['term_clean'] in term_to_y and row['neg_log_p'] > 0:
            xs.append(mod_to_x[row['module']])
            ys.append(term_to_y[row['term_clean']])
            sizes.append(row['gene_ratio'] * 1800 + 15)
            colors.append(row['neg_log_p'])

    if colors:
        sc_p = ax.scatter(xs, ys, s=sizes, c=colors, cmap='YlOrRd',
                          edgecolors='#333', linewidths=0.3, alpha=0.9,
                          vmin=0, vmax=max(colors))
        cbar = plt.colorbar(sc_p, ax=ax, shrink=0.35, pad=0.15, aspect=15)
        cbar.set_label('$-\\log_{10}$(FDR)', fontsize=7.5)
        cbar.ax.tick_params(labelsize=7)

    ax.axvline(1.5, color='#AAA', ls='--', lw=0.6, alpha=0.6)
    ax.text(0.5, len(term_order) + 0.4, 'CSC', ha='center', fontsize=8,
            fontstyle='italic', color='#C0392B', fontweight='bold')
    ax.text(2.5, len(term_order) + 0.4, 'Bipotent', ha='center', fontsize=8,
            fontstyle='italic', color='#2C3E50', fontweight='bold')

    ax.set_xticks(range(len(mods_ordered)))
    ax.set_xticklabels(mod_labels, fontsize=8, fontweight='bold')
    ax.set_yticks(range(len(term_order)))
    ax.set_yticklabels(term_order, fontsize=7)
    ax.set_xlim(-0.5, len(mods_ordered) - 0.5)
    ax.set_ylim(-0.5, len(term_order) + 0.8)
    ax.grid(True, alpha=0.08, linestyle='-')
    ax.set_title('Pathway Enrichment: Focus Modules', fontsize=10, fontweight='bold', pad=15)

    for gr, lb in [(0.05, '5%'), (0.10, '10%'), (0.20, '20%')]:
        ax.scatter([], [], s=gr * 1800 + 15, c='#AAA', edgecolors='#333', linewidths=0.3, label=lb)
    ax.legend(title='Gene Ratio', loc='lower right', bbox_to_anchor=(1.38, 0),
              fontsize=7, title_fontsize=7.5, frameon=True, edgecolor='#CCC')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(f"{FIG_DIR}Fig1B_Pathway_DotPlot.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> Fig1B_Pathway_DotPlot.pdf")
else:
    print("  SKIP: no enrichment results")

# FigS-W1: Module Sizes
print("\n--- FigS-W1: Module Sizes ---")

counts = gene_modules['moduleColors'].value_counts()
counts = counts[counts.index != 'grey']

bar_colors = []
edge_widths = []
for m in counts.index:
    if m in FOCUS:
        bar_colors.append(FOCUS[m]['color'])
        edge_widths.append(2)
    else:
        try: matplotlib.colors.to_rgba(m); bar_colors.append(m)
        except: bar_colors.append('#AAAAAA')
        edge_widths.append(0.5)

fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.5), 4.5))
bars = ax.bar(range(len(counts)), counts.values, color=bar_colors,
              edgecolor='black', linewidth=edge_widths)
ax.set_xticks(range(len(counts)))
ax.set_xticklabels(counts.index, rotation=60, ha='right', fontsize=7)
ax.set_ylabel('Number of Genes', fontweight='bold')
ax.set_title('WGCNA Module Sizes', fontsize=13, fontweight='bold')
for i, (m, v) in enumerate(counts.items()):
    if m in FOCUS:
        ax.text(i, v + 2, f"{FOCUS[m]['label']}\n({v})", ha='center', fontsize=7,
                fontweight='bold', color=FOCUS[m]['color'])
sns.despine()
plt.tight_layout()
plt.savefig(f"{FIG_DIR}FigS_W1_Module_Sizes.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
plt.close()
print("  -> FigS_W1_Module_Sizes.pdf")

# FigS-W2: ME z-score Heatmap
print("\n--- FigS-W2: ME z-score Heatmap ---")

plot_data = zscore_me[sig_mods].T
if plot_data.shape[0] > 1 and plot_data.shape[1] > 1:
    row_link = linkage(plot_data.values, method='ward')
    col_link = linkage(plot_data.values.T, method='ward')
    plot_data = plot_data.iloc[leaves_list(row_link), leaves_list(col_link)]
vmax = min(max(abs(plot_data.values.min()), abs(plot_data.values.max())), 3)

fig, ax = plt.subplots(figsize=(max(6, len(CT_ORDER) * 0.8), max(4, len(sig_mods) * 0.45)))
sns.heatmap(plot_data, cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'ME (z-score)', 'shrink': 0.6}, ax=ax)
ax.set_xlabel('Cell Type', fontweight='bold')
ax.set_ylabel('Module', fontweight='bold')
ax.set_title('Cell-Type-Specific WGCNA Modules', fontsize=13, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=7)
for label in ax.get_yticklabels():
    mod = label.get_text()
    if mod in FOCUS:
        label.set_color(FOCUS[mod]['color'])
        label.set_fontweight('bold')
plt.tight_layout()
plt.savefig(f"{FIG_DIR}FigS_W2_ME_Heatmap.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
plt.close()
print(f"  -> FigS_W2_ME_Heatmap.pdf ({len(sig_mods)} modules)")

# FigS-W3: Focus Module Violins
print("\n--- FigS-W3: Focus Module Violins ---")

fig, axes = plt.subplots(2, 2, figsize=(14, 9))
axes = axes.flatten()
palette = sns.color_palette('Set2', len(celltypes))
ct_pal = dict(zip(celltypes, palette))

for idx, (mod, info) in enumerate(FOCUS.items()):
    ax = axes[idx]
    if mod not in me_df.columns:
        ax.text(0.5, 0.5, f'{info["label"]} not found', ha='center', transform=ax.transAxes)
        continue
    sns.violinplot(data=me_df, x='celltype', y=mod, ax=ax, order=CT_ORDER,
                   palette=ct_pal, inner=None, alpha=0.6, linewidth=0.8, cut=0)
    sns.stripplot(data=me_df, x='celltype', y=mod, ax=ax, order=CT_ORDER,
                  palette=ct_pal, size=2, alpha=0.6, jitter=True)
    ax.set_title(f'{info["label"]} ({mod})', fontsize=11, fontweight='bold', color=info['color'])
    ax.set_xlabel(''); ax.set_ylabel('Module Eigengene')
    ax.tick_params(axis='x', rotation=45)
    ax.axhline(0, color='grey', ls='--', lw=0.5)
plt.suptitle('Focus Module Eigengenes across Cell Types', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}FigS_W3_Focus_Violins.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
plt.close()
print("  -> FigS_W3_Focus_Violins.pdf")

# FigS-W4: Cohen's d Heatmap
print("\n--- FigS-W4: Cohen's d Heatmap ---")

if len(spec_df) > 0:
    pivot_d = spec_df.pivot_table(index='module', columns='celltype',
                                   values='cohens_d', aggfunc='first').fillna(0)
    # Reorder columns
    pivot_d = pivot_d[[c for c in CT_ORDER if c in pivot_d.columns]]
    if pivot_d.shape[0] > 1:
        pivot_d = pivot_d.iloc[leaves_list(linkage(pivot_d.values, method='ward'))]

    fig, ax = plt.subplots(figsize=(max(8, pivot_d.shape[1] * 0.8), max(5, pivot_d.shape[0] * 0.4)))
    sns.heatmap(pivot_d, cmap='YlGnBu', linewidths=0.3, linecolor='white',
                vmin=0, vmax=min(pivot_d.values.max(), 3),
                cbar_kws={'label': "Cohen's d", 'shrink': 0.6},
                ax=ax, annot=True, fmt='.1f', annot_kws={'size': 7})
    ax.set_title("Module Specificity (Effect Size)", fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel('Cell Type', fontweight='bold')
    ax.set_ylabel('Module', fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    for label in ax.get_yticklabels():
        mod = label.get_text()
        if mod in FOCUS:
            label.set_color(FOCUS[mod]['color'])
            label.set_fontweight('bold')
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}FigS_W4_Cohens_d.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W4_Cohens_d.pdf")

# FigS-W5: Module Score on UMAP
print("\n--- FigS-W5: UMAP Module Scores ---")

print("  Loading adata for UMAP...")
adata = sc.read_h5ad(ADATA_PATH)

if 'X_umap' in adata.obsm:
    adata_plot = adata.copy()
    if RAW_LAYER in adata_plot.layers:
        adata_plot.X = adata_plot.layers[RAW_LAYER].copy()
    sc.pp.normalize_total(adata_plot, target_sum=1e4)
    sc.pp.log1p(adata_plot)

    for mod, info in FOCUS.items():
        mod_genes = gene_modules[gene_modules['moduleColors'] == mod].index
        avail = [g for g in mod_genes if g in adata_plot.var_names]
        if len(avail) >= 3:
            sc.tl.score_genes(adata_plot, avail, score_name=f'score_{info["label"]}')

    score_cols = [f'score_{info["label"]}' for mod, info in FOCUS.items()
                  if f'score_{info["label"]}' in adata_plot.obs.columns]

    if score_cols:
        fig, axes = plt.subplots(1, len(score_cols), figsize=(4.5 * len(score_cols), 4))
        if len(score_cols) == 1: axes = [axes]
        for ax, sc_name in zip(axes, score_cols):
            sc.pl.umap(adata_plot, color=sc_name, ax=ax, show=False,
                       title=sc_name.replace('score_', ''), frameon=False, color_map='YlOrRd')
        plt.tight_layout()
        plt.savefig(f"{FIG_DIR}FigS_W5_UMAP_ModuleScores.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
        plt.close()
        print("  -> FigS_W5_UMAP_ModuleScores.pdf")
    del adata_plot
else:
    print("  SKIP: no X_umap")

# FigS-W6: Hub Gene DotPlot ( Cholangiocytes)
print("\n--- FigS-W6: Hub Gene DotPlot ---")

LNCRNA_PAT = r'^(LINC|AC\d|AL\d|AP\d|RN7SL|CTD-|CTC-|RP\d|KB-)'
gene_dict = {}
for mod, info in FOCUS.items():
    sub = hub_df[hub_df['module'] == mod].copy()
    if 'gene' in sub.columns:
        sub = sub[~sub['gene'].str.match(LNCRNA_PAT, na=False)]
        sub = sub.sort_values('kME', ascending=False).head(10)
        genes = [g for g in sub['gene'].tolist() if g in adata.var_names]
        if genes:
            gene_dict[info['label']] = genes

if gene_dict:
    key_types = [t for t in CT_ORDER if any(p in t for p in ['CSC', 'Bipotent', 'Cholangio', 'Hepatocyte'])]
    adata_sub = adata[adata.obs[CELLTYPE_COL].isin(key_types)].copy()
    if RAW_LAYER in adata_sub.layers:
        adata_sub.X = adata_sub.layers[RAW_LAYER].copy()
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    adata_sub.obs[CELLTYPE_COL] = pd.Categorical(
        adata_sub.obs[CELLTYPE_COL], categories=key_types, ordered=True)

    dp = sc.pl.dotplot(adata_sub, var_names=gene_dict, groupby=CELLTYPE_COL,
                       standard_scale='var', show=False, return_fig=True)
    dp.savefig(f"{FIG_DIR}FigS_W6_HubGene_DotPlot.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W6_HubGene_DotPlot.pdf")
    del adata_sub

del adata

# FigS-W7: Soft Threshold Selection (Scale-Free Topology)
print("\n--- FigS-W7: Soft Threshold Selection ---")

# Reconstruct expression matrix for power analysis
mc_path_candidates = [
    "/home/download/csc_article/WGCNA_V2/adata_metacell_merged.h5ad",
    f"{OUTPUT_DIR}adata_metacell.h5ad",
    f"{OUTPUT_DIR}metacell_adata.h5ad",
    f"{OUTPUT_DIR}CSC_PerCT_metacell.h5ad",
]
adata_mc = None
for mp in mc_path_candidates:
    if os.path.exists(mp):
        adata_mc = sc.read_h5ad(mp)
        print(f"  Loaded metacell: {mp}")
        break

if adata_mc is not None:
    adata_mc_w = adata_mc.copy()
    sc.pp.normalize_total(adata_mc_w, target_sum=1e4)
    sc.pp.log1p(adata_mc_w)
    sc.pp.highly_variable_genes(adata_mc_w, n_top_genes=5000)
    hvg = adata_mc_w.var[adata_mc_w.var['highly_variable']].index

    X = adata_mc_w[:, hvg].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    df_expr = pd.DataFrame(X, index=adata_mc_w.obs_names, columns=hvg)

    # Test powers 1-20
    powers = list(range(1, 21))
    r2_vals = []
    mean_k = []

    print("  Computing scale-free fit for powers 1-20...")
    for pw in powers:
        # Adjacency (signed): ((1+cor)/2)^power
        cor_mat = np.corrcoef(df_expr.values.T)  # gene x gene would be huge
        # For efficiency, sample 2000 genes
        n_sample = min(2000, df_expr.shape[1])
        np.random.seed(42)
        idx = np.random.choice(df_expr.shape[1], n_sample, replace=False)
        sub = df_expr.iloc[:, idx].values
        cor_sub = np.corrcoef(sub.T)
        adj = ((1 + cor_sub) / 2) ** pw
        np.fill_diagonal(adj, 0)

        # Connectivity
        k = adj.sum(axis=1)
        mean_k.append(k.mean())

        # Scale-free fit: log(p(k)) ~ log(k)
        hist, bin_edges = np.histogram(k, bins=30)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mask = hist > 0
        if mask.sum() > 2:
            log_k = np.log10(bin_centers[mask])
            log_pk = np.log10(hist[mask] / hist[mask].sum())
            # Linear fit
            slope, intercept = np.polyfit(log_k, log_pk, 1)
            ss_res = np.sum((log_pk - (slope * log_k + intercept)) ** 2)
            ss_tot = np.sum((log_pk - log_pk.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            # Sign: scale-free should have negative slope
            r2_vals.append(-np.sign(slope) * r2)
        else:
            r2_vals.append(0)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel 1: Scale-free fit
    ax1.plot(powers, r2_vals, 'o-', color='#D62728', linewidth=1.5, markersize=5)
    ax1.axhline(0.85, color='red', ls='--', lw=0.8, alpha=0.6)
    ax1.set_xlabel('Soft Threshold (Power)', fontweight='bold')
    ax1.set_ylabel('Scale Free Topology\nModel Fit (signed R²)', fontweight='bold')
    ax1.set_title('Scale-Free Fit Index', fontweight='bold')
    for i, pw in enumerate(powers):
        ax1.annotate(str(pw), (pw, r2_vals[i]), textcoords="offset points",
                     xytext=(0, 8), ha='center', fontsize=7, color='grey')
    # Highlight chosen power
    chosen_idx = powers.index(6)
    ax1.plot(6, r2_vals[chosen_idx], 'o', color='blue', markersize=10, zorder=5)
    ax1.annotate(f'Power={6}', (6, r2_vals[chosen_idx]),
                 textcoords="offset points", xytext=(15, -10),
                 fontsize=9, fontweight='bold', color='blue',
                 arrowprops=dict(arrowstyle='->', color='blue'))
    sns.despine(ax=ax1)

    # Panel 2: Mean connectivity
    ax2.plot(powers, mean_k, 'o-', color='#1F77B4', linewidth=1.5, markersize=5)
    ax2.set_xlabel('Soft Threshold (Power)', fontweight='bold')
    ax2.set_ylabel('Mean Connectivity', fontweight='bold')
    ax2.set_title('Mean Connectivity', fontweight='bold')
    ax2.plot(6, mean_k[chosen_idx], 'o', color='blue', markersize=10, zorder=5)
    sns.despine(ax=ax2)

    plt.suptitle('WGCNA Soft Threshold Selection', fontsize=13, fontweight='bold', y=1.03)
    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}FigS_W7_SoftThreshold.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W7_SoftThreshold.pdf")

    # FigS-W8: Gene Dendrogram + Module Colors
    print("\n--- FigS-W8: Gene Dendrogram + Module Colors ---")

    from scipy.cluster.hierarchy import dendrogram

    # Use sampled genes for tractable dendrogram
    # Get all genes that are in modules (non-grey)
    mod_genes_all = gene_modules[gene_modules['moduleColors'] != 'grey']
    avail_genes = [g for g in mod_genes_all.index if g in df_expr.columns]

    # Sample if too many
    MAX_GENES = 3000
    if len(avail_genes) > MAX_GENES:
        np.random.seed(42)
        # Ensure focus module genes are included
        focus_genes = []
        for mod in FOCUS:
            fg = gene_modules[gene_modules['moduleColors'] == mod].index
            focus_genes.extend([g for g in fg if g in df_expr.columns])
        focus_genes = list(set(focus_genes))
        remaining = [g for g in avail_genes if g not in focus_genes]
        n_extra = MAX_GENES - len(focus_genes)
        if n_extra > 0:
            sampled = list(np.random.choice(remaining, min(n_extra, len(remaining)), replace=False))
        else:
            sampled = []
        avail_genes = focus_genes + sampled
    
    print(f"  Using {len(avail_genes)} genes for dendrogram")

    sub_expr = df_expr[avail_genes].values.T  # genes x samples
    # Correlation distance
    cor = np.corrcoef(sub_expr)
    cor = np.nan_to_num(cor, nan=0)
    dist = 1 - cor
    np.fill_diagonal(dist, 0)
    dist = np.clip(dist, 0, None)

    # Condensed distance
    from scipy.spatial.distance import squareform
    dist_condensed = squareform(dist, checks=False)
    Z = linkage(dist_condensed, method='average')

    # Module colors for these genes
    gene_colors = [gene_modules.loc[g, 'moduleColors'] if g in gene_modules.index else 'grey'
                   for g in avail_genes]

    # Convert module names to actual colors
    import matplotlib.colors as mcolors
    def mod_to_rgb(mod_name):
        if mod_name in FOCUS:
            return FOCUS[mod_name]['color']
        try:
            mcolors.to_rgba(mod_name)
            return mod_name
        except:
            return '#CCCCCC'

    color_array = [mod_to_rgb(c) for c in gene_colors]

    fig, (ax_dendro, ax_colors) = plt.subplots(
        2, 1, figsize=(14, 5),
        gridspec_kw={'height_ratios': [4, 0.4], 'hspace': 0.02})

    # Dendrogram
    dn = dendrogram(Z, ax=ax_dendro, no_labels=True, color_threshold=0,
                    above_threshold_color='#333333', leaf_rotation=0)
    ax_dendro.set_ylabel('1 - Pearson r', fontsize=9)
    ax_dendro.set_title('Gene Dendrogram and Module Colors', fontsize=13, fontweight='bold')
    ax_dendro.spines[['top', 'right', 'bottom']].set_visible(False)
    ax_dendro.tick_params(bottom=False)

    # Color bar below dendrogram (reordered by dendrogram leaves)
    leaf_order = dn['leaves']
    ordered_colors = [color_array[i] for i in leaf_order]

    for i, c in enumerate(ordered_colors):
        ax_colors.axvspan(i - 0.5, i + 0.5, color=c, linewidth=0)
    ax_colors.set_xlim(-0.5, len(ordered_colors) - 0.5)
    ax_colors.set_ylim(0, 1)
    ax_colors.set_yticks([0.5])
    ax_colors.set_yticklabels(['Module'], fontsize=8)
    ax_colors.tick_params(bottom=False, left=False)
    ax_colors.set_xticklabels([])
    ax_colors.spines[['top', 'right', 'bottom', 'left']].set_visible(False)

    plt.savefig(f"{FIG_DIR}FigS_W8_Dendrogram.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W8_Dendrogram.pdf")

    del df_expr, adata_mc_w, adata_mc

    # FigS-W9: Module Eigengene Clustering (ME Dendrogram)
    print("\n--- FigS-W9: ME Clustering Dendrogram ---")

    # Recompute ME means per module
    me_mean = me_df.groupby('celltype')[mod_cols].mean()
    # Correlation between modules (across cell types)
    me_mod_cor = me_mean[mod_cols].corr(method='pearson')
    me_dist = 1 - me_mod_cor
    np.fill_diagonal(me_dist.values, 0)
    me_dist_condensed = squareform(me_dist.values, checks=False)
    me_Z = linkage(me_dist_condensed, method='average')

    fig, ax = plt.subplots(figsize=(max(10, len(mod_cols) * 0.4), 5))
    dn_me = dendrogram(me_Z, labels=mod_cols, ax=ax, leaf_rotation=90,
                       color_threshold=0.2,  # ME_DISS_THRES = 0.2
                       above_threshold_color='#333333',
                       leaf_font_size=7)

    # Merge threshold line
    ax.axhline(0.2, color='red', ls='--', lw=1.2, alpha=0.7)
    ax.text(ax.get_xlim()[1] * 0.95, 0.22, 'Merge threshold\n(MEDissThres=0.2)',
            ha='right', fontsize=8, color='red', fontstyle='italic')

    ax.set_ylabel('1 - Pearson r', fontweight='bold')
    ax.set_title('Module Eigengene Dendrogram', fontsize=13, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    # Color the focus module labels
    xlabels = ax.get_xticklabels()
    for label in xlabels:
        mod = label.get_text()
        if mod in FOCUS:
            label.set_color(FOCUS[mod]['color'])
            label.set_fontweight('bold')
            label.set_fontsize(9)

    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}FigS_W9_ME_Dendrogram.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W9_ME_Dendrogram.pdf")

    # FigS-W10: Module Merge Comparison (Dynamic vs Merged colors)
    print("\n--- FigS-W10: Module Merge Comparison ---")

    # We only have final merged colors. Show module relationships via
    # a combined view: ME dendrogram + heatmap of inter-module correlations
    # This is the standard "eigengene adjacency" heatmap

    # Reorder by dendrogram
    leaf_order_me = dn_me['leaves']
    ordered_mods = [mod_cols[i] for i in leaf_order_me]
    me_adj = (1 + me_mod_cor.loc[ordered_mods, ordered_mods]) / 2  # adjacency

    fig, ax = plt.subplots(figsize=(max(8, len(ordered_mods) * 0.5),
                                    max(7, len(ordered_mods) * 0.45)))

    # Module color annotation
    mod_color_series = pd.Series(
        {m: FOCUS[m]['color'] if m in FOCUS else
         (m if m in mcolors.CSS4_COLORS else '#AAAAAA')
         for m in ordered_mods})

    # Draw heatmap
    mask = np.zeros_like(me_adj.values, dtype=bool)
    sns.heatmap(me_adj, cmap='RdYlBu_r', vmin=0, vmax=1,
                linewidths=0.3, linecolor='white',
                cbar_kws={'label': 'Eigengene Adjacency', 'shrink': 0.5},
                ax=ax, annot=False, square=True)

    ax.set_title('Module Eigengene Adjacency\n(ordered by clustering)',
                 fontsize=12, fontweight='bold')
    ax.set_xticklabels(ordered_mods, rotation=90, fontsize=6)
    ax.set_yticklabels(ordered_mods, rotation=0, fontsize=6)

    # Color focus module labels
    for label in ax.get_xticklabels():
        mod = label.get_text()
        if mod in FOCUS:
            label.set_color(FOCUS[mod]['color'])
            label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        mod = label.get_text()
        if mod in FOCUS:
            label.set_color(FOCUS[mod]['color'])
            label.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(f"{FIG_DIR}FigS_W10_ME_Adjacency.{FIG_FORMAT}", dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigS_W10_ME_Adjacency.pdf")

else:
    print("  SKIP: no metacell h5ad found. Need one of:")
    for mp in mc_path_candidates:
        print(f"    {mp}")
    print("  Save from notebook: adata_metacell.write(f'{OUTPUT_DIR}adata_metacell.h5ad')")

# DONE
print("\n" + "=" * 60)
print("ALL WGCNA FIGURES DONE!")
print("=" * 60)
print(f"Output: {FIG_DIR}")
print("""
  Main Figures:
    Fig1A_Transcriptional_Similarity.pdf  — ME correlation clustermap
    Fig1B_Pathway_DotPlot.pdf             — Hallmark+KEGG enrichment DotPlot

  Supplementary:
    FigS_W1_Module_Sizes.pdf              — module sizes
    FigS_W2_ME_Heatmap.pdf                — ME z-score heatmap
    FigS_W3_Focus_Violins.pdf             — focus module violin
    FigS_W4_Cohens_d.pdf                  — Cohen's d specificity heatmap
    FigS_W5_UMAP_ModuleScores.pdf         — UMAP module scores
    FigS_W6_HubGene_DotPlot.pdf           — hub gene DotPlot
    FigS_W7_SoftThreshold.pdf             — soft threshold curve
    FigS_W8_Dendrogram.pdf                — gene dendrogram + module colors
    FigS_W9_ME_Dendrogram.pdf             — ME dendrogram (with merge threshold)
    FigS_W10_ME_Adjacency.pdf             — ME adjacency heatmap
""")
