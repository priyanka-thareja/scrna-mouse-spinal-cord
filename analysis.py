"""
Single-cell RNA-seq analysis of mouse spinal cord
Dataset: Sathyamurthy et al. (2018) - GSE103892
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, facecolor='white')


# Load data
counts = pd.read_csv('data/counts.txt', sep='\t', index_col=0)
meta = pd.read_csv('data/metadata.txt', sep='\t', index_col=0)

adata = sc.AnnData(counts.T)  # transpose: scanpy wants cells as rows
adata.obs = meta.loc[adata.obs_names]


# QC filtering
adata = adata[adata.obs['cell.type'] != 'discarded'].copy()

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
plt.tight_layout()
plt.savefig('figures/01_qc_violin.png', dpi=150)
plt.close()

adata = adata[adata.obs['n_genes_by_counts'] > 200]   # remove empty droplets
adata = adata[adata.obs['n_genes_by_counts'] < 5000]  # remove doublets
adata = adata[adata.obs['pct_counts_mt'] < 20]        # remove dying cells
sc.pp.filter_genes(adata, min_cells=10)


# Normalization
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig('figures/02_hvg.png', dpi=150)
plt.close()

adata = adata[:, adata.var['highly_variable']]
sc.pp.scale(adata, max_value=10)


# Dimensionality reduction and clustering
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
plt.savefig('figures/03_pca_elbow.png', dpi=150)
plt.close()

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

sc.pl.umap(adata, color='leiden', show=False)
plt.savefig('figures/04_umap_clusters.png', dpi=150, bbox_inches='tight')
plt.close()


# Marker genes
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=5, show=False)
plt.savefig('figures/05_marker_genes.png', dpi=150, bbox_inches='tight')
plt.close()

marker_df = sc.get.rank_genes_groups_df(adata, group=None)
marker_df.to_csv('results/marker_genes.csv', index=False)


# Cell type annotation
broad_types = adata.obs['cell.type'].astype(str).copy()
broad_types[broad_types.str.match('^[DVM][EI]-|^VC-|^VE-|^ME-|^MI-|^neurons')] = 'Neurons'
adata.obs['broad_type'] = broad_types

sc.pl.umap(adata, color='broad_type', show=False)
plt.savefig('figures/06_umap_celltypes.png', dpi=150, bbox_inches='tight')
plt.close()


# Opioid receptor expression
opioid_genes = ['Oprm1', 'Oprd1', 'Oprk1']

sc.pl.umap(adata, color=opioid_genes, use_raw=True, show=False, ncols=3)
plt.savefig('figures/07_opioid_receptors.png', dpi=150, bbox_inches='tight')
plt.close()

sc.pl.dotplot(adata, var_names=opioid_genes, groupby='broad_type', use_raw=True, show=False)
plt.savefig('figures/08_opioid_dotplot.png', dpi=150, bbox_inches='tight')
plt.close()


# Save
adata.write('data/sathyamurthy_final.h5ad')
