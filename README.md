# Single-Cell RNA-seq Analysis of Mouse Spinal Cord

Analysis of single-nucleus RNA-seq data from mouse spinal cord to characterize cell type diversity and opioid receptor expression patterns.

## Dataset

Sathyamurthy et al. (2018) "Massively Parallel Single Nucleus Transcriptional Profiling Defines Spinal Cord Neurons and Their Activity during Behavior" - Cell Reports

- Source: [GSE103892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103892)
- 18,000 cells, 32,962 genes
- After QC: 13,396 cells, 21,501 genes

## Methods

- QC filtering based on gene counts (200-5000) and mitochondrial content (<20%)
- Normalization to 10k counts per cell + log transform
- 2000 highly variable genes selected for clustering
- PCA (15 components) → Leiden clustering (resolution 0.5) → UMAP
- Marker genes identified using Wilcoxon rank-sum test

## Key Findings

- Identified 13 clusters corresponding to major spinal cord cell types: neurons, astrocytes, oligodendrocytes, microglia, vascular cells, and meninges/Schwann cells
- Opioid receptor expression analysis revealed Oprm1 (mu receptor) is enriched in VM-1 and VM-2 ventral motor neuron populations
- Clusters align well with author annotations

## Figures

| Figure | Description |
|--------|-------------|
| 01_qc_violin.png | QC metrics distribution |
| 02_hvg.png | Highly variable gene selection |
| 03_pca_elbow.png | PCA variance explained |
| 04_umap_clusters.png | UMAP colored by Leiden clusters |
| 05_marker_genes.png | Top marker genes per cluster |
| 06_umap_celltypes.png | UMAP colored by cell type |
| 07_opioid_receptors.png | Opioid receptor expression on UMAP |
| 08_opioid_dotplot.png | Opioid receptor expression by cell type |

## Usage
```bash
conda create -n scrna python=3.10
conda activate scrna
pip install scanpy leidenalg matplotlib pandas numpy

python analysis.py
```

## Author

Priyanka Thareja
