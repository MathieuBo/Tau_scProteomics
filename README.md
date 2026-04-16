![Python 3.12](https://img.shields.io/badge/python-3.12-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

# Single-cell proteomics of tau-bearing neurons in Alzheimer's disease (Foiani, Bourdenx, Kraller, et al. 2026)

This repository contains the analysis code for the single-cell and mini-bulk (pooled ~20 cells) proteomics study of laser-capture microdissected (LCM) neurons from human Alzheimer's disease frontal cortex. Neurons were stratified by phospho-tau (pTau) immunofluorescence status prior to microdissection and analysed by DIA-PASEF mass spectrometry on a timsTOF Ultra.


## Data availability

**Raw data and results are embargoed until publication** and are therefore excluded from this repository via `.gitignore`. Processed AnnData objects and raw intensity matrices will be deposited in a public repository upon publication.

## Repository structure

```
├── notebooks/                   # Analysis notebooks (run in order)
│   ├── 0-Preprocessing.ipynb    # AnnData assembly, QC, filtering
│   ├── 1-Figure1.ipynb          # Protein detection, dynamic range, neuronal markers
│   ├── 1-2-TauProfile.ipynb     # Tau peptide detection profile across 2N4R sequence
│   ├── 2-Figure2_minipools.ipynb# Mini-bulk: pseudotime, correlation, GSEA, heatmaps
│   ├── 3-Figure2_singlecells.ipynb # Single-cell: tau correlation, clustering, cell death
│   ├── 4-Figure3_proteostasis.ipynb # Proteasome/lysosome scores, K48-ubiquitin, V-ATPase
│   └── 5-Figure4_breakpoint.ipynb   # Piecewise regression / breakpoint analysis
├── src/                         # Shared Python modules and R scripts
│   ├── plot_all_datasets.py     # Cross-modality protein plotting with ANOVA + Tukey HSD
│   ├── pathway_score.py         # Z-score-based pathway scoring
│   ├── breakpoint_analysis.py   # Robust z-score and piecewise regression per protein
│   └── EWCE.R                   # Cell-type enrichment analysis (requires R, see below)
├── data/                        # Raw and processed data (not tracked)
│   ├── raw_data/                # Perseus exports, metadata, demographics, ubiquitin
│   └── processed_data/          # .h5ad files, EWCE results, breakpoint pickle
├── results/                     # Supplementary tables and gene lists (not tracked)
└── figures/                     # Output figures organised by panel
```

## Dependencies

Tested with Python 3.12.

### Python

| Package | Version | Purpose |
|---------|---------|---------|
| anndata | 0.11.4 | Data structure (AnnData / .h5ad) |
| scanpy | 1.11.1 | Single-cell analysis toolkit |
| pandas | 2.2.3 | Data manipulation |
| numpy | 2.2.5 | Numerical computing |
| scipy | 1.15.2 | Statistics (correlation, FDR, clustering, interpolation) |
| statsmodels | 0.14.4 | ANOVA, post-hoc tests |
| matplotlib | 3.10.1 | Plotting |
| seaborn | 0.13.2 | Statistical visualisation |
| PyComplexHeatmap | 1.8.2 | Annotated clustered heatmaps |
| adjustText | 1.3.0 | Non-overlapping text labels |
| palantir | 1.4.1 | Pseudotime inference (diffusion maps) |
| gseapy | 1.1.11 | Gene set enrichment analysis (pre-ranked GSEA, Enrichr) |
| scikit-learn | 1.5.2 | Dimensionality reduction |
| piecewise-regression | 1.5.0 | Breakpoint / segmented regression |
| joblib | 1.4.2 | Parallel computation |
| tqdm | 4.67.1 | Progress bars |
| openpyxl | 3.1.5 | Excel file reading |

### R

`src/EWCE.R` performs cell-type enrichment analysis and requires R with the following packages: [EWCE](https://github.com/NathanSkene/EWCE), MAGMA.Celltyping, reticulate, SingleCellExperiment, scKirby.

## Setup

Requires Python >= 3.11. Use [uv](https://docs.astral.sh/uv/) for environment setup:

```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[notebooks]"
```

## Running the analysis

Start a Jupyter session from the project root:

```bash
jupyter lab
```
Then open notebooks in `notebooks/` and run them sequentially. The preprocessing notebook (`0-Preprocessing.ipynb`) generates the `.h5ad` files consumed by all downstream notebooks.


## Author
Mathieu Bourdenx — UK Dementia Research Institute at UCL
