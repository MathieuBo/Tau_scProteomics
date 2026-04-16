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

### Python

| Package | Purpose |
|---------|---------|
| anndata | Data structure (AnnData / .h5ad) |
| scanpy | Single-cell analysis toolkit |
| pandas | Data manipulation |
| numpy | Numerical computing |
| scipy | Statistics (correlation, FDR, clustering, interpolation) |
| statsmodels | ANOVA, post-hoc tests |
| matplotlib | Plotting |
| seaborn | Statistical visualisation |
| PyComplexHeatmap | Annotated clustered heatmaps |
| adjustText | Non-overlapping text labels |
| palantir | Pseudotime inference (diffusion maps) |
| gseapy | Gene set enrichment analysis (pre-ranked GSEA, Enrichr) |
| scikit-learn | Dimensionality reduction |
| piecewise-regression | Breakpoint / segmented regression |
| joblib | Parallel computation |
| tqdm | Progress bars |
| openpyxl | Excel file reading |

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
