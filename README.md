# Single-cell proteomics of tau-bearing neurons in Alzheimer's disease

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
├── src/                         # Shared Python modules
│   ├── plot_all_datasets.py     # Cross-modality protein plotting with ANOVA + Tukey HSD
│   ├── pathway_score.py         # Z-score-based pathway scoring
│   └── breakpoint_analysis.py   # Robust z-score and piecewise regression per protein
├── data/                        # Raw and processed data (not tracked)
│   ├── raw_data/                # Perseus exports, metadata, demographics, ubiquitin
│   └── processed_data/          # .h5ad files, EWCE results, breakpoint pickle
├── results/                     # Supplementary tables and gene lists (not tracked)
└── figures/                     # Output figures organised by panel (not tracked)
```

## Setup

Requires Python >= 3.11. Use [uv](https://docs.astral.sh/uv/) for environment setup:

```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[notebooks]"
```

## Running the analysis

Notebooks are numbered and should be run sequentially. The preprocessing notebook (`0-Preprocessing.ipynb`) generates the `.h5ad` files consumed by all downstream notebooks.

From within a notebook, shared utilities are imported via:

```python
import sys
sys.path.append('..')
from src.pathway_score import calc_pathway_score
```

## Authors

Mathieu Bourdenx — UKRI Future Leaders Fellow, UK Dementia Research Institute at UCL
