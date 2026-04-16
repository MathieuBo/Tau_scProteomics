# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell and mini-bulk (pooled ~20 cells) proteomics of laser-capture microdissected (LCM) neurons from human Alzheimer's disease frontal cortex. Neurons are stratified by phospho-tau (pTau) status (positive vs. negative) using immunofluorescence prior to LCM. DIA-PASEF on timsTOF Ultra.

## Data Architecture

- **Raw data**: Perseus-exported tab-separated matrices (`data/raw_data/perseus_filtered_imputed_*.txt`) with multi-row headers (file paths, sample IDs, batch, tau_status, patient ID, modality). Intensities are log2-transformed.
- **Metadata**: `metadata_ptau.csv` (per-injection, pTau epitope quantification), `scLCMdemographics.xlsx` (patient demographics), `AT8_K48_raw_data_excel.xlsx` (immunofluorescence quantification).
- **Ubiquitin data**: `ubi_data_new.csv` — K48-ubiquitinated peptide intensities (precursor-level, UniMod:121).
- **Processed data**: AnnData `.h5ad` files in `data/processed_data/`. Observations = cells/pools, variables = proteins (indexed by Uniprot ID, gene names in `var['GeneName']`). `all_prot_bk.pkl` is a large (~2.4 GB) pickle of breakpoint analysis results.
- Protein intensities are **log2-scaled** throughout. Back-transform with `2 ** intensity` when needed.

## Analysis Pipeline

Notebooks run sequentially:
- `0-Preprocessing` — AnnData assembly from Perseus exports, contaminant removal, demographic merge, ubiquitin/pTau merge, MAPT QC filter, writes `.h5ad`
- `1-Figure1` — protein detection, dynamic range, neuronal markers, EWCE, tau peptides
- `1-2-TauProfile` — smoothed tau peptide detection profile across 2N4R sequence
- `2-Figure2_minipools` — pseudotime (Palantir), correlation with pseudotime, pre-ranked GSEA, hierarchical clustering heatmaps
- `3-Figure2_singlecells` — tau-level correlation, volcano, clustered heatmap, cell death pathway heatmap
- `4-Figure3_proteostasis` — proteasome/V-ATPase subunit correlation matrices, pathway scores, K48Ub-tau correlation, IF validation
- `5-Figure4_breakpoint` — piecewise regression per protein, trajectory classification (early/late, up/down)

Shared modules in `src/`: `plot_all_datasets.py` (cross-modality violin + ANOVA), `pathway_score.py` (z-score pathway scoring), `breakpoint_analysis.py` (robust z-score + piecewise regression).

## Key Conventions

- **AnnData is the central data structure** — not raw DataFrames.
- **Tau status**: `obs['TauStatus']` — values `'positive'`/`'negative'` (strings).
- **Color palette**: Sapphire/Dandelion — `pos='#FEC200'` (gold), `neg='#0C6F9F'` (blue).
- **Two modalities**: "20x Cells" (pooled) and "1x Cell" (single-cell).
- **Statistics**: Pearson correlation + BH FDR for genome-wide screens; two-way ANOVA + Tukey HSD for cross-modality comparisons; piecewise regression (BIC model selection) for breakpoint analysis.
- **Figures**: saved as both PNG (300 dpi) and PDF, organised in `figures/` subdirectories by panel.
- **Importing src from notebooks**: `sys.path.append('..')` then `from src.module import func`.

## Environment Setup

```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[notebooks]"
```
