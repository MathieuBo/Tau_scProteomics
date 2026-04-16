import os
import numpy as np
from scipy.stats import zscore
import scanpy as sc

def calc_pathway_score(adata, pathway_genes):
    """
    Calculate the pathway score for a given pathway
    """
    data_path = adata[:, adata.var['GeneName'].isin(pathway_genes)].copy()
    pathway_zs = zscore(data_path.X.T, axis=1) # Median/MAD zscoring
    pathway_score = np.mean(pathway_zs, axis=0)
    return pathway_score


if __name__ == "__main__":
    
    processed_data_path = '../data/processed_data'
    adata = sc.read_h5ad(os.path.join(processed_data_path, 'sc_processed.h5ad'))
    pathway_genes = ['ATP6V0A1', 'ATP6V0A2', 'ATP6V0A3', 'ATP6V0A4', 'ATP6V0A5', 'ATP6V0A6', 'ATP6V0A7', 'ATP6V0A8', 'ATP6V0A9', 'ATP6V0A10']
    pathway_score = calc_pathway_score(adata, pathway_genes)
    print(pathway_score)