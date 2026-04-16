import numpy as np
from scipy.stats import median_abs_deviation
import piecewise_regression

def robust_zscore(data, axis=None):
    """
    Compute robust z-scores using median and MAD.
    
    Parameters:
    -----------
    data : array-like
        Input data
    axis : int or None
        Axis along which to compute (None for entire array)
    
    Returns:
    --------
    z_scores : array
        Robust z-scores
    """
    median = np.median(data, axis=axis, keepdims=True)
    mad = np.median(np.abs(data - median), axis=axis, keepdims=True)

    # Scale factor for consistency with standard normal distribution
    mad_scaled = 1.4826 * mad

    # Avoid division by zero
    mad_scaled = np.where(mad_scaled == 0, 1, mad_scaled)

    z_robust = (data - median) / mad_scaled

    return z_robust


def bkpoint_analysis(p, adata):

    res = {}
    res['GeneName'] = p

    # Subset data to pathway proteins
    data_path = adata[:, p].copy()

    # Calculate pathway score
    pathway_zs = robust_zscore(data_path.X.T, axis=1) # Median/MAD zscoring

    # Winsorization
    med = np.median(data_path.X)
    mad = median_abs_deviation(data_path.X, scale='normal')
    limit_lower = med - (3 * mad)
    limit_upper = med + (3 * mad)

    mask_keep = (data_path.X >= limit_lower) & (data_path.X <= limit_upper)

    data_path = data_path[mask_keep, :]

    data_path.obs[p] = pathway_zs.reshape(-1,1)[mask_keep]

    # Sort data
    expression_df_sorted = data_path.obs.sort_values('TauLevels').reset_index(drop=True)
    # Run model selection
    ms = piecewise_regression.ModelSelection(expression_df_sorted['TauLevels'].values, expression_df_sorted[p].values, max_breakpoints=2, n_boot=1000)

    # Choose best model
    bics = []
    best_bic = np.inf
    best_n_breakpoints = None

    for mod in ms.model_summaries:

        bics.append(mod['bic'])
        if mod['bic'] is not None:
            if mod['bic'] < best_bic:
                best_bic = mod['bic']
                best_n_breakpoints = mod['n_breakpoints']
    
    res['n_breakpoint'] = best_n_breakpoints

    if best_n_breakpoints > 0:
        res['model'] = ms.models[best_n_breakpoints - 1]
    
    return res

