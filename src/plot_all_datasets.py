import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp as mc


def plot_a_protein(pool, single, protein_name, label_name = None, filename = None, figure_path = None):


    colors = {'pos': '#FEC200', 
         'neg': '#0C6F9F'}

    # Create a DF for plotting
    pool_data = pd.DataFrame([pool.X[:, pool.var['GeneName'] == protein_name].flatten(), pool.obs['TauStatus']]).T
    pool_data.columns = [protein_name, 'TauStatus']
    pool_data['Modality'] = '20x Cells'
    
    sc_data = pd.DataFrame([single.X[:, single.var['GeneName'] == protein_name].flatten(), single.obs['TauStatus']]).T
    sc_data.columns = [protein_name, 'TauStatus']
    sc_data['Modality'] = '1x Cell'

    # Concat DF
    all_data = pd.concat([pool_data, sc_data])
    all_data[protein_name] = all_data[protein_name].astype('float64')
    all_data['Combined'] = all_data['Modality'] + '_' + all_data['TauStatus']
    
    # Plot
    plt.figure(figsize=(3,3), constrained_layout=True)
    sns.violinplot(x=all_data['Modality'], y=all_data[protein_name], density_norm='width', 
                   hue=all_data['TauStatus'], hue_order=['negative', 'positive'], palette=[colors['neg'], colors['pos']])
    sns.despine()
    if label_name:
        plt.ylabel(f'{label_name} (Log2 Intensity)')
    else:
        plt.ylabel(f'{protein_name} (Log2 Intensity)')
    plt.xlabel('')
    plt.legend(title='pTau')
    plt.savefig(os.path.join(figure_path, f'{filename}.pdf'))
    plt.savefig(os.path.join(figure_path, f'{filename}.png'), dpi=300, bbox_inches='tight')
    plt.show()

    # Stats
    model = ols(formula=f'{protein_name} ~ TauStatus * Modality ', data=all_data).fit()
    aov_table = sm.stats.anova_lm(model, typ=2)
    print(aov_table)

    comp = mc.MultiComparison(all_data[protein_name], all_data['Combined'])
    post_hoc_res = comp.tukeyhsd()
    print(post_hoc_res.summary())


if __name__ == '__main__':
    
    # Load data
    processed_data_path = os.path.join('..', 'data', 'processed_data')
    figure_path = os.path.join('..', 'figures', 'figure1')
    
    pool = ad.read_h5ad(os.path.join(processed_data_path, 'pool_processed.h5ad'))
    single = ad.read_h5ad(os.path.join(processed_data_path, 'sc_processed.h5ad'))

    # Plot all datasets
    plot_a_protein(pool, single, protein_name='MAPT', label_name='Tau', filename='Fig1e_MAPT_levels', figure_path=figure_path)