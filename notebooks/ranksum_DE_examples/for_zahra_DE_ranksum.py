import numpy as np
import pandas as pd
import scipy, sys, sklearn.decomposition
import statsmodels.api as sm
import functions

#### This is an example script utilising the Mann Whitney Ranksum test implement in scipy.
#### The groups being tested here are the in vitro vs in vivo DC1 cells

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

data         = pd.read_csv('/location/of/expression.tsv', sep='\t', index_col=0)
annotations  = pd.read_csv('/location/of/annotations.tsv', sep='\t', index_col=0)
genes        = pd.read_csv('/location/of/myeloid_atlas_genes_v7.1.tsv', sep='\t', index_col=0)

cut_data    = functions.transform_to_percentile(data.loc[genes.inclusion.values])

# Select only DC1 samples

annotations = annotations.loc[annotations.tier1=='DC1'] #### Select only DC1 cells for example
cut_data    = cut_data[annotations.index.values]

pvals        = np.array([])
delta_median = np.array([])

# Define dataframe to keep results in. Also keep the mean and std of each group for the hell of it

df_output = pd.DataFrame(index=gene_list, columns= ['P val', 'In vitro mean', 'In vivo mean', 'In vitro std', 'In vivo std'])

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in cut_data.index.values:

    in_vitro_sample = cut_data.loc[i_gene, annotations.tier2=='In vitro'].values
    in_vivo_sample  = cut_data.loc[i_gene, annotations.tier2=='In vivo'].values
   
    # Perform the test here  
 
    pvals = np.append(pvals, scipy.stats.mannwhitneyu(in_vitro_sample.flatten(), in_vivo_sample.flatten(), alternative='two-sided')[1])
    #delta_median = np.append(delta_median, np.abs(np.median(in_vitro_sample)-np.median(in_vivo_sample)))

    df_output.loc[i_gene, 'In vitro mean'] = in_vitro_sample.mean()
    df_output.loc[i_gene, 'In vitro std'] = in_vitro_sample.std()

    df_output.loc[i_gene, 'In vivo mean'] = in_vivo_sample.mean()
    df_output.loc[i_gene, 'In vivo std'] = in_vivo_sample.std()

# Include the gene symbol names for ease of viewing

df_output = df_output.merge(genes['symbol'], how='left', left_index=True, right_index=True)

#Correct for multiple testing

#corrected_pvals = sm.stats.fdrcorrection(pvals)
corrected_pvals = sm.stats.multipletests(pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)

df_output['P val'] = corrected_pvals
df_output.to_csv('/location/to/save.tsv', sep='\t')


# The following is a template for making and saving the dot plots for each gene. Paul used to do this (takes a while though).
'''
import seaborn as sns
import matplotlib.pyplot as pyplot
sns.set(style="ticks")

for i_gene in gene_list:

    ensembl_id = genes.index.values[genes.symbol.values==i_gene]
    df_violin  = pd.DataFrame(columns=['Expression', 'Cell Property', 'Cell Colour'])
    platform_std = np.sqrt((genes.Platform_VarFraction.values*all_ranked_data.std(axis=1).values)[genes.symbol.values==i_gene][0])
    platform_mean = all_ranked_data.loc[ensembl_id].mean(axis=1)[0]

    for i_sample_type, i_sample_ids in zip(sample_types, sample_ids):

        if genes.loc[genes.symbol.values==i_gene].inclusion.values:    
            i_type_vals = cut_data.loc[ensembl_id, i_sample_ids].values.flatten().astype(float)
        else:
            i_type_vals = all_ranked_data.loc[ensembl_id, i_sample_ids].values.flatten().astype(float)

        tier1_vals    = np.full(shape=i_type_vals.shape[0], fill_value=i_sample_type)
        df_violin     = df_violin.append(pd.DataFrame(columns=['Expression', 'Cell Property'], data = np.array([i_type_vals, tier1_vals]).transpose()), ignore_index=True)

    df_violin.Expression = df_violin.Expression.astype(float)


    fig, ax = pyplot.subplots(1,1, figsize=(13.5,7.0))
    ax = sns.swarmplot(x="Cell Property", y="Expression", size=9, order = sample_types, color='blue', data=df_violin)

    ax.fill_between(x=[ax.get_xlim()[0], ax.get_xlim()[1]], y1=platform_mean-platform_std, y2=platform_mean+platform_std, color='lightgrey', alpha=0.5)

    pyplot.ylabel(i_gene+' Expression')
    pyplot.ylim(0.0, 1.0)
    pyplot.show()
    #pyplot.savefig('/Users/pwangel/PlotlyWorkspace/iMac_plots/%s/%s/%s.eps' %(folder, fname,i_gene))
    ax.remove()
    pyplot.close()
'''


