import numpy as np
import pandas as pd
import scipy, sys, sklearn.decomposition
import statsmodels.api as sm
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

gene_list      = ['KIT',
                  'JAM2',
                  'VCAM1',
                  'CD226',
                  'SELP'
                  ]


data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_annotations.tsv', sep='\t', index_col=0)
genes_s4m      = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes_v7.1.tsv', sep='\t', index_col=0)
genes_varPart  = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t', index_col=0)
genes = genes_s4m.merge(genes_varPart['Platform_VarFraction'], how='left', left_index=True, right_index=True)

cut_data    = functions.transform_to_percentile(data.loc[genes.loc[genes.inclusion.values].index.values.astype(str), annotations.index.values])
all_ranked_data = functions.transform_to_percentile(data.loc[:, annotations.index.values])

fname      = 'macrophages'
folder     = 'myeloid_pluripotent_macrophages'

sample_types = []
sample_ids   = []

# Create a list of sample types and the ids for each sample type
for i_sample_type in [['in vitro', 'myeloid', 'macrophage'], ['in vitro', 'pluripotent stem cell', 'macrophage']]:

    sample_types.append(i_sample_type[0]+'/'+i_sample_type[1]+'/'+i_sample_type[2])
    sel_samples =  (annotations.tier1==i_sample_type[0]) & (annotations.tier2==i_sample_type[1]) & (annotations.celltype==i_sample_type[2])

    sample_ids.append(annotations.loc[sel_samples].index.values)
    
    print(i_sample_type, "N samples = %d" %annotations.loc[sel_samples].shape[0])
    print("Sample platforms = ", np.unique(annotations.loc[sel_samples].Platform_Category.values))
    print("N per respective platform = ", np.unique(annotations.loc[sel_samples].Platform_Category.values, return_counts=True)[1])

annotations = annotations.loc[np.concatenate(sample_ids)]
cut_data    = cut_data[np.concatenate(sample_ids)]
cut_unfiltered_data = all_ranked_data[np.concatenate(sample_ids)]

pvals        = np.array([])
delta_median = np.array([])

df_output = pd.DataFrame(index=gene_list, columns=np.append(np.append(sample_types, [i+' Mean' for i in sample_types]), [i+' Std' for i in sample_types]))

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    for i_sample_type, i_sample_ids in zip(sample_types, sample_ids):

        ensembl_id = genes.index.values[genes.symbol.values==i_gene]

        if ensembl_id in cut_data.index.values:
                
            ranks_of_sample  = cut_data.loc[ensembl_id, i_sample_ids].values
            ranks_not_sample = cut_data.loc[ensembl_id, ~np.in1d(cut_data.columns.values, i_sample_ids)].values
    
            pvals = np.append(pvals, scipy.stats.mannwhitneyu(ranks_of_sample.flatten(), ranks_not_sample.flatten(), alternative='two-sided')[1])
            delta_median = np.append(delta_median, np.abs(np.median(ranks_of_sample)-np.median(ranks_not_sample)))

        else:

            ranks_of_sample  = cut_unfiltered_data.loc[ensembl_id, i_sample_ids].values
            ranks_not_sample = cut_unfiltered_data.loc[ensembl_id, ~np.in1d(cut_unfiltered_data.columns.values, i_sample_ids)].values
    
            pvals = np.append(pvals, scipy.stats.mannwhitneyu(ranks_of_sample.flatten(), ranks_not_sample.flatten(), alternative='two-sided')[1])
            delta_median = np.append(delta_median, np.abs(np.median(ranks_of_sample)-np.median(ranks_not_sample)))

        df_output.loc[i_gene, i_sample_type+' Mean'] = ranks_of_sample.mean()
        df_output.loc[i_gene, i_sample_type+' Std'] = ranks_of_sample.std()


#Correct for multiple testing
corrected_pvals = sm.stats.fdrcorrection(pvals)
corrected_pvals = sm.stats.multipletests(pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)

#Loop through again and place into dataframe
i_pval    = 0

for i_gene in gene_list:

    for i_sample_type, i_sample_ids in zip(sample_types, sample_ids):

        ensembl_id = genes.index.values[genes.symbol.values==i_gene]

        df_output.loc[i_gene, i_sample_type] = corrected_pvals[1][i_pval]

        i_pval += 1                

df_output.loc[~df_output.iloc[:,0].isnull()].to_csv('/Users/pwangel/PlotlyWorkspace/iMac_plots/%s/%s_pvals.tsv' %(folder, fname), sep='\t')

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



