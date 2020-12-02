import numpy as np
import pandas as pd
import scipy, sys, sklearn.mixture
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

gene_list      = [
                  'KLF4',
                  'KLF2',
                  'STAT3',
                  'TCF3',
                  'TBX3',
                  'GBX2',
                  'SALL4',
                  'ESRRB',
                  'TFCP2L1',
                  'POU5F1',
                  'NANOG',
                  'SOX2'
                  #'MEKERK'
                  ]


data           = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)
genes_s4m      = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t', index_col=0)
genes_conversion  = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])
genes = genes_s4m.merge(genes_conversion, how='left', left_index=True, right_index=True)

cut_data    = functions.transform_to_percentile(data.loc[genes.loc[genes.inclusion.values].index.values.astype(str), annotations.index.values])
all_ranked_data = functions.transform_to_percentile(data.loc[:, annotations.index.values])

sel_samples = np.ones(shape=annotations.shape[0]).astype(bool)

annotations = annotations.loc[sel_samples]
cut_data    = cut_data.loc[:,sel_samples]
cut_unfiltered_data = all_ranked_data.loc[:,sel_samples]

gmm = sklearn.mixture.GaussianMixture(n_components=2)

df_output = pd.DataFrame(index=gene_list, columns=['Mean 1', 'Mean 2', 'Cov 1', 'Cov 2'])

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    ensembl_id = genes.index.values[genes.symbol.values==i_gene]

    gmm.fit(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))

    df_output.loc[i_gene, 'Mean 1'] = gmm.means_[0,0]
    df_output.loc[i_gene, 'Mean 2'] = gmm.means_[1,0]
    df_output.loc[i_gene, 'Cov 1'] = gmm.covariances_[0,0,0]
    df_output.loc[i_gene, 'Cov 2'] = gmm.covariances_[1,0,0]

