import numpy as np
import pandas as pd
import scipy, sys, sklearn.mixture
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

#Example script using the dot plot function 

# Reading in data

dataset_list = ['/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7124/source/processed.hg38.91.20180713-201511/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7135/source/processed.hg38.91.20180716-154848/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7240/source/processed.hg38.91.20180718-124924/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7253/source/processed.hg38.91.20180815-171954/gene_count_frags.txt']

# These are the genes to plot
gene_list      = [
'ENSG00000228592',
'ENSG00000211772',
'ENSG00000080503',
'ENSG00000180730',
'ENSG00000274021',
'ENSG00000180730',
'ENSG00000105696',
'ENSG00000196542',
'ENSG00000211751',
'ENSG00000117154'
                  ]


data = pd.DataFrame()
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='outer', left_index=True, right_index=True)

#data           = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)
lizzis_anno    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcell_annotations.tsv', sep='\t')
annotations = annotations.merge(lizzis_anno[['chip_id', 'LM_Group_COLOR']], how='inner', left_on='chip_id', right_on='chip_id')
#genes_s4m      = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t', index_col=0)
genes_conversion  = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])
#genes_conversion  = genes_conversion.loc[main_ensembl_ids]
#genes = genes_s4m.merge(genes_conversion, how='left', left_index=True, right_index=True)
#data = data.merge(genes_conversion['symbol'], how='outer', left_index=True, right_index=True).set_index('symbol')

annotations = annotations.loc[np.in1d(annotations.Dataset.values.astype(int), [7124, 7135, 7240, 7253])]
annotations = annotations.loc[np.in1d(annotations.LM_Group_COLOR, ['naive', 'primed'])]
data = data[annotations.chip_id] #Not sure if the samples are in the right order

# Using CPM transformation
data = np.log2(1e6*data/data.sum()+1) 

# Dot plot function defined in functions.py
functions.plot_dot_plots(data, annotations, cell_property='LM_Group_COLOR', cell_colour_by='Dataset', gene_list=gene_list, output_dir = '/users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/')
