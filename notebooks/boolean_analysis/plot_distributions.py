import numpy as np
import pandas as pd
import scipy, sys, sklearn.mixture
sys.path.append('/Users/pwangel/pauls_atlas_building/notebooks/')
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

main_ensembl_ids = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/ensembl_hg38.91_chromosome.csv').ensembl_gene_id.values.astype(str)

dataset_list = ['/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7124/source/processed.hg38.91.20180713-201511/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7135/source/processed.hg38.91.20180716-154848/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7240/source/processed.hg38.91.20180718-124924/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7253/source/processed.hg38.91.20180815-171954/gene_count_frags.txt']
                #'/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6884/source/processed.hg38.91.20180627-112849/gene_count_frags.txt']

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
                  'POU5F1', #multiple ensemble ids for this one..
                  'NANOG',
                  'SOX2'
                  #'MEKERK'
                  ]


#data = pd.DataFrame(index=main_ensembl_ids)
data = pd.DataFrame()
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='outer', left_index=True, right_index=True)

#data           = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)
lizzis_anno = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcell_annotations.tsv', sep='\t')
annotations = annotations.merge(lizzis_anno[['chip_id', 'LM_Group_COLOR']], how='inner', left_on='chip_id', right_on='chip_id')
genes_s4m      = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t', index_col=0)
genes_conversion  = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])
genes_conversion  = genes_conversion.loc[main_ensembl_ids]
genes = genes_s4m.merge(genes_conversion, how='left', left_index=True, right_index=True)

annotations = annotations.loc[np.in1d(annotations.Dataset.values.astype(int), [7124, 7135, 7240, 7253])]#, 6884, 7253])]
annotations = annotations.loc[np.in1d(annotations.LM_Group_COLOR, ['naive', 'primed'])]
data = data[annotations.chip_id] #Not sure if the samples are in the right order

data = functions.transform_to_percentile(data)

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    ensembl_id = genes_conversion.index.values[genes_conversion.symbol.values==i_gene]

    fig = Figure()
    for i_dataset in annotations.Dataset.unique():
        for i_type, i_colour in zip(['naive', 'primed'], ['red', 'blue']):

            sel  = (annotations.LM_Group_COLOR==i_type) & (annotations.Dataset==i_dataset)
            fig.add_trace(Histogram(x=data.loc[ensembl_id, sel.values].values[0], name=str(i_dataset)+' '+str(i_type), marker_color=i_colour,    
            xbins=dict(
                       start=0.0,
                       end=1.0,
                       size=0.05)
            ))

    #fig.update_layout(barmode='overlay')
    fig.update_layout(barmode='stack')
    fig.update_traces(opacity=0.5)
    plot(fig, auto_open=False, filename='/users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcells_gene_distributions_%s.html' %i_gene)
