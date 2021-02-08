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
                  'SOX2',
                  'TFE3', 
                  'KLF17', 
                  #'SSEA3', why can't I find this? 
                  #'TRA-1-60', why can't I find this?
                  'CD24', 
                  #'SSEA4', why can't I find this? 
                  'NLGN4X', 
                  'F11R', 
                  'EPCAM', 
                  'OTX2', 
                  'ZIC2', 
                  'ZIC3', 
                  'ZIC5', 
                  'KLF5', 
                  'DPPA3', 
                  'MAEL',
                  'UTF1', 
                  'ZFP57',
                  'DNMT3L'
                  #'MEKERK'
                  ]


#data = pd.DataFrame(index=main_ensembl_ids)
data = pd.DataFrame()
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='outer', left_index=True, right_index=True)

#data        = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/RNASeq_pluripotent_annotations.tsv', sep='\t', index_col=0)
genes       = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])

annotations.index = [i+';'+j for i,j in zip(annotations.chip_id.values.astype(str), annotations.Dataset.values.astype(int).astype(str))]
#data        = data[annotations.index]
data        = data[annotations.chip_id]

data = functions.transform_to_percentile(data)

gmm = sklearn.mixture.GaussianMixture(n_components=2)

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    ensembl_id = genes.index.values[genes.symbol.values==i_gene]

    if i_gene == 'POU5F1':
        ensembl_id = ['ENSG00000204531']
    if i_gene == 'ZFP57':
        ensembl_id = ['ENSG00000204644']

    gmm.fit(data.loc[ensembl_id].values.reshape(-1, 1))
    prediction = gmm.predict(data.loc[ensembl_id].values.reshape(-1, 1))
    n_ambiguous_samples = ((gmm.predict_proba(data.loc[ensembl_id].values.reshape(-1, 1))[:,0] >= 0.25) & \
                                         (gmm.predict_proba(data.loc[ensembl_id].values.reshape(-1, 1))[:,0] < 0.75)).sum()

    #delta_mean = data.loc[ensembl_id,prediction==0].values.mean()-data.loc[ensembl_id,prediction==1].values.mean()
    #std_sum    = data.loc[ensembl_id,prediction==0].values.std()+data.loc[ensembl_id,prediction==1].values.std()
    delta_mean = data.loc[ensembl_id,annotations.LM_Group_COLOR.values=='naive'].values.mean()-data.loc[ensembl_id,annotations.LM_Group_COLOR.values=='primed'].values.mean()
    std_sum    = data.loc[ensembl_id,annotations.LM_Group_COLOR.values=='naive'].values.std()+data.loc[ensembl_id,annotations.LM_Group_COLOR.values=='primed'].values.std()

    print(i_gene, ensembl_id, n_ambiguous_samples, np.abs(delta_mean)/std_sum)

    fig = Figure()
    fig.add_trace(Histogram(x=data.loc[ensembl_id, prediction==0].values[0], name='Group 0',    
            xbins=dict(
                       start=0.0,
                       end=1.0,
                       size=0.025)
            ))

    fig.add_trace(Histogram(x=data.loc[ensembl_id, prediction==1].values[0], name='Group 1',    
            xbins=dict(
                       start=0.0,
                       end=1.0,
                       size=0.025)
            ))
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.5)
    #plot(fig)