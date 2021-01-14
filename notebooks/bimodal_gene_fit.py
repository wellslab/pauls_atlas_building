import numpy as np
import pandas as pd
import scipy, sys, sklearn.mixture
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

main_ensembl_ids = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/ensembl_hg38.91_chromosome.csv').ensembl_gene_id.values.astype(str)

dataset_list = ['/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7124/source/processed.hg38.91.20180713-201511/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7135/source/processed.hg38.91.20180716-154848/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7240/source/processed.hg38.91.20180718-124924/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7253/source/processed.hg38.91.20180815-171954/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6884/source/processed.hg38.91.20180627-112849/gene_count_frags.txt']

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


data = pd.DataFrame(index=main_ensembl_ids)
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='left', left_index=True, right_index=True)

#data           = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)
lizzis_anno    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcell_annotations.tsv', sep='\t')
annotations = annotations.merge(lizzis_anno[['chip_id', 'LM_Group_COLOR']], how='inner', left_on='chip_id', right_on='chip_id')
genes_s4m      = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t', index_col=0)
genes_conversion  = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])
genes_conversion  = genes_conversion.loc[main_ensembl_ids]
genes = genes_s4m.merge(genes_conversion, how='left', left_index=True, right_index=True)

annotations = annotations.loc[np.in1d(annotations.Dataset.values.astype(int), [7124, 7135, 7240, 6884, 7253])]
annotations = annotations.loc[np.in1d(annotations.LM_Group_COLOR, ['naive', 'primed'])]
data = data[annotations.chip_id] #Not sure if the samples are in the right order

# Add some single cell data

sc_data     = pd.read_csv('/Users/pwangel/Data/Single_Cell/Han/aggregated_by_cluster_100_0pt0.tsv', sep='\t', index_col=0)
sc_annotations = pd.read_csv('/Users/pwangel/Data/Single_Cell/Han/aggregated_by_cluster_metadata_100_0pt0.tsv', sep='\t', index_col=0)
sc_annotations['LM_Group_COLOR'] = sc_annotations.celltype.values

data = data.merge(sc_data, how='inner', left_index=True, right_index=True).fillna(0.0)
annotations = pd.concat([annotations, sc_annotations])

# Choose normalisation or transformation
#data = np.log2(1.e6*data/data.sum()+1)
data = functions.transform_to_percentile(data)

cut_data    = data
all_ranked_data = data

#cut_data    = functions.transform_to_percentile(data.loc[genes.loc[genes.inclusion.values].index.values.astype(str), annotations.index.values])
#all_ranked_data = functions.transform_to_percentile(data.loc[:, annotations.index.values])

sel_samples = np.ones(shape=annotations.shape[0]).astype(bool)

annotations = annotations.loc[sel_samples]
cut_data    = cut_data.loc[:,sel_samples]
cut_unfiltered_data = all_ranked_data.loc[:,sel_samples]

gmm = sklearn.mixture.GaussianMixture(n_components=2)

df_output = pd.DataFrame(index=gene_list, columns=['Mean 1', 'Mean 2', 'Cov 1', 'Cov 2', 'Ambiguous Samples', 'Naive', 'Primed'])

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    ensembl_id = genes_conversion.index.values[genes_conversion.symbol.values==i_gene]

    gmm.fit(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))

    df_output.loc[i_gene, 'Mean 1'] = gmm.means_[0,0]
    df_output.loc[i_gene, 'Mean 2'] = gmm.means_[1,0]
    df_output.loc[i_gene, 'Cov 1'] = gmm.covariances_[0,0,0]
    df_output.loc[i_gene, 'Cov 2'] = gmm.covariances_[1,0,0]
    df_output.loc[i_gene, 'Ambiguous Samples'] = ((gmm.predict_proba(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))[:,0] >= 0.25) & \
                                         (gmm.predict_proba(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))[:,0] < 0.75)).sum()

    df_output.loc[i_gene, 'Naive'] = gmm.predict_proba(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))[annotations.LM_Group_COLOR=='naive',0].mean()
    df_output.loc[i_gene, 'Primed'] = gmm.predict_proba(cut_unfiltered_data.loc[ensembl_id].values.reshape(-1, 1))[annotations.LM_Group_COLOR=='primed',0].mean()

data = data.loc[data.sum(axis=1)>10,:]
pca        = sklearn.decomposition.PCA(n_components=3)
print("PCA variance:")
output = pca.fit_transform(data.transpose())
#output = pca.fit_transform((functions.transform_to_percentile(data).div(functions.transform_to_percentile(data).std(axis=1), axis=0)).transpose())
#output = pca.fit_transform((functions.transform_to_percentile(data)).transpose())
print(pca.explained_variance_ratio_)

data_to_plot = []
for i_type in np.unique(annotations.LM_Group_COLOR.values.astype(str)):

    sel  = (annotations.LM_Group_COLOR==i_type)
    hover = np.core.defchararray.add(annotations.generic_sample_type.values[sel].astype(str), np.full(shape=annotations.loc[sel].shape[0], fill_value='<br>'))
    hover = np.core.defchararray.add(hover, annotations.Handle.values[sel].astype(str))
    data_to_plot.append(Scatter3d(x=output[sel,0], y=output[sel,1], z=output[sel,2], opacity=0.9, mode='markers', name=str(i_type), text=hover, marker=dict(size=5)))

fig = Figure(data=data_to_plot, 
        layout=Layout(title="Pluripotency PCA", xaxis=dict(title='Component 1'), yaxis=dict(title='Component 2')))
plot(fig, auto_open=False, filename='/users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/naive_sc_comp1_2_3.html')

genes_conversion = genes_conversion.loc[data.index]
genes_conversion['PC2 loading'] = pca.components_[1,:]
print(genes_conversion.loc[np.in1d(genes_conversion.symbol, gene_list),:])

print(stop)
annotations['Platform_Category'] = annotations.Dataset.values.astype(int).astype(str)
annotations.index=data.columns.values
data = data.loc[data.sum(axis=1)>10,:]
dataset_varFrac = functions.calculate_platform_dependence(data.loc[data.sum(axis=1)>10,:], annotations)
dataset_varFrac = dataset_varFrac.merge(genes_conversion, how='inner', left_index=True, right_index=True)

output = pca.fit_transform(functions.transform_to_percentile(data.loc[dataset_varFrac.Platform_VarFraction.values < 0.15]).transpose())
print(pca.explained_variance_ratio_)

data_to_plot = []
for i_type in np.unique(annotations.LM_Group_COLOR.values.astype(str)):

    sel  = (annotations.LM_Group_COLOR==i_type)
    hover = np.core.defchararray.add(annotations.generic_sample_type.values[sel].astype(str), np.full(shape=annotations.loc[sel].shape[0], fill_value='<br>'))
    hover = np.core.defchararray.add(hover, annotations.Handle.values[sel].astype(str))
    data_to_plot.append(Scatter3d(x=output[sel,0], y=output[sel,1], z=output[sel,2], opacity=0.9, mode='markers', name=str(i_type), text=hover, marker=dict(size=5)))

fig = Figure(data=data_to_plot, 
        layout=Layout(title="Pluripotency PCA", xaxis=dict(title='Component 1'), yaxis=dict(title='Component 2')))
plot(fig, auto_open=False, filename='/users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/cut_naive_sc_comp1_2_3.html')