import numpy as np
import pandas as pd 
import sys, os, scipy.stats, sklearn.decomposition

def convert_symbols_to_ensembl(dataframe):

    '''
    Convert expression matrix from gene symbol indexed to ensembl id index

    Parameters:
    ----------

    dataframe
        Dataframe containing expression values, index as gene symbols

    Returns:
    -----------

    dataframe
        Expression dataframe with variables converted to ensembl ids

    '''

    conversion = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', header=None, names=['ensembl', 'symbol'])
    conversion.drop_duplicates(subset='symbol', keep=False, inplace=True)

    dataframe  = dataframe.merge(conversion, how='left', left_index=True, right_on='symbol')

    dataframe.set_index(['ensembl'], inplace=True)
    dataframe.drop(columns=['symbol'], inplace=True)
    dataframe = dataframe.loc[~dataframe.index.isnull()]

    return dataframe


# Generate an atlas using pluripotent stem cell data
# Annotations expected to be updated in the near future

# Read in and find cut
data = pandas.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pandas.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)

# Temporary hack to include som externally processed data and see how the atlas would look
df1 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE114873/GSE114873_Cufflinks_output_FPKM_and_Gene_Names.txt', sep='\t', index_col=0)[['condition', 'replicate', 'raw_frags']]
df1['sample_id'] = [i+'_'+j for i,j in zip(df1.condition.to_list(),df1.replicate.to_list())]
df1.drop(labels=['condition', 'replicate'], axis=1, inplace=True)
df1=df1.pivot(columns='sample_id')
df1.columns = df1.columns.droplevel()
df1 = convert_symbols_to_ensembl(df1)

df2 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE123055/GSE123055_counts.txt', sep='\t', index_col=0)
df2.drop(labels=['Unnamed: 52'], axis=1, inplace=True)

ext_annotations = pd.DataFrame(index=df1.columns.to_list()+df2.columns.to_list())
ext_annotations['Platform_Category']='RNASeq'
ext_annotations['Dataset'] = ['GSE114873' for i in range(df1.shape[1])] +  ['GSE123055' for i in range(df2.shape[1])]
ext_annotations['sample_type'] = 'Unannotated'
ext_annotations['display_metadata']  = [str(i_dataset)+'<br>'+str(i_sample) for i_dataset, i_sample in zip(ext_annotations.Dataset.values, ext_annotations.index.values)]

datasets_to_keep = pd.read_csv('/Users/pwangel/Downloads/Pluripotency datasets in stemformatics - existing stemformatics data.tsv', sep='\t')['Dataset'].values
annotations = annotations.loc[np.in1d(annotations.Dataset, datasets_to_keep)]

annotations = pd.concat([annotations, ext_annotations])
data = data.merge(df1, how='inner', left_index=True, right_index=True)
data = data.merge(df2, how='inner', left_index=True, right_index=True)
data = data[annotations.index]

data = functions.transform_to_percentile(data)
genes = functions.calculate_platform_dependence(data, annotations)

pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())

functions.plot_pca(pca_coords, annotations,pca, \
                   labels=['sample_type', 'Platform_Category', 'Dataset']+list(nadias_annotations.keys()), colour_dict=blood_atlas_colours)
