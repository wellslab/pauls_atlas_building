#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy


# In[2]:


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


# In[3]:


data        = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)

annotations['display_metadata'] = annotations.generic_sample_type
annotations = annotations.loc[annotations.Platform_Category!='Illumina V2']


# Temporary hack to include some externally processed data and see how the atlas would look

# In[4]:


df1 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE114873/GSE114873_Cufflinks_output_FPKM_and_Gene_Names.txt', sep='\t', index_col=0)[['condition', 'replicate', 'raw_frags']]
df1['sample_id'] = [i+'_'+j for i,j in zip(df1.condition.to_list(),df1.replicate.to_list())]
df1.drop(labels=['condition', 'replicate'], axis=1, inplace=True)
df1=df1.pivot(columns='sample_id')
df1.columns = df1.columns.droplevel()
df1 = convert_symbols_to_ensembl(df1)

df2 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE123055/GSE123055_counts.txt', sep='\t', index_col=0)
df2.drop(labels=['Unnamed: 52'], axis=1, inplace=True)

df3 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE131551/GSE131551_human_bulk.raw_count_matrix.tsv', sep='\t', index_col=0)
df3.drop(labels='geneID', axis=1, inplace=True)


# In[5]:


ext_annotations = pd.DataFrame(index=df1.columns.to_list()+df2.columns.to_list()+df3.columns.to_list())
ext_annotations['Platform_Category']='RNASeq'
ext_annotations['Dataset'] = ['GSE114873' for i in range(df1.shape[1])] +  ['GSE123055' for i in range(df2.shape[1])] + ['GSE131551' for i in range(df3.shape[1])]
ext_annotations['generic_sample_type'] = 'Unannotated'
ext_annotations['display_metadata']  = [str(i_dataset)+'<br>'+str(i_sample) for i_dataset, i_sample in zip(ext_annotations.Dataset.values, ext_annotations.index.values)]


# In[6]:


datasets_to_keep = pd.read_csv('/Users/pwangel/Downloads/Pluripotency datasets in stemformatics - existing stemformatics data.tsv', sep='\t')['Dataset'].values
annotations = annotations.loc[np.in1d(annotations.Dataset, datasets_to_keep)]

annotations = pd.concat([annotations, ext_annotations])
data = data.merge(df1, how='inner', left_index=True, right_index=True)
data = data.merge(df2, how='inner', left_index=True, right_index=True)
data = data.merge(df3, how='inner', left_index=True, right_index=True)
data = data[annotations.index]


# In[7]:


data = functions.transform_to_percentile(data)


# In[8]:


#genes = functions.calculate_platform_dependence(data, annotations)
#genes.to_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes_with_ext.tsv', sep='\t')
#genes = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t')
genes = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes_with_ext.tsv', sep='\t') 


# In[9]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())


# In[10]:


pca_coords_ext = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25][ext_annotations.index]).transpose())


# In[11]:


functions.plot_pca([pca_coords, pca_coords_ext], [annotations, ext_annotations],pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset'], colour_dict={}, pcs=[1,2,3], out_file='/Users/pwangel/Downloads/pluripotent_atlas_with_external.html')


# Now try to 'zoom in' on the pluripotent cells (isolate them by applying k means clustering)

# In[12]:


kmeans = sklearn.cluster.KMeans(n_clusters=4).fit(pca_coords)
annotations['K Means'] = kmeans.labels_
ext_annotations['K Means'] = annotations['K Means'].loc[ext_annotations.index]


# In[13]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset', 'K Means'], colour_dict={}, pcs=[1,2,3], out_file='/Users/pwangel/Downloads/pluripotent_kmeans_atlas_with_external.html')


# In[14]:


pluripotent_annotations = annotations.loc[annotations['K Means']==0]
pluripotent_data = data[pluripotent_annotations.index]
#pluripotent_genes = functions.calculate_platform_dependence(pluripotent_data, pluripotent_annotations)
#pluripotent_genes.to_csv('/Users/pwangel/Downloads/pluripotent_only_atlas_genes_with_ext.tsv', sep='\t')
pluripotent_genes = pd.read_csv('/Users/pwangel/Downloads/pluripotent_only_atlas_genes_with_ext.tsv', sep='\t')


# In[15]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(pluripotent_data.loc[pluripotent_genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(pluripotent_data.loc[pluripotent_genes.Platform_VarFraction.values<=0.2]).transpose())


# In[16]:


functions.plot_pca(pca_coords, pluripotent_annotations,pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset'], colour_dict={}, pcs=[1,2,3], out_file='/Users/pwangel/Downloads/pluripotent_only_atlas_with_external.html')


# In[ ]:




