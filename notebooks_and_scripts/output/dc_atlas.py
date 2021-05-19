#!/usr/bin/env python
# coding: utf-8

# This is a notebook to include new dc samples in the blood and myeloid atlases.

# In[1]:


import numpy as np
import pandas as pd
import sklearn, json
import functions
from general_processing.process_functions import convert_symbols_to_ensembl, remove_microarray_duplicates


# Reading in blood data first. Annoyingly have to switch the format of the sample ids

# In[2]:


data           = pd.read_csv('/Users/pwangel/dc_atlas/DCAtlasFiles/DCAtlasFiles_v1.3/expression.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/dc_atlas/DCAtlasFiles/DCAtlasFiles_v1.3/samples.tsv', sep='\t', index_col=0)
data = data[annotations.index]
annotations['Dataset'] = [i.split('_')[0] for i in annotations.index.values]

colours_df = pd.read_json('/Users/pwangel/dc_atlas/DCAtlasFiles/DCAtlasFiles_v1.3/colours.json')
colour_dict = {}
[colour_dict.update(i_dict) for i_dict in colours_df.colours.values]

print(data.shape, annotations.shape)


# In[98]:


ext_data1       = pd.read_csv('/Users/pwangel/Data/External_Data/GSE136731/GSE136731_expression.tsv', sep='\t', index_col=0)
ext_annotations1 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE136731/GSE136731_annotation.tsv', sep='\t', index_col=0)
ext_annotations1['Dataset'] = 'GSE136731'

ext_data2       = pd.read_csv('/Users/pwangel/Data/External_Data/GSE87494/GSE87494_expression.tsv', sep='\t', index_col=0)
ext_annotations2 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE87494/GSE87494_annotation.tsv', sep='\t', index_col=0)
ext_annotations2['Dataset'] = 'GSE87494'

ext_data3       = pd.read_csv('/Users/pwangel/Data/External_Data/GSE144435/GSE144435_expression.tsv', sep='\t', index_col=0)
ext_annotations3 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE144435/GSE144435_annotation.tsv', sep='\t', index_col=0)
ext_annotations3['Platform_Category'] = 'RNAseq'
ext_annotations3['Dataset'] = 'GSE144435'

ext_data4       = pd.read_csv('/Users/pwangel/Data/External_Data/GSE151086/GSE151086_expression.tsv', sep='\t', index_col=0)
ext_annotations4 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE151086/GSE151086_annotation.tsv', sep='\t', index_col=0)
ext_annotations4['Platform_Category'] = 'RNAseq'
ext_annotations4['Dataset'] = 'GSE151086'

ext_data5       = pd.read_csv('/Users/pwangel/Data/External_Data/GSE151073/GSE151073_expression.tsv', sep='\t', index_col=0)
ext_annotations5 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE151073/GSE151073_annotation.tsv', sep='\t', index_col=0)
ext_annotations5['Platform_Category'] = 'RNAseq'
ext_annotations5['Dataset'] = 'GSE151073'

ext_data6 = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/Microarray_Data/datasets/dataset7002.gct', sep='\t', index_col=0, skiprows=2)
ext_data6 = ext_data6.drop(['Description'], axis=1)
ext_data6_map= pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/Microarray_Data/probe_mappings/mapping_16.txt', delim_whitespace=True, header=None)
ext_data6_map.columns = ['NAME', 'Gene']
ext_data6_map  = ext_data6_map.set_index('NAME')
ext_data6_map = ext_data6_map.loc[np.intersect1d(ext_data6_map.index.values, ext_data6.index.values)]
ext_data6 = remove_microarray_duplicates(ext_data6, ext_data6_map)
ext_data6.columns = ['7002_'+str(i) for i in ext_data6.columns]
ext_annotations6 = pd.read_csv('/Users/pwangel/Downloads/7002.txt', sep='\t', index_col=0)
ext_annotations6['Dataset'] = 7002
ext_data6 = ext_data6[ext_annotations6.index]

ext_data = pd.DataFrame(index=data.index)
for i_ext in [ext_data1, ext_data2, ext_data3, ext_data4, ext_data5, ext_data6]:
    ext_data = ext_data.merge(i_ext, how='left', left_index=True, right_index=True)
ext_annotations = pd.concat([ext_annotations1, ext_annotations2, ext_annotations3, ext_annotations4, ext_annotations5, ext_annotations6])

print(ext_data.shape, ext_annotations.shape)


# In[4]:


#data = data.merge(ext_data, how='inner', left_index=True, right_index=True)
#annotations = pd.concat([annotations, ext_annotations])


# In[5]:


weird_index = annotations.loc[(annotations['Platform Category']=='Illumina V4')&(annotations['Sample Source']=='in vivo')].index
annotations.loc[weird_index, 'Platform Category'] = 'Illumina V4 2'


# In[7]:


data = functions.transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already, in the above we have already read a previously calculated version into the gene dataframe

# In[8]:


annotations.rename(columns={'Platform Category':'Platform_Category'}, inplace=True)
genes = functions.calculate_platform_dependence(data, annotations)


# In[9]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=1.0]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=1.0]).transpose())


# In[47]:


annotations['display_metadata'] = annotations.index
functions.plot_pca(pca_coords, annotations,pca,                    labels=['Cell Type', 'Sample Source', 'Platform_Category', 'Dataset'], colour_dict=colour_dict)#, out_file='/Users/pwangel/Downloads/dc_atlas.html')


# In[11]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.15]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.15]).transpose())


# Plot the pca

# In[12]:


annotations['display_metadata'] = annotations.index
functions.plot_pca(pca_coords, annotations,pca,                    labels=['Cell Type', 'Sample Source', 'Platform_Category', 'Dataset'], colour_dict=colour_dict)#, out_file='/Users/pwangel/Downloads/dc_atlas.html')


# In[99]:


ext_data = ext_data.loc[genes.loc[genes.Platform_VarFraction.values<=0.15].index]
ext_annotations['N Zeroes'] = ext_data.isnull().sum().values
ext_data.fillna(0.0)
annotations['N Zeroes'] = np.nan
ext_coords = pca.transform(functions.transform_to_percentile(ext_data).transpose())
#ext_annotations.rename(columns={'Platform Category':'Platform_Category'}, inplace=True)


# In[100]:


ext_annotations['Sample Source'] = ext_annotations['Sample Source'].str.replace('_', '')
ext_annotations['Platform_Category'] = ext_annotations['Platform_Category'].str.replace('_', '')


# In[102]:


functions.plot_pca([pca_coords, ext_coords], [annotations, ext_annotations],pca,                    labels=['Cell Type', 'Sample Source', 'Platform_Category', 'Dataset', 'N Zeroes'], colour_dict=colour_dict)#, out_file='/Users/pwangel/Downloads/dc_atlas.html')


# In[ ]:




