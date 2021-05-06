#!/usr/bin/env python
# coding: utf-8

# This is a notebook to include new dc samples in the blood and myeloid atlases.

# In[21]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy


# In[22]:


ensembl_to_symbol = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])


# In[23]:


blood_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/blood_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}


# Reading in blood data first. Annoyingly have to switch the format of the sample ids

# In[24]:


data           = pd.read_csv('/Users/pwangel/Downloads/blood_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/blood_atlas_samples_v7.1.tsv', sep='\t', index_col=0)
annotations = annotations.merge(    pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_annotations.tsv', sep='\t', index_col=0)['Platform_Category'],    how='left', left_index=True, right_index=True)
data.columns = [i.split(';')[1]+'_'+i.split(';')[0] for i in data.columns.values.astype(str)]
annotations.index = [i.split(';')[1]+'_'+i.split(';')[0] for i in annotations.index.values.astype(str)]
data = data[annotations.index]

print(data.shape, annotations.shape)

ext_data       = pd.read_csv('/Users/pwangel/Downloads/dc_atlas_expression_v1.3.tsv', sep='\t', index_col=0)

ext_annotations = pd.read_csv('/Users/pwangel/Downloads/dc_bloodatlas_nadia_moresamples.txt', sep='\t', index_col=0)
ext_annotations.index = [i.split(';')[1]+'_'+i.split(';')[0] for i in ext_annotations.index.values.astype(str)]
ext_annotations = ext_annotations.loc[[i.split('_')[0]!='3378' for i in ext_annotations.index.astype(str)]] #Apparently dataset 3378 is duplicated from 6612 and does not belong

print(ext_data.shape, ext_annotations.shape)

ext_annotations = ext_annotations.merge(    pd.read_csv('/Users/pwangel/Downloads/dc_atlas_samples_v1.3.tsv', sep='\t', index_col=0)['Platform Category'],     how='inner', left_index=True, right_index=True)
ext_data = ext_data[ext_annotations.index]

print(ext_data.shape, ext_annotations.shape)
print(np.intersect1d(ext_annotations.index, annotations.index))


# In[25]:


print(ext_annotations.columns)
print(annotations.columns)
ext_annotations.rename(columns={"Platform Category":"Platform_Category"}, inplace=True)


# In[26]:


data = data.merge(ext_data, how='inner', left_index=True, right_index=True)
annotations = pd.concat([annotations, ext_annotations])


# In[27]:


data = functions.transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already, in the above we have already read a previously calculated version into the gene dataframe

# In[28]:


genes = functions.calculate_platform_dependence(data, annotations)
genes['inclusion'] = (genes.Platform_VarFraction <=0.2)
genes = genes.merge(ensembl_to_symbol, how='left', left_index=True, right_index=True)
genes.index.name='ensembl'
genes.to_csv('/Users/pwangel/Downloads/blood_atlas_genes_v2.tsv', sep='\t') 
#genes = pd.read_csv('/Users/pwangel/Downloads/blood_atlas_genes_v2.tsv', sep='\t', index_col=0)
annotations.to_csv('/Users/pwangel/Downloads/blood_atlas_samples_v2.tsv', sep='\t')
data.to_csv('/Users/pwangel/Downloads/blood_atlas_expression_v2.tsv', sep='\t')
data.loc[genes.loc[genes.inclusion].index].to_csv('/Users/pwangel/Downloads/blood_atlas_expression_v2.filtered.tsv', sep='\t')


# In[29]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pd.DataFrame(data=pca_coords, index=annotations.index, columns = ['PCA'+str(i) for i in range(1,11)]).to_csv('/Users/pwangel/Downloads/blood_atlas_coordinates_v2.tsv', sep='\t')


# Plot the pca

# In[30]:


annotations['display_metadata'] = annotations.index
functions.plot_pca(pca_coords, annotations,pca,                    labels=['Cell Type', 'Sample Source', 'Progenitor Type', 'Platform_Category'], colour_dict=blood_atlas_colours, out_file='/Users/pwangel/Downloads/blood_atlas_with_ext_dc.html')


# In[31]:


myeloid_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/imac_atlas_colours.tsv', sep='\t').set_index('Sample Source')
myeloid_atlas_colours = {key:value[0] for key, value in zip(myeloid_atlas_colours.index.values, myeloid_atlas_colours.values)}


# In[32]:


data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_samples_v7.1 (3).tsv', sep='\t', index_col=0)
data.columns = [i.split(';')[1]+'_'+i.split(';')[0] for i in data.columns.values.astype(str)]
annotations.index = [i.split(';')[1]+'_'+i.split(';')[0] for i in annotations.index.values.astype(str)]
data = data[annotations.index]

print(data.shape, annotations.shape)

ext_data       = pd.read_csv('/Users/pwangel/Downloads/dc_atlas_expression_v1.3.tsv', sep='\t', index_col=0)

ext_annotations = pd.read_csv('/Users/pwangel/Downloads/dc_myeloidatlas_nadia.txt', sep='\t', index_col=0)
ext_annotations.index = [i.split(';')[1]+'_'+i.split(';')[0] for i in ext_annotations.index.values.astype(str)]
ext_annotations = ext_annotations.loc[[i.split('_')[0]!='3378' for i in ext_annotations.index.astype(str)]] #Apparently dataset 3378 is duplicated from 6612 and does not belong

print(ext_data.shape, ext_annotations.shape)

ext_annotations = ext_annotations.merge(    pd.read_csv('/Users/pwangel/Downloads/dc_atlas_samples_v1.3.tsv', sep='\t', index_col=0)['Platform Category'],     how='inner', left_index=True, right_index=True)
ext_data = ext_data[ext_annotations.index]

print(ext_data.shape, ext_annotations.shape)
print(np.intersect1d(ext_annotations.index, annotations.index))


# In[33]:


print(ext_annotations.columns)
print(annotations.columns)
annotations.rename(columns={"Platform Category":"Platform_Category"}, inplace=True)


# In[34]:


data = data.merge(ext_data, how='inner', left_index=True, right_index=True)
annotations = pd.concat([annotations, ext_annotations])


# In[35]:


data = functions.transform_to_percentile(data)


# In[36]:


genes = functions.calculate_platform_dependence(data, annotations)
genes['inclusion'] = (genes.Platform_VarFraction <=0.2)
genes = genes.merge(ensembl_to_symbol, how='left', left_index=True, right_index=True)
genes.index.name='ensembl'
genes.to_csv('/Users/pwangel/Downloads/myeloid_atlas_genes_v2.tsv', sep='\t') 
#genes = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes_v2.tsv', sep='\t', index_col=0)
annotations.to_csv('/Users/pwangel/Downloads/myeloid_atlas_samples_v2.tsv', sep='\t')
data.to_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v2.tsv', sep='\t')
data.loc[genes.loc[genes.inclusion].index].to_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v2.filtered.tsv', sep='\t')


# In[37]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pd.DataFrame(data=pca_coords, index=annotations.index, columns = ['PCA'+str(i) for i in range(1,11)]).to_csv('/Users/pwangel/Downloads/myeloid_atlas_coordinates_v2.tsv', sep='\t')


# In[38]:


annotations['display_metadata'] = annotations.index
functions.plot_pca(pca_coords, annotations,pca,                    labels=['Cell Type', 'Sample Source', 'Progenitor Type', 'Platform_Category'], colour_dict=blood_atlas_colours, out_file='/Users/pwangel/Downloads/myeloid_atlas_with_ext_dc.html')

