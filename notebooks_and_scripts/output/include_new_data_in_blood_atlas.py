#!/usr/bin/env python
# coding: utf-8

# This is an example notebook walking through the construction of the atlas

# In[2]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy
from general_processing.processing_functions import convert_symbols_to_ensembl, transform_to_percentile 


# In[2]:


blood_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/blood_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}


# Reading in data, including nadias annotations, excel spreadsheet with multiple tabs

# In[3]:


data           = pd.read_csv('/Users/pwangel/Downloads/blood_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/Downloads/blood_atlas_samples_v7.1.tsv', sep='\t', index_col=0)
#genes          = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t', index_col=0)

ext_data       = pd.read_csv('/Users/pwangel/Downloads/dc_atlas_expression.tsv', sep='\t', index_col=0)
ext_annotations = pd.read_csv('/Users/pwangel/Downloads/dc_bloodatlas_nadia.txt', sep='\t', index_col=0) 
ext_platform_anno = pd.read_csv('/Users/pwangel/Downloads/dc_myeloidatlas_nadia.txt', sep='\t', index_col=0) 

print(ext_annotations.shape)
ext_annotations = ext_annotations.merge(ext_platform_anno['Platform_Category'], how='inner', left_index=True, right_index=True)
print(ext_annotations.shape)


# In[13]:


annotations.columns


# In[23]:


i_anno = annotations.sample(frac=1.0, replace=True)
i_data = data[i_anno.index]
temp = i_data.copy().transpose()
temp['blah'] = i_anno['Cell Type']


# In[24]:


temp


# In[5]:


data = data[annotations.index]
annotations.index = [i.split(";")[1]+"_"+i.split(";")[0] for i in annotations.index]
data.columns = annotations.index
ext_annotations.index = [i.split(";")[1]+"_"+i.split(";")[0] for i in ext_annotations.index]
print(np.intersect1d(data.columns, ext_data.columns).shape)


# In[7]:


ext_data = ext_data.loc[:,ext_annotations.index]
print(ext_data.shape, ext_annotations.shape)
print(np.intersect1d(annotations.index.values, ext_annotations.index.values).shape)
ext_annotations.to_csv('/Users/pwangel/Downloads/new_dc_blood.tsv', sep='\t')


# In[47]:


ext_annotations


# In[48]:


data = data.merge(ext_data, how='inner', left_index=True, right_index=True)
annotations = pd.concat([annotations, ext_annotations])


# In[49]:


print(annotations.shape)
print(data.shape)


# In[50]:


data = transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already, in the above we have already read a previously calculated version into the gene dataframe

# In[51]:


#genes = functions.calculate_platform_dependence(data, annotations)
#genes.to_csv('/Users/pwangel/Downloads/temp_ext_blood_atlas_genes.tsv', sep='\t') 
genes = pd.read_csv('/Users/pwangel/Downloads/temp_ext_blood_atlas_genes.tsv', sep='\t', index_col=0) 


# In[52]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())


# Plot the pca

# In[54]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['celltype', 'Platform_Category', 'Dataset'], colour_dict=blood_atlas_colours, out_file='/Users/pwangel/Downloads/blood_atlas_with_ext_dc.html')


# In[7]:


functions.plot_gene_platform_dependence_distribution(data, annotations, genes)


# Make a graph of the threshold lowering process using the Kruskal Wallis H Test

# In[ ]:


functions.plot_KW_Htest(data, annotations, genes)

