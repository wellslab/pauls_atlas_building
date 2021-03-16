#!/usr/bin/env python
# coding: utf-8

# This is an example notebook walking through the construction of the atlas

# In[40]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy


# In[41]:


blood_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/imac_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}


# Reading in data, including nadias annotations, excel spreadsheet with multiple tabs

# In[50]:


data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_annotations.tsv', sep='\t', index_col=0)
genes          = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t', index_col=0)

ext_data       = pd.read_csv('/Users/pwangel/Downloads/DC_expression_matrix.txt', sep='\t', index_col=0)
ext_annotations = pd.read_csv('/Users/pwangel/Downloads/DC_samples_matrix.txt', sep='\t', index_col=0) 


# In[51]:


print(ext_annotations.columns)
print(annotations.columns)
ext_annotations.rename(columns={"dataset_id":"Dataset", "platform_category":"Platform_Category", "cell_type":"celltype"}, inplace=True)


# In[53]:


ext_annotations.index = [str(i)+';'+str(j) for i,j in zip(ext_annotations.external_source_id.values,ext_annotations.Dataset.values.astype(int))]
ext_data.columns = ext_annotations.index
ext_annotations = ext_annotations.loc[np.setdiff1d(ext_annotations.index, annotations.index)]
ext_data = ext_data.loc[:,ext_annotations.index]
print(np.intersect1d(annotations.index.values, ext_annotations.index.values).shape)


# In[54]:


data = data.merge(ext_data, how='inner', left_index=True, right_index=True)
annotations = pd.concat([annotations, ext_annotations])


# In[55]:


annotations.shape


# In[56]:


data = functions.transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already, in the above we have already read a previously calculated version into the gene dataframe

# In[57]:


genes = functions.calculate_platform_dependence(data, annotations)
genes.to_csv('/Users/pwangel/Downloads/temp_ext_myeloid_atlas_genes.tsv', sep='\t') 


# In[58]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())


# Plot the pca

# In[59]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['celltype', 'Platform_Category', 'Dataset'], colour_dict=blood_atlas_colours)


# In[7]:


functions.plot_gene_platform_dependence_distribution(data, annotations, genes)


# Make a graph of the threshold lowering process using the Kruskal Wallis H Test

# In[ ]:


functions.plot_KW_Htest(data, annotations, genes)

