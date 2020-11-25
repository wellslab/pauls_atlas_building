#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy


# In[13]:


data        = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)
genes       = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t', index_col=0)

annotations['display_metadata'] = annotations.generic_sample_type


# In[3]:


data = functions.transform_to_percentile(data)


# In[4]:


#genes = functions.calculate_platform_dependence(data, annotations)
genes = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_genes.tsv', sep='\t') 


# In[5]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())


# In[14]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset'], colour_dict={})


# In[8]:


annotations.columns


# In[ ]:




