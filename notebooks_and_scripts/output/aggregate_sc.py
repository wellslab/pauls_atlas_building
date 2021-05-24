#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt


# In[2]:


from general_processing.SC_functions import aggregate_clusters
from general_processing.process_functions import convert_symbols_to_ensembl


# In[3]:


data = sc.read_h5ad('/Users/pwangel/Data/Single_Cell/Gut_Colon/Colon_cell_atlas_scrublet.h5ad')


# In[4]:


test1, test2 = aggregate_clusters(data, 'donor', n_samples_per_cluster=100, average_library_size_per_cluster=1.e6)


# In[5]:


test1


# In[6]:


test1 = convert_symbols_to_ensembl(test1)


# In[7]:


test1


# In[ ]:




