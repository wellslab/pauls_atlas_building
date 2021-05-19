#!/usr/bin/env python
# coding: utf-8

# This is a notebook used to analysis constclust results. The majority of constclust is run on hpc.

# In[4]:


from constclust import cluster, reconcile, plotting

import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt

from bokeh.io import show, output_notebook
output_notebook()

from sklearn.metrics import adjusted_rand_score


# Reading in constclust data

# In[2]:


params = pd.read_pickle("/Users/pwangel/Data/Single_Cell/Gut_Colon/params.pkl")
clusterings = pd.read_pickle("/Users/pwangel/Data/Single_Cell/Gut_Colon/clusterings.pkl")

data = sc.read_h5ad('/Users/pwangel/Data/Single_Cell/Gut_Colon/Colon_cell_atlas_scrublet.h5ad')


# In[17]:


sc.pl.umap(data, color=['donor', 'cell_type'], legend_loc='on data')


# In[3]:


clusterings.index = clusterings.index.astype(str)


# In[5]:


get_ipython().run_cell_magic('time', '', 'rec = reconcile(params, clusterings)\nrec')


# In[23]:


plotting.edge_weight_distribution(rec)


# In[6]:


get_ipython().run_cell_magic('time', '', 'comps = rec.get_components(0.98)')


# In[11]:


comps


# In[7]:


show(
    comps.filter(min_solutions=50, min_intersect=10).plot_hierarchies(data.obsm['X_umap']))


# In[26]:


comps.describe().sort_values("n_solutions", ascending=False)  # Looking at some statistics


# In[27]:


comps.describe().plot.scatter("n_solutions", "n_intersect")


# In[37]:


comps[:5].plot_components(data)


# In[35]:


comps[0].settings


# In[ ]:




