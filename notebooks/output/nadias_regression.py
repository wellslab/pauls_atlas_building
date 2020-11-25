#!/usr/bin/env python
# coding: utf-8

# In[39]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[40]:


blood_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/imac_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}


# Reading in data, including nadias annotations, excel spreadsheet with multiple tabs

# In[41]:


data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_annotations.tsv', sep='\t', index_col=0)
genes          = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t', index_col=0)
nadias_annotations = pd.read_excel('/Users/pwangel/Downloads/Regression Analysis Updated.xlsx', sheet_name=None)


# Reformat nadias annotations so each variable is represented in one column (for easier plotting)

# In[42]:


for i_key in nadias_annotations.keys():
    i_sheet = nadias_annotations[i_key]
    i_sheet.index = [str(i_dataset)+";"+str(int(i_sample)) for i_dataset, i_sample in                      zip(i_sheet['Sample'].values, i_sheet['Dataset ID'].values)]
    i_sheet = i_sheet.drop(['Dataset ID', 'Sample', 'Cell type'], axis=1)
    i_sheet[i_key] = 'Unannotated'
    for i_col in i_sheet.drop([i_key], axis=1).columns.values:
        i_sheet[i_key].loc[i_sheet[i_col]=='X'] = i_col
    annotations = annotations.merge(i_sheet[i_key], how='left', left_index=True, right_index=True)
    
annotations.fillna('Unannotated', inplace=True)    


# Make some combinations of Nadias annotations, i.e. TPO and Hypoxic vs ... etc etc

# In[43]:


for i_key in range(len(nadias_annotations.keys())):
    for j_key in range(i_key+1, len(nadias_annotations.keys())):
        i_name = list(nadias_annotations.keys())[i_key]+'__'+list(nadias_annotations.keys())[j_key]
        annotations[i_name] = [i+'__'+j for i,j in zip(annotations[list(nadias_annotations.keys())[i_key]],annotations[list(nadias_annotations.keys())[j_key]])]
annotations.replace(to_replace='Unannotated__Unannotated', value='Unannotated', inplace=True)


# In[44]:


data = functions.transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already

# In[7]:


#genes = functions.calculate_platform_dependence(data, annotations)
#genes.to_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t') 


# In[45]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())


# In[47]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['celltype', 'Platform_Category', 'Dataset']+list(annotations.columns[31:36].values), colour_dict=blood_atlas_colours)


# In[ ]:


functions.plot_KW_Htest(data, annotations, genes)


# In[ ]:




