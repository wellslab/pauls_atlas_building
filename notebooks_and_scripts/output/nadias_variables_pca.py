#!/usr/bin/env python
# coding: utf-8

# Make PCA graphs but colour only certain samples and colour them according to data provided by Nadia. This data is provided as ms excel sheets which must be read in

# In[115]:


import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
from general_processing.processing_functions import convert_symbols_to_ensembl, transform_to_percentile 


# In[116]:


blood_atlas_colours = pd.read_csv('/Users/pwangel/Data/Metadata_dumps/imac_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}


# Reading in data, including nadias annotations, excel spreadsheet with multiple tabs

# In[117]:


data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_annotations.tsv', sep='\t', index_col=0)
genes          = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t', index_col=0)
nadias_annotations = pd.read_excel('/Users/pwangel/Downloads/Regression Analysis Latest.xlsx', sheet_name=None, header=0)


#Want to update the cell type variable (actually called the new one cell_type)
new_annotations = pd.read_csv('/Users/pwangel/Downloads/myeloid Tier update - samples_myeloid_merged.tsv', sep='\t')
new_annotations.index = [str(i_dataset)+";"+str(int(i_sample)) for i_dataset, i_sample in                      zip(new_annotations['sample_id'].values, new_annotations['dataset_id'].values)]
annotations = annotations.merge(new_annotations['cell_type'], left_index=True, right_index=True)


# Reformat nadias annotations so each variable is represented in one column (for easier plotting). PLEASE DO NOT HAVE THE SHEET NAME THE SAME AS A COLUMN NAME

# In[118]:


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

# In[119]:


for i_key in range(len(nadias_annotations.keys())):
    for j_key in range(i_key+1, len(nadias_annotations.keys())):
        i_name = list(nadias_annotations.keys())[i_key]+'__'+list(nadias_annotations.keys())[j_key]
        annotations[i_name] = [i+'__'+j for i,j in zip(annotations[list(nadias_annotations.keys())[i_key]],annotations[list(nadias_annotations.keys())[j_key]])]
annotations.replace(to_replace='Unannotated__Unannotated', value='Unannotated', inplace=True)


# In[120]:


data = transform_to_percentile(data)


# Only need to compute gene variance fraction if not done already

# In[121]:


#genes = functions.calculate_platform_dependence(data, annotations)
#genes.to_csv('/Users/pwangel/Downloads/myeloid_atlas_genes.tsv', sep='\t') 


# In[122]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())


# In[123]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['cell_type','Dataset']+list(nadias_annotations.keys()), colour_dict=blood_atlas_colours)


# This section is showing microglia only

# In[124]:


for i_col in list(annotations.columns[31:38].values):
    annotations.loc[annotations.cell_type!='microglia', i_col] = 'Unannotated'


# In[125]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['cell_type', 'Platform_Category', 'Dataset']+list(annotations.columns[31:38].values), colour_dict=blood_atlas_colours)


# In[ ]:




