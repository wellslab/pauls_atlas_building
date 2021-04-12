#!/usr/bin/env python
# coding: utf-8

# In[30]:


import pandas as pd
import numpy as np

meta = pd.read_csv('/Users/pwangel/Downloads/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t', index_col=0)


# This notebooks extracts data from gtex according to sample tissue of origin. Data is from https://www.gtexportal.org/home/datasets. This script is necessary as gtex provides data as one single giant aggregate matrix which is difficult to read into memory. Here we are just breaking it up into smaller bits. 

# Sabrinas list of cell types (tissues) she wants

# In[27]:


categories = [[
'Brain - Putamen (basal ganglia)',
'Brain - Nucleus accumbens (basal ganglia)',
'Brain - Caudate (basal ganglia)',
'Brain - Cerebellum',
'Brain - Cerebellar Hemisphere',
'Brain - Anterior cingulate cortex (BA24)',
'Brain - Frontal Cortex (BA9)',
'Brain - Cortex',
'Brain - Hypothalamus',
'Brain - Hippocampus',
'Brain - Amygdala',
'Brain - Substantia nigra',
'Brain - Spinal cord (cervical c-1)',
'Nerve - Tibial'
],[
'Artery - Tibial',
'Artery - Coronary',
'Artery - Aorta',
'Heart - Atrial Appendage'
'Lung'
],[
'Small Intestine - Terminal Ileum',
'Colon - Sigmoid',
'Colon - Transverse',
'Esophagus - Muscularis',
'Esophagus - Gastroesophageal Junction',
'Stomach'
],[
'Pituitary',
'Thyroid',
'Adrenal Gland',
'Liver'
]]
names = ['CNS', 'Cardiovascular', 'ENS', 'Endocrine']


# In[70]:


for i_cat, i_name in zip(categories, names):
    samples = meta.loc[np.in1d(meta['SMTSD'].values, i_cat)&(meta['SMAFRZE']=='RNASEQ')].index.astype(str).tolist()
    data = pd.read_csv('/Users/pwangel/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', usecols=samples, skiprows=2)
    data.to_csv('/Users/pwangel/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm_%s.gct' %i_name, sep='\t')


# In[ ]:




