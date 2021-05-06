#!/usr/bin/env python
# coding: utf-8

# This notebook is an overview and walkthough of the nascent pluripotency atlas. 
# 
# The prerequisites for creating the atlas are a) gene expression data for each sample, and b) metadata for each sample containing at a minimum the experimental platform the sample was measured on. These two files are usually called 'data' and 'annotations' respectively. 
# 
# This notebook also reads in some externally processed data, which we hope to process ourselves and include in the near future. These external datasets are recognisable by the GSE in their filename.

# In[1]:


import numpy as np
import pandas as pd
import sklearn
import functions
from general_processing.process_functions import convert_symbols_to_ensembl 


# Reading in expression data and metadata (annotations.tsv). For plotting, a 'display_metadata' field is required in the annotations dataframe, so I have used a temporary column here, the 'generic_sample_type'. To see the format of these dataframes just have a look at the example .tsv files I have here.

# In[60]:


data        = pd.read_csv('../data/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('../data/pluripotent_annotations.tsv', sep='\t', index_col=0)
data = data[annotations.index]

lizzis_anno    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcell_annotations.tsv', sep='\t', index_col=0)
#lizzis_anno.index = [str(i)+";"+str(j) for i, j in zip(lizzis_anno.chip_id, lizzis_anno.Dataset.astype(int))]
annotations = annotations.merge(lizzis_anno[['LM_Group_COLOR']], how='left', left_index=True, right_index=True)
annotations.LM_Group_COLOR = annotations[['LM_Group_COLOR']].fillna("Unannotated")
annotations.index = [str(int(i))+'_'+str(j) for i,j in zip(annotations.Dataset, annotations.chip_id)]
data.columns = annotations.index

day_anno = pd.read_csv('/Users/pwangel/Downloads/pluripotent_anno_feeder_day_year.tsv', sep='\t', index_col=0)
annotations = annotations.merge(day_anno[['Day', 'Year','Feeder']], how='left', left_index=True, right_index=True)

annotations['display_metadata'] = annotations.generic_sample_type
#annotations = annotations.loc[annotations.Platform_Category!='Illumina V2']


# In[56]:


annotations.shape


# This is a temporary hack to include some externally processed data and see how the atlas would look. They look good so I recommend they be processed with the stemformatics pipeline and included. I have not placed the external data into my git repo, it would be simpler just to process with s4m pipeline.

# In[61]:


df1 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE114873/GSE114873_Cufflinks_output_FPKM_and_Gene_Names.txt', sep='\t', index_col=0)[['condition', 'replicate', 'raw_frags']]
df1['sample_id'] = [i+'_'+j for i,j in zip(df1.condition.to_list(),df1.replicate.to_list())]
df1.drop(labels=['condition', 'replicate'], axis=1, inplace=True)
df1=df1.pivot(columns='sample_id')
df1.columns = df1.columns.droplevel()
df1 = convert_symbols_to_ensembl(df1)

df2 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE123055/GSE123055_counts.txt', sep='\t', index_col=0)
df2.drop(labels=['Unnamed: 52'], axis=1, inplace=True)

#This dataset is weird, not using it
#df3 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE131551/GSE131551_human_bulk.raw_count_matrix.tsv', sep='\t', index_col=0)
#df3.drop(labels='geneID', axis=1, inplace=True)

df3 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE137295/GSE137295_All_Raw_CountData.csv', index_col=0)

df4 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE167089/GSE167089_counts.tsv', sep='\t',index_col=0)


# As these external datasets do not have metadata, I created my own temporary external annotations dataframe.

# In[62]:


ext_annotations = pd.DataFrame(index=df1.columns.to_list()+df2.columns.to_list()+df3.columns.to_list()+df4.columns.to_list())
ext_annotations['Platform_Category']='RNASeq'
ext_annotations['Dataset'] = ['GSE114873' for i in range(df1.shape[1])] + ['GSE123055' for i in range(df2.shape[1])]                           + ['GSE137295' for i in range(df3.shape[1])] + ['GSE167089' for i in range(df4.shape[1])]
ext_annotations['generic_sample_type'] = 'Unannotated'
ext_annotations['display_metadata']  = [str(i_dataset)+'<br>'+str(i_sample) for i_dataset, i_sample in zip(ext_annotations.Dataset.values, ext_annotations.index.values)]


# Refine the data to only include those datasets appearing in the collabrorative annotations google docs file shared by the the biologists working on the project. These are the datasets of biological interest.
# 
# Also include the externally processed data in this step.

# In[63]:


#This tsv has the datasets of biological interest
datasets_to_keep = pd.read_csv('../data/Pluripotency datasets in stemformatics - existing stemformatics data.tsv', sep='\t')['Dataset'].values
annotations = annotations.loc[np.in1d(annotations.Dataset, datasets_to_keep)]

annotations = pd.concat([annotations, ext_annotations], sort=True)
data = data.merge(df1, how='inner', left_index=True, right_index=True)
data = data.merge(df2, how='inner', left_index=True, right_index=True)
data = data.merge(df3, how='inner', left_index=True, right_index=True)
data = data.merge(df4, how='inner', left_index=True, right_index=True)
data = data[annotations.index]


# Before we make the atlas, we will do some simple analysis on the data. In particular I would like to see the distribution of experimental platform. It would also be nice to see the cell type distributions, but without the proper metadata it is difficult

# In[66]:


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(annotations[['Dataset', 'Platform_Category']].drop_duplicates().groupby('Platform_Category').size())


# In[67]:


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(annotations.groupby(['Platform_Category']).size())


# Now to actually make the atlas. First step in the atlas two step process: transform expression values to percentile values.

# In[64]:


data = functions.transform_to_percentile(data)


# Second step: model the influence of platform upon expression for each gene. As this can take a while, I often save the results and just read them in rather than recompute them. In this case the results are saved in 'pluripotent_atlas_genes_with_ext.tsv'.

# In[28]:


#genes = functions.calculate_platform_dependence(data, annotations)
#genes.to_csv('../data/pluripotent_atlas_genes_with_ext.tsv', sep='\t')
#genes = pd.read_csv('../data/pluripotent_atlas_genes.tsv', sep='\t')
genes = pd.read_csv('../data/pluripotent_atlas_genes_with_ext.tsv', sep='\t') 


# Run the PCA on the expression data of the filtered, transformed genes. The value of the gene filter threshold is 0.25. I have not looked closely at this value. Perhaps a higher value would allow more components into the PCA. 

# In[65]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25]).transpose())


# I generate a separate set of coordinates for the external data as I would like to project them by themselves. The way this works in that a list of data/annotations dataframes is passed to the plot function. The first set of data/annotations is the base and subsequent sets are projected on. The plot the pca is saved as a .html in the <out_file> location.

# In[68]:


pca_coords_ext = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.25][ext_annotations.index]).transpose())

#First dataframes in the list of the base coordinates, following dataframes are projected on
functions.plot_pca([pca_coords, pca_coords_ext], [annotations, ext_annotations],pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset'], colour_dict={}, pcs=[1,2,3],                    out_file='/Users/pwangel/Downloads/pluripotent_atlas_with_external.html')


# Now try to 'zoom in' on the pluripotent cells (isolate them by applying k means clustering). This is a fairly rough way to identify the samples that are relevant to the 'naive' vs 'primed' analysis. I want stem cells only, no differentiated samples, this is best cleared up by the biological annotations but k means will do for now.

# In[69]:


kmeans = sklearn.cluster.KMeans(n_clusters=4).fit(pca_coords)
annotations['K Means'] = kmeans.labels_
ext_annotations['K Means'] = annotations['K Means'].loc[ext_annotations.index]


# Plot the PCA again but now with the kmeans clusters, so we can identify the biology of each cluster.

# In[70]:


functions.plot_pca(pca_coords, annotations,pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset', 'K Means', 'LM_Group_COLOR'],                    colour_dict={"Unannotated":'grey'}, pcs=[1,2,3],                    out_file='/Users/pwangel/Downloads/pluripotent_kmeans_atlas_with_external.html')


# Essentially a repeat of the previous process to generate an atlas, but now selecting only stem cells first. WARNING: the number label of the stemcells group in the kmeans clustering is random each time you run it. 

# In[71]:


pluripotent_annotations = annotations.loc[(annotations['K Means']==0) | np.in1d(annotations.Dataset.values, [7253, 7240, 7124, 7135, 'GSE137295', 'GSE123055', 'GSE114873', 'GSE167089'])]
pluripotent_data = data[pluripotent_annotations.index]
#pluripotent_genes = functions.calculate_platform_dependence(pluripotent_data, pluripotent_annotations)
#pluripotent_genes.to_csv('../data/pluripotent_only_atlas_genes_with_ext.tsv', sep='\t')
pluripotent_genes = pd.read_csv('../data/pluripotent_only_atlas_genes_with_ext.tsv', sep='\t')


# In[72]:


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(pluripotent_annotations[['Dataset', 'Platform_Category']].drop_duplicates().groupby('Platform_Category').size())


# In[73]:


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(pluripotent_annotations.groupby(['Platform_Category']).size())


# In[74]:


pluripotent_annotations.generic_sample_type.value_counts()


# In[23]:


with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(pluripotent_annotations.groupby(['Platform_Category','generic_sample_type']).size())


# Save this subset for future fun

# In[14]:


pluripotent_annotations.to_csv('../data/pluripotent_RNASeq_annotations.tsv', sep='\t')


# Run the PCA again, using a threshold of 0.2. The general pattern with these thresholds is that the narrower the biological range of samples is, the stricter the threshold must be. Lower biological variation means a greater fraction of platform variation.

# In[75]:


pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(pluripotent_data.loc[pluripotent_genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(pluripotent_data.loc[pluripotent_genes.Platform_VarFraction.values<=0.2]).transpose())


# In[76]:


annotations.Year.unique()


# Plot the stemcell only PCA (and save it).

# In[82]:


plot_pca(pca_coords, pluripotent_annotations,pca,                    labels=['generic_sample_type', 'Platform_Category', 'Dataset', 'LM_Group_COLOR', 'Day', 'Year', 'Feeder'],                    colour_dict={"Unannotated":"grey"}, pcs=[1,2,3],                    out_file='/Users/pwangel/Downloads/pluripotent_only_atlas_with_external.html')


# In[80]:


from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *
import plotly.figure_factory as ff
import plotly.io


# In[81]:


more_colour_list = ['red', 'green', 'black', 'yellow', 'brown', 'coral', 'sandybrown', 'darkorange', 'olive', 'limegreen', 'lightseagreen']

def plot_pca(pca_coords, annotations, pca, labels, colour_dict, pcs=[1,2,3], out_file=None):

    '''
    Generate PCA plot of Atlas

    Parameters
    ----------

    pca_coords
        coordinates of each sample in pca space, n_samplesxn_dimensions sized np array. Can be a list of such arrays

    annotations
        Metadata datafrane, index as samples, columns as properties 
        Must have a 'Platform_Category' column
        Can be a list of such dataframes

    pca
        sklearn pca, already been fitted

    labels
        list of metadata labels to plot, e.g. ['Dataset', 'celltype', 'Platform_Category']

    colour_dict
        colours to used to plot in a dictionary format

    '''

    if isinstance(pca_coords, np.ndarray):
        pca_coords = [pca_coords]
        annotations = [annotations]

    visibility_df = pd.DataFrame(columns=['type', 'label'])

    counter = 0
    for i_label in labels:
        for i_annotation in annotations:
            if (i_annotation[i_label].dtype==np.float64()) and (i_label!='Dataset'):
                visibility_df.loc[counter] = [i_label+'Cont', i_label]
                counter+=1
                visibility_df.loc[counter] = [i_label+'Nan', i_label]
                counter+=1
            else:
                for i_type in i_annotation[i_label].unique().astype(str):
                    visibility_df.loc[counter] = [i_type, i_label]
                    counter+=1

    extended_colour_dict = {visibility_df.type.unique()[i_key]:more_colour_list[i_key%len(more_colour_list)] for i_key in                             range(visibility_df.type.unique().shape[0]) if visibility_df.type.unique()[i_key] not in colour_dict.keys()}
    colour_dict.update(extended_colour_dict)

    fig = Figure()

    button_list = []
    for i_label in labels:

        visibility_list = (visibility_df.label.values==i_label).tolist()

        for i in range(len(annotations)):
            i_annotation = annotations[i]
            i_pca_coords = pca_coords[i]

            if (i==0) and len(annotations)>1:
                opac = 0.3
            else:
                opac = 0.9

            if (i_annotation[i_label].dtype==np.float64()) and (i_label!='Dataset'):

                sel_not_nan = ~pd.isnull(i_annotation[i_label]).values   
                fig.add_trace(Scatter3d(x=i_pca_coords[sel_not_nan,pcs[0]-1], y=i_pca_coords[sel_not_nan,pcs[1]-1], z=i_pca_coords[sel_not_nan,pcs[2]-1], 
                    mode='markers', text=[str(i)+'<br>'+str(j) for i, j in zip(i_annotation[i_label].values[sel_not_nan], i_annotation.display_metadata.values[sel_not_nan])], 
                    opacity=opac, name=i_label, visible=False, 
                    marker=dict(size=5, color=i_annotation[i_label].values[sel_not_nan],
                    colorscale = 'viridis', colorbar=dict(thickness=20))))

                fig.add_trace(Scatter3d(x=i_pca_coords[~sel_not_nan,pcs[0]-1], y=i_pca_coords[~sel_not_nan,pcs[1]-1], z=i_pca_coords[~sel_not_nan,pcs[2]-1], 
                    mode='markers', text=i_annotation.display_metadata.values[~sel_not_nan], 
                    opacity=0.3, name=i_label, visible=False, 
                    marker=dict(size=5, color='grey')))
            else:
                for i_type in i_annotation[i_label].unique().astype(str):
    
                    sel = i_annotation[i_label].values.astype(str) == i_type 
                    i_colour = colour_dict[i_type]

                    fig.add_trace(Scatter3d(x=i_pca_coords[sel,pcs[0]-1], y=i_pca_coords[sel,pcs[1]-1], z=i_pca_coords[sel,pcs[2]-1], 
                        mode='markers', text=i_annotation.display_metadata.values[sel], 
                        opacity=opac, name=i_type, visible=False, 
                        marker=dict(size=5, color=i_colour)))
        
        button_list.append(dict(label=i_label,
                                method="update",
                                args=[{"visible": visibility_list},
                                    {"title": i_label}]))
      
        fig.update_layout(
            updatemenus=[dict(active=0,buttons=button_list,)],
            scene = dict(xaxis_title='PC%d %%%.2f' %(pcs[0], pca.explained_variance_ratio_[pcs[0]-1]*100),
                         yaxis_title='PC%d %%%.2f' %(pcs[1], pca.explained_variance_ratio_[pcs[1]-1]*100),
                         zaxis_title='PC%d %%%.2f' %(pcs[2], pca.explained_variance_ratio_[pcs[2]-1]*100))
            )

    if out_file is not None:
        plot(fig, auto_open=False, filename=out_file) 
    else:
        iplot(fig)


# In[ ]:




