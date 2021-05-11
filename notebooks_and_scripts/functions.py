import numpy as np
import pandas as pd
import gc, os, psutil, sys, scipy, sklearn
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sklearn.decomposition
from sklearn.cluster import KMeans

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *
import plotly.figure_factory as ff
import plotly.io

def transform_to_percentile(dataframe):

    '''
    Apparently this is properly called the spearman rank

    Parameters:
    ----------

    dataframe
        Dataframe containing expression values, index as variables (genes), columns as samples

    Returns:
    -----------

    transformed_dataframe
        Dataframe with expression as rank (percentile) values

    '''

    transformed_dataframe = (dataframe.shape[0] - dataframe.rank(axis=0, ascending=False, na_option='bottom')+1)/(1+dataframe.shape[0])

    return transformed_dataframe


def calculate_platform_dependence(data, annotations):

    '''
    Calculates the fraction of variance due to platform

    Parameters:
    ----------

    data
        Dataframe containing expression values, index as variables (genes), columns as samples

    annotations
        Dataframe containing metadata, index as samples, requires a column 'Platform_Category'

    Returns:
    ----------

    output_df
        Dataframe contains fraction of variance due to platform

    '''
 
    import warnings
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    warnings.simplefilter('ignore', ConvergenceWarning)

    output_df = pd.DataFrame(index=data.index, columns=['Platform_VarFraction'])
    data = data.transpose()
    data['Platform_Category'] = annotations['Platform_Category'].values

    for i_gene in data.columns.values[:-1]:

        md  = smf.mixedlm("%s ~ Platform_Category" %str(i_gene), data=data, groups = data['Platform_Category'])
        mdf = md.fit()
        output_df.loc[i_gene, 'Platform_VarFraction'] = mdf.fittedvalues.std()**2/(mdf.fittedvalues.std()**2+mdf.resid.std()**2)

    data.drop(labels = 'Platform_Category', axis=1, inplace=True)
    data = data.transpose()

    return output_df

def calculate_celltype_dependence(data, annotations, gene=None, pval=False):

    '''
    Calculates the fraction of variance due to celltype and platform. Copied from Yidi Deng's work. BUT NOT WORKING AT THE MOMENT.

    Parameters:
    ----------

    data
        Dataframe containing expression values, index as variables (genes), columns as samples

    annotations
        Dataframe containing metadata, index as samples, requires a column 'Platform_Category'

    gene
        either a single gene to analyse, or a list, defaults to doing all genes in data

    Returns:
    ----------

    output_df
        Dataframe contains fraction of variance due to platform

    '''

    if isinstance(gene, str):
        genes_to_iterate = [gene]
    elif gene is None:
        genes_to_iterate = data.index.values
    else:
        genes_to_iterate = gene

    output_df = pd.DataFrame(index=genes_to_iterate, columns=['Platform_VarFraction', 'celltype_VarFraction', 'celltype_pval'])
    temp_data = data.copy().transpose().merge(annotations[['Platform_Category', 'celltype']], how='left', left_index=True, right_index=True)

    temp_data["group"] = 1                                                                                                            
    vcf = {"celltype": "0 + C(celltype)", "Platform_Category": "0 + C(Platform_Category)"} 
    vcf_nocelltype = {"Platform_Category": "0 + C(Platform_Category)"}

    for i in range(len(genes_to_iterate)):

        i_gene = genes_to_iterate[i]
        print(i, i_gene) 

        md = sm.MixedLM.from_formula("%s ~ 1" %str(i_gene), groups="group",                                                    
                                vc_formula=vcf, re_formula=None, data=temp_data)  
        mdf = md.fit(reml=True,method="powell")

        output_df.loc[i_gene, 'Platform_VarFraction'] = mdf.vcomp[0]/(mdf.vcomp.sum()+mdf.scale)
        output_df.loc[i_gene, 'celltype_VarFraction'] = mdf.vcomp[1]/(mdf.vcomp.sum()+mdf.scale)

        if pval:

            md_nocelltype  = sm.MixedLM.from_formula("%s ~ 1" %str(i_gene), groups="group",                                                    
                                vc_formula=vcf_nocelltype, re_formula=None, data=temp_data)  
            mdf_nocelltype = md.fit(reml=True,method="powell")

            output_df.loc[i_gene, 'celltype_pval'] = sm.stats.anova_lm(mdf, mdf_nocelltype)['PR(>F)'] 

    if not pval:
        output_df.drop(labels=['celltype_pval'], axis=1)

    return output_df

def resample_clustering(data, annotations, resample_strategy, n_resamples=200, n_clusters_list=[3,4], ids_to_isolate=None, base_output_file=None):
    
    ''' 
    This cumbersome function performs either jack-knife or bootstrap resampling. 
    WARNING: At each iteration it will recalculate gene platform dependence. So it takes 1 minutes to platform dependence, and 200 resamples are performed, 
    it will take 200x1 = 200 minutes (a bit long...).
 
    Parameters
    ---------

    data
        Expression values, index as genes, columns as samples

    annotations
        Metadata, index as samples, columns as properties 
        Must have a 'Platform_Category' column and a 'Dataset' column

    resample_strategy
        Either 'bootstrap', 'jackknife' or None. If it is None then it will simply perform the clustering algorithm n_resamples times.

    n_resamples
        Number of bootstrap resampling iterations if resample_stragety is 'bootstrap'

    n_clusters_list
        List of cluster parameters, each value is tested for clustering stability independently

    id_to_isolate
        Will only perform clustering upon this subset of samples

    base_output_file
        File to save the base clustering values to

    Returns
    ----------

    results
        Numpy list of arrays, each array contains the H index for each cluster

    retained_genes_list
        Numpy list of arrays, each array contains the retained genes for each iteration of the resampling

    '''
    print("Performing resampling upon data of shape:", data.shape)

    # Initial search for platform dependent genes
 
    base_platform_dependence = calculate_platform_dependence(data, annotations)  
 
    base_genes  = base_platform_dependence.index.values[base_platform_dependence.Platform_VarFraction.values<=0.2]
    base_data   = transform_to_percentile(data.loc[base_genes].copy())
    pca         = sklearn.decomposition.PCA(n_components=3, svd_solver='full')
    base_output = pd.DataFrame(pca.fit_transform(base_data.transpose()), index=base_data.columns, columns = ['PC1', 'PC2', 'PC3'])
    
    print("Utlising %d genes as baseline expression data\n" %base_genes.shape[0])

    retained_genes_list = [base_genes]

    resamples_clusters_list = []
    for i_clusters in n_clusters_list:

        if ids_to_isolate is None:
            clustering_df = pd.DataFrame(KMeans(n_clusters=i_clusters).fit_predict(base_output), index=base_data.columns, columns=['Base'])
        else:
            clustering_df = pd.DataFrame(KMeans(n_clusters=i_clusters).fit_predict(base_output.loc[ids_to_isolate]), index=base_output.loc[ids_to_isolate].index, columns=['Base'])

        #This is optional to save the base clusters
        if base_output_file is not None:
            clustering_df.to_csv(base_output_file+"_%d" %int(i_clusters))
        resampled_clusters_list.append(clustering_df)    

    if (resample_strategy=='bootstrap') or (resample_strategy is None):
        iterations = np.arange(n_resamples)
    elif resample_strategy=='jackknife':
        iterations = np.arange(annotations['Dataset'].unique().shape[0])

    #### Do Jacknife
    print("Starting resampling\n")

    for i_iter in iterations:
           
        if resample_strategy=='jackknife':
 
            print("Omitting dataset %d" %annotations['Dataset'].unique()[i_iter])
            i_annotations = annotations.copy().loc[annotations['Dataset'] != annotations['Dataset'].unique()[i_iter]]
            i_data        = transform_to_percentile(data[i_annotations.index.values].copy())
            i_varPart  = calculate_platform_dependence(i_data, i_annotations) 
            i_cut_data = transform_to_percentile(i_data.loc[i_varPart.loc[i_varPart['Platform_VarFraction']<=0.2].index.values])  
         
        elif resample_strategy=='bootstrap':
 
            print("Bootstrap resampling number %d" %i_iter)
            i_annotations = annotations.copy().sample(frac=1.0, replace=True)
            i_data        = transform_to_percentile(data[i_annotations.index.values].copy())
            i_varPart  = calculate_platform_dependence(i_data, i_annotations) 
            i_cut_data = transform_to_percentile(i_data.loc[i_varPart.loc[i_varPart['Platform_VarFraction']<=0.2].index.values])  

        elif resample_strategy is None:
 
            print("Non resampling number %d" %i_iter)
            i_annotations = annotations.copy()
            i_cut_data    = base_data.copy()
            i_varPart  = base_platform_dependence.copy()
    
        sys.stdout.flush()
 
        i_output   = pd.DataFrame(pca.fit_transform(i_cut_data.transpose()), index=i_cut_data.columns, columns = ['PC1', 'PC2', 'PC3'])
    
        for i in range(len(n_clusters_list)):
  
            if ids_to_isolate is None:
                i_clustering_df = pd.DataFrame(KMeans(n_clusters=i).fit_predict(i_output).astype(int), index=i_output.index, columns=[i_iter])
            else:
                i_clustering_df = pd.DataFrame(KMeans(n_clusters=i).fit_predict(i_output.loc[ids_to_isolate]).astype(int), index=i_output.loc[ids_to_isolate].index, columns=['i_iter'])

            i_clustering_df.drop_duplicates(inplace=True)
            resampled_clusters_list[i] = resampled_clusters_list[i].merge(i_clustering_df, how='left', left_index=True, right_index=True)
    
        gc.collect()    
    
        retained_genes_list.append(i_varPart.loc[i_varPart['Platform_VarFraction']<=0.2].index.values)

    results = [calc_H_index(resampled_clusters_list[i]) for i in range(len(n_clusters_list))]
    
    print('Done\n')

    return results, retained_genes_list

def calc_H_index(clustering_dataframe):

    '''
    Calculate H index of from a dataframe of resampled iterations

    Parameters
    ---------

    clustering_dataframe
        The first column of the dataframe should be the 'true' reference set
        Agnostic as to the class labels used and they don't have to be consistent across resampling iterations

    Returns
    ---------
        Numpy array containing the H index for each cluster in the true reference set

    '''

    h_index_per_cluster = []

    for i_cluster in sorted(clustering_dataframe.iloc[:,0].unique()):

        samples_in_cluster_i = clustering_dataframe.loc[clustering_dataframe.iloc[:,0].values == i_cluster].index.values
        max_similarity_list  = []

        for i_resample in range(1, clustering_dataframe.shape[1]):

            # Should only compare to samples that were present in the resample (samples may be ommitted in a given iteration of resampling)
            samples_to_compare   = np.intersect1d(samples_in_cluster_i, clustering_dataframe.iloc[:, i_resample].dropna().index.values) 

            temp_similarity = []

            for j_cluster in clustering_dataframe.iloc[:,i_resample].dropna().unique():

                samples_in_cluster_j = clustering_dataframe.loc[clustering_dataframe.iloc[:,i_resample].values == j_cluster].index.values

                temp_similarity.append(np.intersect1d(samples_in_cluster_j, samples_to_compare).shape[0]/np.union1d(samples_in_cluster_j, samples_to_compare).shape[0])

            max_similarity_list.append(max(temp_similarity))

        h_index = max([h for h in np.arange(0,1.0,0.001) if (np.array(max_similarity_list)>=h).astype(int).sum()/(clustering_dataframe.shape[1]-1)>=h])

        h_index_per_cluster.append(h_index)    

    return np.array(h_index_per_cluster)

def KruskalWallisHTest(coords, annotations):

    '''
    Uses the Kruskal Wallis H Test as a proxy for the dependence of a single principal component on platform.

    Parameters
    ----------

    coords
        Numpy array, coordinates of each sample along a principal component.

    annotations
        Dataframe with sample metadata, must have 'Platform_Category' as a column

    Returns
    ----------

    kruskal
        Float, the result of the Kruskal Wallis H Test

    '''

    groups = {}

    for i_platform in annotations['Platform_Category'].unique():
        sel = annotations['Platform_Category']==i_platform
        groups[i_platform] = coords[sel]

    args = groups.values()
    
    return scipy.stats.kruskal(*args)[0]  


def plot_KW_Htest(data, annotations, varPart_df):

    '''
    Varies the platform dependence threshold from 0.02->1.0 
    Calculates the strenth of the platform dependence upon a the first 10 principal components as the threshold changes
    Plots results as a heatmap

    Parameters
    ----------

    data
        Expression values, index as genes, columns as samples

    annotations
        Metadata, index as samples, columns as properties 
        Must have a 'Platform_Category' column

    varPart_df
        Dataframe contains fraction of variance due to platform under the column 'Platform_VarFraction'

    Returns
    ----------

    platform_kruskal
        Dataframe containing the KW H test results for each of the first 10 componenets, for each of the platform variance thresholds iterated over

    '''

    print("Assessing platform dependence for principal components with varying threshold.")

    # Loop through cut values
    n_components   = 10
    n_thresholds   = 50

    threshold_list = np.arange(1, n_thresholds)/n_thresholds
    platform_kruskal = pd.DataFrame(data = np.empty(shape=(threshold_list.shape[0], n_components)), index=threshold_list, columns=np.arange(1,n_components+1))
    n_genes_list = []

    for i_thresh in range(threshold_list.shape[0]):

        cut_thresh    = threshold_list[i_thresh]
        sel_varPart   = varPart_df.Platform_VarFraction.values <= cut_thresh
        genes_to_keep = data.loc[sel_varPart].index.values

        print("Analysing threshold of %f (%d genes)" %(cut_thresh, genes_to_keep.shape[0]))

        transformed_data = transform_to_percentile(data.loc[genes_to_keep].copy())

        pca    = sklearn.decomposition.PCA(n_components=n_components)
        output = pca.fit_transform(transformed_data.transpose())

        for i_component in range(n_components): 
            platform_kruskal.iloc[i_thresh, i_component] = KruskalWallisHTest(output[:,i_component], annotations)/output.shape[0]

    fig = Figure(data=Heatmap(z=platform_kruskal.values, x=np.arange(1,n_components+1), y=threshold_list, colorscale = 'Viridis'), 
                layout=Layout(title="KW H test",xaxis_title="Component",yaxis_title="Platform Variance Fraction Threshold", yaxis_nticks=10, width=700, height=700,
                              autosize = False))
    iplot(fig)

    return platform_kruskal

    
def plot_gene_platform_dependence_distribution(varPart_df):

    '''
    Parameters
    ----------

    varPart_df
        Dataframe contains fraction of variance due to platform under the column 'Platform_VarFraction'
    '''

    #Distribution of platform dependence variance fraction
    fig = Figure(data=[Histogram(x=varPart_df.Platform_VarFraction.values, xbins=dict(start=0.0,end=1.0,size=0.02))], 
                layout=Layout(xaxis_title="Platform Variance Fraction Threshold", yaxis_title="N Genes", width=700, height=700))
    iplot(fig)

def plot_gene_platform_dependence_cumulative_distribution(varPart_df):

    '''
    Parameters
    ----------

    varPart_df
        Dataframe contains fraction of variance due to platform under the column 'Platform_VarFraction'
    '''

    #Number of genes passing cut
    vals, bins = np.histogram(varPart_df.Platform_VarFraction.values, bins=50)
    fig = Figure(data=[Scatter(x=0.5*(bins[:-1]+bins[1:]), y=np.cumsum(vals))], 
                layout=Layout(xaxis_title="Platform Variance Fraction Threshold", yaxis_title="Number of genes passing cut", width=700, height=700))
    fig.update_xaxes(autorange="reversed")
    iplot(fig)


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

    extended_colour_dict = {visibility_df.type.unique()[i_key]:more_colour_list[i_key%len(more_colour_list)] for i_key in \
                            range(visibility_df.type.unique().shape[0]) if visibility_df.type.unique()[i_key] not in colour_dict.keys()}
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


def plot_dot_plots(dataframe, annotations, cell_property, cell_colour_by, gene_list, output_dir, plot_hist=False):

    import seaborn as sns
    import matplotlib.pyplot as pyplot
    sns.set(style="ticks")
    
    for i_gene in gene_list:
   
        print(i_gene) 
        df_violin  = pd.DataFrame(columns=['Expression', 'Cell Property', 'Cell Colour'])
        df_violin['Expression'] = dataframe.loc[i_gene].values.flatten().astype(float)
        df_violin['Cell Property'] = annotations[cell_property].values.flatten()
        df_violin['Cell Colour'] = annotations[cell_colour_by].values.flatten()

        if plot_hist==True:
            fig, (ax0, ax1) = pyplot.subplots(2,1, figsize=(13.5,12.0))
            ax1.set_xlim(df_violin['Expression'].min(), df_violin['Expression'].max())
        else:
            fig, ax0 = pyplot.subplots(1,1, figsize=(13.5,6.0))

        sns.swarmplot(x="Expression", y="Cell Property", size=9, hue='Cell Colour', data=df_violin, ax=ax0, palette="Set2")
        if plot_hist==True:
            sns.histplot(data=df_violin, x="Expression", ax=ax1)
    
        ax0.set_title(i_gene+' Expression')
        ax0.set_xlim(df_violin['Expression'].min(), df_violin['Expression'].max())
        #pyplot.show()
        pyplot.savefig(output_dir+'/%s_dotplot.pdf' %i_gene)
        pyplot.close()
