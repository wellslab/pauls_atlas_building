import scanpy as sc
import numpy as np
import pandas as pd

data_location = '/scratch/yf0/pa5933/Data/'

# This is a collection of functions used to read in single cell data. It is organised in this way as there is no real standardised format, and most datasets are different. Therefore I have a function for each dataset

def read_GSE132044():

    min_library_size = 5000
    min_genes = 1000
    
    folder = data_location+'/Single_Cell/Ding/'
    fname  = data_location+'/Single_Cell/Ding/counts.read.txt'
    
    #Read in and form clusters
    data  = sc.read_mtx(fname)
    data = sc.AnnData(X=data.X.transpose())
    data.uns['min_library_size'] = min_library_size
    data.uns['min_genes'] = min_genes
    data.uns['folder'] = folder
    
    samples = pd.read_csv(folder+"/cells.read.txt", skiprows=0, header=None, sep='\t')
    genes = pd.read_csv(folder+"/genes.read.txt", skiprows=0, header=None, sep='\t')
    genes = [item.split("_")[0] for item in genes[0].values.astype(str)]
    
    name_map = pd.read_csv(folder+"/map.CPM.names.Count.names.txt", sep='\t', index_col=0)
    samples_meta = pd.read_csv(folder+"/meta.txt", sep="\t", index_col=0)
    
    name_map = name_map.merge(samples_meta, how='left', left_index=True, right_index=True)
    samples  = samples.merge(name_map, how='left', left_on=0, right_index=True).set_index(0) 
    
    #The authors do not provide meta data for all samples, presumably they are excluded from the analysis for a good reason
    sel         = ~samples.Method.isnull().values & ~samples.CellType.isnull().values & (np.array(data.X.sum(1)).reshape(-1)>=min_library_size) \
                  & (np.array(data.X.astype(bool).sum(axis=1)).reshape(-1) >= min_genes)

    return None


def read_gut_atlas_colon():

    return sc.read_h5ad(data_location+'/Single_Cell/Gut_Colon/Colon_cell_atlas.h5ad')


def aggregate_clusters(anndata, key, n_samples_per_cluster=20, average_library_size_per_cluster=None):

    import random
    df_output = pd.DataFrame(index=anndata.var_names.values)
    df_annotations = pd.DataFrame()
    metadata = []

    for i_cluster in np.unique(anndata.obs[key].values):

        print("Processing cluster %s" %str(i_cluster))
        sel = anndata.obs[key].values==i_cluster
        columns = list(anndata.obs_names.values[sel])

        i_df = pd.DataFrame.sparse.from_spmatrix(anndata.X[sel,:], columns=anndata.var_names.values, index=columns).transpose()

        #If we are targetting an average library size per cluster, set n_samples_per_cluster to be equivalent to this
        if average_library_size_per_cluster is not None:
            n_samples_per_cluster = int(i_df.shape[1]*average_library_size_per_cluster/i_df.sum().sum())

        i=0
        while len(columns)>n_samples_per_cluster:
            selectedColumns = random.sample(columns, n_samples_per_cluster)
            df_output["%s_%s" %(str(i_cluster),str(i))] = i_df[selectedColumns].sum(axis=1).values
            metadata.append(i_cluster)
            columns  = set(columns).difference(set(selectedColumns))
            i += 1

    df_annotations[key] = metadata
    df_annotations.index=df_output.columns

    return df_output, df_annotations            



