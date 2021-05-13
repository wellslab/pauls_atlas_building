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


