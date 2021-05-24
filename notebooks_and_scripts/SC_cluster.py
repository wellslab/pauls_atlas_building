import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import general_processing.SC_reader

# Read in a bunch of single cell datasets and cluster them. Aggregate clusters into 

# Read in here

gut_colon_data = general_processing.SC_reader.read_gut_atlas_colon()

# Find and remove doublets

sce.pp.scrublet(gut_colon_data)
gut_colon_data.write_h5ad("/scratch/yf0/pa5933/Data/Single_Cell/Gut_Colon/Colon_cell_atlas_scrublet.h5ad")

# Cluster using constclust

from constclust import cluster, reconcile, plotting

sc.tl.leiden(gut_colon_data)

gut_colon_params, gut_colon_clusterings = cluster(
    gut_colon_data,
    n_neighbors=np.linspace(15, 90, 6, dtype=np.int),
    resolutions=np.geomspace(0.01, 100, 40),
    random_state=[0, 1, 2],
    n_procs=6
)

pd.to_pickle(gut_colon_params, "/scratch/yf0/pa5933/Data/Single_Cell/Gut_Colon/gut_colon_params.pkl")
pd.to_pickle(gut_colon_clusterings, "/scratch/yf0/pa5933/Data/Single_Cell/Gut_Colon/gut_colon_clusterings.pkl")

# Aggregate into psuedo bulk samples. Take pickled results to local machine to plot and aggregate into bulk samples 
