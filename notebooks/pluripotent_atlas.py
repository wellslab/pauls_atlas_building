import numpy, pandas, sys, os, scipy.stats, sklearn.decomposition
sys.path.append("/Users/pwangel/Gene_Analysis")

# Generate an atlas using pluripotent stem cell data
# Annotations expected to be updated in the near future

# Read in and find cut
data = pandas.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pandas.read_csv('/Users/pwangel/Downloads/pluripotent_annotations.tsv', sep='\t', index_col=0)

data = functions.transform_to_percentile(data)
genes = functions.calculate_platform_dependence(data, annotations)

pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())
pca_coords = pca.transform(functions.transform_to_percentile(data.loc[genes.Platform_VarFraction.values<=0.2]).transpose())

functions.plot_pca(pca_coords, annotations,pca, \
                   labels=['celltype', 'Platform_Category', 'Dataset']+list(nadias_annotations.keys()), colour_dict=blood_atlas_colours)
