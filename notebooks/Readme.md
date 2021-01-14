## Folder structure

functions.py

Contains basic functions used to build the atlas. the dataframe or data object is a pandas dataframe with genes as rows and samples as columns. The annotations dataframe has samples as rows and variables as columns:

	transform_to_percentile(dataframe)
		Transforms expression values in dataframe to percentile values

	calculate_platform_dependence(data, annotations)
		Uses univariate linear modelling to determine the what fraction of the variance of each gene is attributable to platform. See blood atlas paper. 

	calculate_celltype_dependence(data, annotations, gene=None)	
		Uses univariate linear modelling to determine the what fraction of the variance of each gene is attributable to platform or celltype. See blood atlas paper.

	resample_clustering(data, annotations, resample_strategy, n_resamples=200, n_clusters_list=[3,4])
		Iterates randomised resampling and performs the atlas generation on each iteration. Uses H index to measure stability of clusters found (clustering algos are hard coded in for now).

	calc_H_index(clustering_dataframe)
		H index calculation function used in resample_clustering() functions. clustering_dataframe has particular format.

	KruskalWallisHTest(coords, annotations)
		Performs Kruskal Wallis H Test (grouped by platform) upon principal components. 

	plot_KW_Htest(data, annotations, varPart_df)
		Generate a plot the Kruskal Wallis H Test for the first 10 principal components while scanning a range of platform variance threshold cuts.

	plot_gene_platform_dependence_distribution(varPart_df)
		Plots the distribution of platform variance fraction

	plot_gene_platform_dependence_cumulative_distribution(varPart_df)
		Plots the cumulative distribution of platform variance fraction

	plot_pca(pca_coords, annotations, pca, labels, colour_dict)
		Plots the 3d atlas. Variables to plot and colour dictionary to use are optional - pass a list of variables into labels e.g. ['Dataset', 'celltype', 'Platform_Category']. 
 
