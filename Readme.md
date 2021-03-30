## Folder structure

This repo contains the basic scripts Paul uses to run atlases. It does not contain scripts to pull data or annotations from stemformatics.

data/ contains example data files for the pluriptency atlas

notebooks_and_scripts/ contains Python scripts and .ipynb notebooks that are used to generate an atlas.

Many of the core functions are found in the functions.py file. There are multiple .ipynb notebooks using these functions.

Notable functions are:

transform_to_percentile()
    Transforms an expression matrix into percentile values.

calculate_platform_dependence()
    Runs univariate variance modelling on each gene in an expression matrix to calculate the fraction of it's variance that depends upon platform.

plot_pca()
    Plots the atlas using pca coordinates and sample annotations.

