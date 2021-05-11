import pandas as pd
import numpy as np

main_ensembl_ids = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/ensembl_hg38.91_chromosome.csv').ensembl_gene_id.values.astype(str)

def remove_microarray_duplicates(data, maps):

    '''
    Removes duplicates probe-gene maps from microarray.
    If there are multiple probes for a single gene, the probe with the highest mean is chosen.
    If there are multiple genes for a single probe, the probe is discarded.

    Parameters:
    ----------

    data
        Dataframe containing expression values, index as probe, columns as samples

    maps
        Dataframe containing probe-gene maps, index is probe, column is gene

    Returns:
    -----------

    data
        Expression dataframe with variables converted to gene

    '''

    print("\nGoing to remove duplicates from dataset containing %d probes\n" %len(data))

    maps = maps.loc[numpy.intersect1d(maps.index.values, data.index.values)]
    data = data.assign(Mean=data.mean(axis=1))

    print("> Dataset has %d (%d unique) probes before removing the devious false haplotype mapping ids" %(len(maps), len(maps.index.unique())))

    maps = maps.loc[numpy.in1d(maps.Gene.values.astype(str), main_ensembl_ids)]
    print("> Mapping has been reduced to %d (%d unique) probes" %(maps.shape[0], len(maps.index.unique().values)))

    maps = maps.loc[maps.index.value_counts()[maps.index.value_counts()==1].index.values]
    print("> Mapping has been reduced to %d (%d unique) after cutting multimapping genes" %(maps.shape[0], len(maps.index.unique().values)))

    data = data.merge(maps, left_index=True, right_index=True, how='inner').set_index('Gene')
    data = data.sort_values(by=['Mean'], ascending=False)
    data = data.loc[~data.index.duplicated(keep='first')]
    data = data.drop(columns=['Mean'])

    print("> Dataset left with %d genes after mapping genes with high mean expression probes\n" %(len(data.index.values)))

    return data


# Convert indices from/to ensembl and symbols, will drop rows with indices that don't convert

def convert_ensembl_to_symbols(dataframe):

    '''
    Convert expression matrix from ensembl gene indexed to gene symbol indexed 

    Parameters:
    ----------

    dataframe
        Dataframe containing expression values, index as ensembl ids

    Returns:
    -----------

    dataframe
        Expression dataframe with variables converted to gene symbols

    '''

    conversion = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', header=None, names=['ensembl', 'symbol'])
    conversion.drop_duplicates(subset='ensembl', keep=False, inplace=True)

    dataframe  = dataframe.merge(conversion, how='left', left_index=True, right_on='ensembl')

    dataframe.set_index(['symbol'], inplace=True)
    dataframe.drop(columns=['ensembl'], inplace=True)
    dataframe = dataframe.loc[~dataframe.index.isnull()]

    return dataframe

def convert_symbols_to_ensembl(dataframe):

    '''
    Convert expression matrix from gene symbol indexed to ensembl id index

    Parameters:
    ----------

    dataframe
        Dataframe containing expression values, index as gene symbols

    Returns:
    -----------

    dataframe
        Expression dataframe with variables converted to ensembl ids

    '''

    conversion = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', header=None, names=['ensembl', 'symbol'])
    conversion.drop_duplicates(subset='symbol', keep=False, inplace=True)

    dataframe  = dataframe.merge(conversion, how='left', left_index=True, right_on='symbol')

    dataframe.set_index(['ensembl'], inplace=True)
    dataframe.drop(columns=['symbol'], inplace=True)
    dataframe = dataframe.loc[~dataframe.index.isnull()]

    return dataframe

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

def retrieve_expression_data(dataframe):

    '''
    Makes http get requests to stemformatics api to retrieve expression data based upon datasets in the dataframe

    Parameters:
    ----------

    dataframe
        Dataframe containing annotations, requires a 'Dataset' column
    '''

    import requests, io

    expression = pd.DataFrame()

    for i_dataset in dataframe.Dataset.unique():
        url = 'https://api-dev.stemformatics.org/datasets/%d/expression?as_file=True&key=raw' %i_dataset
        req = requests.get(url, verify=False) #Lol, should we trust s4m?
        i_df = pd.read_csv(io.StringIO(req.text), sep='\t', header=0, index_col=0)
        #i_df.columns = [str(int(i_dataset))+'_'+str(i_name) for i_name in i_df.columns.values] #This ought to be the right format for s4m going forward

    return expression


