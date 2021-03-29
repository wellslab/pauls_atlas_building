import numpy as np
import pandas as pd
import scipy, sys, sklearn, sklearn.mixture, sklearn.cluster
sys.path.append('/Users/pwangel/pauls_atlas_building/notebooks/')
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

main_ensembl_ids = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/ensembl_hg38.91_chromosome.csv').ensembl_gene_id.values.astype(str)

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

dataset_list = [ #Get rid of this when the data is better curated...
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6884/source/processed.hg38.91.20180627-112849/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7046/source/processed.hg38.91.20180627-121841/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6608/source/processed.hg38.91.20181010-180343/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7264/source/processed.hg38.91.20180807-153815/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6639/source/processed.hg38.91.20180815-162300/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6730/source/processed.hg38.91.20180807-163949/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6896/source/processed.hg38.91.20180528-101935/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7061/source/processed.hg38.91.20180717-162353/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/6855/source/processed.hg38.91.20180628-124614/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7124/source/processed.hg38.91.20180713-201511/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7159/source/processed.hg38.91.20180718-154044/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7188/source/processed.hg38.91.20180716-213651/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7240/source/processed.hg38.91.20180718-124924/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7239/source/processed.hg38.91.20181123-134435/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7268/source/processed.hg38.91.20180815-170348/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7171/source/processed.hg38.91.20180727-152501/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7135/source/processed.hg38.91.20180716-154848/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7339/source/processed.hg38.91.20181008-172616/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7228/source/processed.hg38.91.20180806-170827/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7350/source/processed.hg38.91.20181102-183434.TRIMMED/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7249/source/processed.hg38.91.20180920-124838/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7254/source/processed.hg38.91.20180614-122055.MAIN/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7275/source/processed.hg38.91.20180905-124455/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7253/source/processed.hg38.91.20180815-171954/gene_count_frags.txt'
                #'/Users/pwangel/Data/External_Data/GSE114873/GSE114873_Cufflinks_output_FPKM_and_Gene_Names.txt', 
                #'/Users/pwangel/Data/External_Data/GSE123055/GSE123055_counts.txt'
                ]

genes_df = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/REIN_genes.tsv', sep='\t', index_col=0)
genes_network = pd.read_csv('/Users/pwangel/Downloads/REIN_network.txt', sep='\t', skiprows=5, header=0, names=['Gene A','Gene B', 'Relation'])

#data = pd.DataFrame(index=main_ensembl_ids)
data = pd.DataFrame()
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='outer', left_index=True, right_index=True)

#some external data
df1 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE114873/GSE114873_Cufflinks_output_FPKM_and_Gene_Names.txt', sep='\t', index_col=0)[['condition', 'replicate', 'raw_frags']]
df1['sample_id'] = [i+'_'+j for i,j in zip(df1.condition.to_list(),df1.replicate.to_list())]
df1.drop(labels=['condition', 'replicate'], axis=1, inplace=True)
df1 = df1.pivot(columns='sample_id')
df1.columns = df1.columns.droplevel()
df1 = convert_symbols_to_ensembl(df1)

df2 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE123055/GSE123055_counts.txt', sep='\t', index_col=0)
df2.drop(labels=['Unnamed: 52'], axis=1, inplace=True)

df3 = pd.read_csv('/Users/pwangel/Data/External_Data/GSE137295/GSE137295_All_Raw_CountData.csv', index_col=0)

#### THS IS TEMPORARY UNTIL NEW DATASETS ARE PROCESSED!!!!
data = data.merge(df1, how='left', left_index=True, right_index=True)
data = data.merge(df2, how='left', left_index=True, right_index=True)
data = data.merge(df3, how='left', left_index=True, right_index=True)
data.fillna(0.0, inplace=True)

#data        = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/pluripotent_RNASeq_annotations.tsv', sep='\t', index_col=0)
lizzi_anno  = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/stemcell_annotations.tsv', sep='\t', index_col=0)
annotations = annotations.merge(lizzi_anno['LM_Group_COLOR'], how='left', left_index=True, right_index=True)
experiment_anno = pd.read_csv('/Users/pwangel/Downloads/RNASeq_only_pluripotent_annotations.tsv', sep='\t', index_col=0)
experiment_anno.index = [i+';'+j for i,j in zip(experiment_anno.chip_id.values.astype(str), experiment_anno.Dataset.values.astype(int).astype(str))]
annotations = annotations.merge(experiment_anno[['Experiment', 'Time', 'Initial Condition']], how='left', left_index=True, right_index=True)
#annotations.Dataset = annotations.Dataset.astype(float).astype(int).astype(str)

genes       = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])
gene_list   = np.intersect1d(genes.loc[np.intersect1d(data.index, genes.index)].symbol.values, genes_df.index.values)
annotations['chip_id'] = [i.split(';')[0] for i in annotations.index.values.astype(str)]
annotations = annotations.loc[(annotations.Platform_Category=='RNASeq') & (annotations.Dataset!='7275.0')] #This dataset is mistakenly in, it is annotated endoderm
data        = data[annotations.chip_id]

data = functions.transform_to_percentile(data)

# Run pca

pca        = sklearn.decomposition.PCA(n_components=10, svd_solver='full')
pca.fit(functions.transform_to_percentile(data.transpose()))
pca_coords = pca.transform(data.transpose())

functions.plot_pca(pca_coords, annotations,pca, labels=['generic_sample_type', 'Platform_Category', 'Dataset'], colour_dict={}, \
                   pcs=[1,2,3], out_file='/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/RNASeq_only_pluripotent.html')

#### Apply k means clustering to divide genes on/off state

kmeans = sklearn.cluster.KMeans(n_clusters=2)
#data_output = pd.DataFrame(index=gene_list, columns=['Bimodal val', 'Low_Expr', 'High_Expr', 'Low_Std', 'High_Std'])

for i_gene in gene_list:

    ensembl_id = genes_df.loc[genes_df.index.values==i_gene].Ensembl

    kmeans.fit(data.loc[ensembl_id].values.reshape(-1, 1))
    prediction = kmeans.predict(data.loc[ensembl_id].values.reshape(-1, 1))

    delta_mean = data.loc[ensembl_id,prediction==0].values.mean()-data.loc[ensembl_id,prediction==1].values.mean()
    std_sum    = data.loc[ensembl_id,prediction==0].values.std()+data.loc[ensembl_id,prediction==1].values.std()

    low_exp_group = np.argmin([data.loc[ensembl_id,prediction==0].values.mean(),data.loc[ensembl_id,prediction==1].values.mean()])

    if pd.isnull(genes_df.loc[i_gene, 'Lower Threshold']):
        genes_df.loc[i_gene, 'Lower Threshold'] = np.round(data.loc[ensembl_id,prediction==low_exp_group].values.mean()+data.loc[ensembl_id,prediction==low_exp_group].values.std(), decimals=2)
    if pd.isnull(genes_df.loc[i_gene, 'Upper Threshold']):
        genes_df.loc[i_gene, 'Upper Threshold'] = np.round(data.loc[ensembl_id,prediction!=low_exp_group].values.mean()-data.loc[ensembl_id,prediction!=low_exp_group].values.std(), decimals=2)
    genes_df.loc[i_gene, 'Bimodal Score']   = np.round(((data.loc[ensembl_id].values<genes_df.loc[i_gene, 'Lower Threshold']) | \
                                             (data.loc[ensembl_id].values>genes_df.loc[i_gene, 'Upper Threshold'])).sum()/data.shape[1], decimals=2)

    prediction[:] = 0
    prediction[(data.loc[ensembl_id].values<genes_df.loc[i_gene, 'Lower Threshold']).flatten()] = -1
    prediction[(data.loc[ensembl_id].values>genes_df.loc[i_gene, 'Upper Threshold']).flatten()] = 1

    # Make a plot for each gene so we can examine the bimodality
    functions.plot_dot_plots(data, annotations.merge(pd.DataFrame(index=annotations.index, data=prediction, columns=['Prediction']), how='inner', left_index=True, right_index=True), \
              cell_property='Dataset', cell_colour_by='Prediction', gene_list=ensembl_id, output_dir='/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/', plot_hist=True)

genes_df.to_csv('/Users/pwangel/PlotlyWorkspace/combine_data/naive_stemcells/REIN_genes.tsv', sep='\t')
#only want genes that show some bimodality
#genes_df = genes_df.loc[genes_df['Bimodal Score']<0.2]

condition_dataframe = pd.DataFrame(index=genes_df.index)
for i_dataset in annotations.Dataset.unique():
    for j_condition in annotations.loc[annotations.Dataset==i_dataset]['Experiment'].unique():
        for k_celltype in annotations.loc[(annotations.Dataset==i_dataset)&(annotations.Experiment==j_condition)].LM_Group_COLOR.unique():
            i_name = '_'.join([str(i_dataset),str(j_condition),str(k_celltype)])
            print(i_name)
            sel    = (annotations.Dataset.values == i_dataset) & (annotations.Experiment.values==j_condition) & (annotations.LM_Group_COLOR.values==k_celltype)
            condition_dataframe[i_name] = np.nan

            for i_gene in genes_df.index.values:

                ensembl_id = genes_df.loc[genes_df.index.values==i_gene].Ensembl

                if data.loc[ensembl_id, sel].values.mean()<genes_df.loc[i_gene, 'Lower Threshold']:
                    condition_dataframe.loc[i_gene, i_name] = 0
                elif data.loc[ensembl_id, sel].values.mean()>genes_df.loc[i_gene, 'Upper Threshold']:
                    condition_dataframe.loc[i_gene, i_name] = 1
                else:
                    condition_dataframe.loc[i_gene, i_name] = np.nan

genes_network = genes_network.loc[np.in1d(genes_network['Gene A'].values, genes_df.loc[genes_df.Inclusion.values].index) & \
                                  np.in1d(genes_network['Gene B'].values, genes_df.loc[genes_df.Inclusion.values].index)]
genes_network.to_csv('/Users/pwangel/Downloads/temp/temp.REIN', sep='\t', index=False)

f_experiments  = open('/Users/pwangel/Downloads/temp/temp.REIN', 'a') 
f_culture_cond = open('/Users/pwangel/Downloads/temp/temp_conditions.txt', 'a')
f_experiments.write('\n')

#Print out results into text files for loading into REIN
experiment_counter=0
for i_dataset in annotations.Dataset.unique():
    sel_conditions = (annotations.Dataset==i_dataset)&~annotations.Experiment.isnull()&(annotations.Time>0)
    for j_condition, initial_condition in zip(annotations.loc[sel_conditions].groupby(['Experiment', 'Initial Condition'], dropna=False).size().reset_index()['Experiment'].values, \
                                              annotations.loc[sel_conditions].groupby(['Experiment', 'Initial Condition'], dropna=False).size().reset_index()['Initial Condition'].values):

        print(j_condition)
        print(annotations.loc[sel_conditions].groupby(['Experiment', 'Initial Condition'], dropna=False).size().reset_index())
        # These are the name as appearing in the condition matrix
        condition_name = '_'.join([str(i_dataset),str(j_condition),str('naive')])

        f_experiments.write('#Experiment%d[20] |= $Experiment%d_%s "";\n' %(experiment_counter,experiment_counter, j_condition))
        if initial_condition!='None':
            annotations_row = annotations.loc[(annotations.Dataset==i_dataset)&(annotations['Experiment']==initial_condition)]
            initial_condition_name = '_'.join([str(i_dataset),str(initial_condition),str(annotations_row.LM_Group_COLOR.unique()[0])])
            f_experiments.write('#Experiment%d[0] |= $Experiment%d_%s "";\n' %(experiment_counter,experiment_counter, initial_condition))
        f_experiments.write('\n')

        f_culture_cond.write('$Experiment%d_%s :=\n' %(experiment_counter, j_condition))
        f_culture_cond.write('{\n')
        for i_gene in condition_dataframe.loc[~condition_dataframe.loc[:,condition_name].isnull()].index.values:
            if i_gene!=condition_dataframe.loc[~condition_dataframe.loc[:,condition_name].isnull()].index.values[-1]:
                f_culture_cond.write(' %s = %d and\n' %(i_gene, condition_dataframe.loc[i_gene,condition_name]))
            else:
                f_culture_cond.write(' %s = %d\n' %(i_gene, condition_dataframe.loc[i_gene,condition_name]))
        f_culture_cond.write('};\n')
        f_culture_cond.write('\n')

        if initial_condition!='None':
            f_culture_cond.write('$Experiment%d_%s :=\n' %(experiment_counter, initial_condition))
            f_culture_cond.write('{\n')
            for i_gene in condition_dataframe.loc[~condition_dataframe.loc[:,initial_condition_name].isnull()].index.values:
                if i_gene!=condition_dataframe.loc[~condition_dataframe.loc[:,initial_condition_name].isnull(),initial_condition_name].index.values[-1]:
                    f_culture_cond.write(' %s = %d and\n' %(i_gene, condition_dataframe.loc[i_gene,initial_condition_name]))
                else:
                    f_culture_cond.write(' %s = %d\n' %(i_gene, condition_dataframe.loc[i_gene,initial_condition_name]))
            f_culture_cond.write('};\n')
            f_culture_cond.write('\n')

        experiment_counter+=1

f_experiments.write('\n')

f_experiments.close()
f_culture_cond.close()


f = open("/Users/pwangel/Downloads/temp/temp.REIN", "a")
for tempfile in ['/Users/pwangel/Downloads/temp/temp_conditions.txt']:
    f.write(open(tempfile).read())
f.close()

