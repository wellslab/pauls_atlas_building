import numpy as np
import pandas as pd
import scipy, sys, sklearn.mixture, sklearn.cluster
sys.path.append('/Users/pwangel/pauls_atlas_building/notebooks/')
import functions

from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *

main_ensembl_ids = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/ensembl_hg38.91_chromosome.csv').ensembl_gene_id.values.astype(str)

dataset_list = ['/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7124/source/processed.hg38.91.20180713-201511/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7135/source/processed.hg38.91.20180716-154848/gene_count_frags.txt',
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7240/source/processed.hg38.91.20180718-124924/gene_count_frags.txt', 
                '/Users/pwangel/Data/ensembl_hg38.91/RNASeq_Data/7253/source/processed.hg38.91.20180815-171954/gene_count_frags.txt']

gene_list      = [
                  'KLF4',
                  'KLF2',
                  'STAT3',
                  'TCF3',
                  'TBX3',
                  'GBX2',
                  'SALL4',
                  'ESRRB',
                  'TFCP2L1',
                  'POU5F1', #multiple ensemble ids for this one..
                  'NANOG',
                  'SOX2',
                  'TFE3', 
                  'KLF17', 
                  #'SSEA3', why can't I find this? 
                  #'TRA-1-60', why can't I find this?
                  'CD24', 
                  #'SSEA4', why can't I find this? 
                  'NLGN4X', 
                  'F11R', 
                  'EPCAM', 
                  'OTX2', 
                  'ZIC2', 
                  'ZIC3', 
                  'ZIC5', 
                  'KLF5', 
                  'DPPA3', 
                  'MAEL',
                  'UTF1', 
                  'ZFP57',
                  'DNMT3L'
                  #'MEKERK'
                  ]


#data = pd.DataFrame(index=main_ensembl_ids)
data = pd.DataFrame()
for i_fname in dataset_list:
    data = data.merge(pd.read_csv(i_fname, sep='\t'), how='outer', left_index=True, right_index=True)

#data        = pd.read_csv('/Users/pwangel/Downloads/pluripotent_atlas_data.tsv', sep='\t', index_col=0)
annotations = pd.read_csv('/Users/pwangel/Downloads/RNASeq_pluripotent_annotations.tsv', sep='\t', index_col=0)
genes       = pd.read_csv('/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', index_col=0, names=['symbol'])

annotations.index = [i+';'+j for i,j in zip(annotations.chip_id.values.astype(str), annotations.Dataset.values.astype(int).astype(str))]
data        = data[annotations.chip_id]

data = functions.transform_to_percentile(data)

kmeans = sklearn.cluster.KMeans(n_clusters=2)

data_output = pd.DataFrame(index=gene_list, columns=['Bimodal val', 'Low_Expression', 'High_Expression', 'Low_Std', 'High_Std'])

#Loop through the sample types and genes to find differentially expressed genes
for i_gene in gene_list:

    ensembl_id = genes.index.values[genes.symbol.values==i_gene]

    if i_gene == 'POU5F1':
        ensembl_id = ['ENSG00000204531']
    if i_gene == 'ZFP57':
        ensembl_id = ['ENSG00000204644']

    
    kmeans.fit(data.loc[ensembl_id].values.reshape(-1, 1))
    prediction = kmeans.predict(data.loc[ensembl_id].values.reshape(-1, 1))

    delta_mean = data.loc[ensembl_id,prediction==0].values.mean()-data.loc[ensembl_id,prediction==1].values.mean()
    std_sum    = data.loc[ensembl_id,prediction==0].values.std()+data.loc[ensembl_id,prediction==1].values.std()

    low_exp_group = np.argmin([data.loc[ensembl_id,prediction==0].values.mean(),data.loc[ensembl_id,prediction==1].values.mean()])

    data_output.loc[i_gene, 'Bimodal val'] = np.abs(delta_mean)/std_sum
    data_output.loc[i_gene, 'Low_Expression'] = data.loc[ensembl_id,prediction==low_exp_group].values.mean()
    data_output.loc[i_gene, 'High_Expression'] = data.loc[ensembl_id,prediction!=low_exp_group].values.mean()
    data_output.loc[i_gene, 'Low_Std'] = data.loc[ensembl_id,prediction==low_exp_group].values.std()
    data_output.loc[i_gene, 'High_Std'] = data.loc[ensembl_id,prediction!=low_exp_group].values.std()

#only want genes that show some bimodality
data_output = data_output.loc[data_output['Bimodal val']>0.9]

condition_dataframe = pd.DataFrame(index=data_output.index)
for i_dataset in annotations.Dataset.unique():
    for j_condition in annotations.loc[annotations.Dataset==i_dataset]['Experiment'].unique():
        for k_celltype in annotations.loc[(annotations.Dataset==i_dataset)&(annotations.Experiment==j_condition)].LM_Group_COLOR.unique():
            i_name = '_'.join([str(i_dataset),str(j_condition),str(k_celltype)])
            sel    = (annotations.Dataset.values == i_dataset) & (annotations.Experiment.values==j_condition) & (annotations.LM_Group_COLOR.values==k_celltype)
            condition_dataframe[i_name] = np.nan

            for i_gene in data_output.index.values:

                ensembl_id = genes.index.values[genes.symbol.values==i_gene]

                if i_gene == 'POU5F1':
                    ensembl_id = ['ENSG00000204531']
                if i_gene == 'ZFP57':
                    ensembl_id = ['ENSG00000204644']

                if np.abs(data.loc[ensembl_id, sel].values.mean()-data_output.loc[i_gene].Low_Expression)/data_output.loc[i_gene].Low_Std < \
                   np.abs(data.loc[ensembl_id, sel].values.mean()-data_output.loc[i_gene].High_Expression)/data_output.loc[i_gene].High_Std:
                    condition_dataframe.loc[i_gene, i_name] = 0
                else:
                    condition_dataframe.loc[i_gene, i_name] = 1

#Open two files to place experimental conditions into, we will append one to the other afterwards
f_experiments  = open('/Users/pwangel/Downloads/temp/temp_experiments.txt', 'a') 
f_culture_cond = open('/Users/pwangel/Downloads/temp/temp_conditions.txt', 'a')

#Print out results into text files for loading into REIN
experiment_counter=0
for i_dataset in annotations.Dataset.unique():
    sel_conditions = (annotations.Dataset==i_dataset)&~pd.isnull(annotations.Experiment)&(annotations.Time>0)
    for j_condition, initial_condition in zip(annotations.loc[sel_conditions].groupby(['Experiment', 'Initial Condition']).size().reset_index()['Experiment'].values, \
                                              annotations.loc[sel_conditions].groupby(['Experiment', 'Initial Condition']).size().reset_index()['Initial Condition'].values):

        print(j_condition)
        f_experiments.write('#Experiment%d[20] |= $%s "";\n' %(experiment_counter,j_condition))
        if initial_condition!='None':
            f_experiments.write('#Experiment%d[0] |= $%s "";\n' %(experiment_counter, initial_condition))
        f_experiments.write('\n')

        f_culture_cond.write('$%s :=\n' %j_condition)
        f_culture_cond.write('{\n')
        for i_gene in condition_dataframe.index.values:
            f_culture_cond.write(' %s = %d\n' %(i_gene, condition_dataframe.loc[i_gene,'_'.join([str(i_dataset),str(j_condition),str('naive')])]))
        f_culture_cond.write('};\n')
        f_culture_cond.write('\n')

        experiment_counter+=1

for initial_condition in annotations['Initial Condition'].dropna().unique():

    if initial_condition!='None':
        f_culture_cond.write('$%s :=\n' %(initial_condition))
        f_culture_cond.write('{\n')
        for i_gene in condition_dataframe.index.values:
            annotations_row = annotations.loc[annotations['Experiment']==initial_condition]
            f_culture_cond.write(' %s = %d\n' %(i_gene, condition_dataframe.loc[i_gene,'_'.join([str(i_dataset),str(annotations_row.Experiment.unique()[0]),str(annotations_row.LM_Group_COLOR.unique()[0])])]))
        f_culture_cond.write('};\n')
        f_culture_cond.write('\n')

f_experiments.write('\n')

f_experiments.close()
f_culture_cond.close()

f = open("/Users/pwangel/Downloads/temp/observations.txt", "w")
for tempfile in ['/Users/pwangel/Downloads/temp/temp_experiments.txt', '/Users/pwangel/Downloads/temp/temp_conditions.txt']:
    f.write(open(tempfile).read())
f.close()

f_model = open('/Users/pwangel/Downloads/temp/model.txt', 'a')

f_model.write('directive updates sync;\n')
f_model.write('directive length 20;\n')
f_model.write('directive uniqueness interactions;\n')
f_model.write('directive limit 0;\n')
f_model.write('directive regulation legacy;\n')
'[](0..15); '.join([str(i_gene) for i_gene in gene_list])+';\n'
for i_gene in gene_list:
    for j_gene in gene_list:
        direction = np.random.choice(['positive', 'negative', 'positive optional', 'negative optional'])
        f_model.write('%s\t%s\t%s;\n' %(i_gene, j_gene, direction))

f_model.close()
