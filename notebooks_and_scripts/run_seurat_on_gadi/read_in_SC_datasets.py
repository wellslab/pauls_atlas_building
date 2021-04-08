import numpy, scipy, sys, os, pandas, sklearn, sklearn.decomposition, scanpy, gc

#Set parameters for cell types and output dataframes
min_library_size = 5000
min_genes        = 1000

df_metadata  = pandas.DataFrame()
df_data_out = pandas.DataFrame()

data_locations = '/scratch/yf0/pa5933/Data'

#Madissoon et al data

folder = data_location+'/Single_Cell/Madissoon/'
fname  = data_location+'/Single_Cell/Madissoon/matrix.mtx'

data  = scanpy.read_mtx(fname)
data = scanpy.AnnData(X=data.X.transpose())

samples = pandas.read_csv(folder+"/cells.tsv", skiprows=0, sep='\t')
genes = pandas.read_csv(folder+"/features.tsv", skiprows=0, header=None, sep='\t', index_col=0)
genes = genes.index.values

print("Raw Madissoon data shape",data.shape,samples.shape,genes.shape)
print("Using classified data")

cls_samples = pandas.read_csv(folder+'/tissue_sensitivity_hca_dcp_uuid_mapping_per_cell.csv', sep=',')
cls_samples.index = numpy.core.defchararray.add(cls_samples.barcode.values.astype(str), cls_samples['cell_suspension.provenance.document_id'].values.astype(str))
print(cls_samples.shape, numpy.unique(cls_samples.index.values).shape)

samples.index = numpy.core.defchararray.add(samples.barcode.values.astype(str), samples['cell_suspension.provenance.document_id'].values.astype(str))
samples = samples.merge(cls_samples, how='left', left_index=True, right_index=True)

#Remove samples with low library size and no cell type classification
sel     = (numpy.array(data.X.sum(1)).reshape(-1)>=min_library_size) & (~samples.Celltypes.isnull().values) & (numpy.array(data.X.astype(bool).sum(axis=1)).reshape(-1) >= min_genes) & (numpy.core.defchararray.count(samples.Celltypes.values.astype(str), 'Epi')==0) & \
        ~numpy.in1d(samples.Celltypes.values.astype(str), ['Glands_duct','Alveolar_Type1','Alveolar_Type2','Blood_vessel','Ciliated','Fibroblast','Lymph_vessel','Muscle_cells'])
samples = samples.loc[sel]
data    = scanpy.AnnData(X=data.X[sel,:])

print("Found %d samples above library size=%d" %(sel.sum(), min_library_size))

sys.stdout.flush()

for i_donor in numpy.unique(samples['Donor'].values.astype(str)):
    for i_tissue in numpy.unique(samples['organ'].values.astype(str)):

        i_sel = (samples['Donor'].values.astype(str)==i_donor) & (samples['organ'].values.astype(str)==i_tissue)

        if (i_sel.sum()>100) and (i_tissue!='nan'):

            print("Analysing data from donor %s and tissue %s" %(i_donor, i_tissue))

            i_data = scanpy.AnnData(X=data.X[i_sel,:]) 
            i_samples = samples.loc[i_sel]

            i_data = i_data.to_df().transpose()
            i_data.index=genes

            i_name = 'donor_'+str(i_donor)+'_'+i_tissue.replace(' ', '_')

            i_metadata = pandas.DataFrame(index=i_data.columns)
            i_metadata['celltype'] = i_samples.Celltypes.values
            i_metadata['Platform'] = '10X v2'
            i_metadata['Dataset'] = 'Madissoon'
            i_metadata['Donor'] = str(i_donor)
            i_metadata['Tissue'] = str(i_tissue)
            i_metadata['Disease'] = 'Healthy'

            df_metadata = df_metadata.merge(i_metadata.transpose(), how='outer', left_index=True, right_index=True)
            df_data_out = df_data_out.merge(i_data, how='outer', left_index=True, right_index=True)

            print(i_data.shape, i_metadata.shape)

            sys.stdout.flush()

del data
#Ding et al data

folder = data_location+'/Single_Cell/Ding/'
fname  = data_location+'/Single_Cell/Ding/counts.read.txt'

#Read in and form clusters
data = scanpy.read_mtx(fname)
data = scanpy.AnnData(X=data.X.transpose())

samples = pandas.read_csv(folder+"/cells.read.txt", skiprows=0, header=None, sep='\t')
genes = pandas.read_csv(folder+"/genes.read.txt", skiprows=0, header=None, sep='\t')
genes = [item.split("_")[0] for item in genes[0].values.astype(str)]

name_map = pandas.read_csv(folder+"/map.CPM.names.Count.names.txt", sep='\t', index_col=0)
samples_meta = pandas.read_csv(folder+"/meta.txt", sep="\t", index_col=0)

name_map = name_map.merge(samples_meta, how='left', left_index=True, right_index=True)
samples  = samples.merge(name_map, how='left', left_on=0, right_index=True).set_index(0)

#The authors do not provide meta data for all samples, presumably they are excluded from the analysis for a good reason
sel  = ~samples.Method.isnull().values & ~samples.CellType.isnull().values & (numpy.array(data.X.sum(1)).reshape(-1)>=min_library_size) \
    & (numpy.array(data.X.astype(bool).sum(axis=1)).reshape(-1) >= min_genes)
data    = scanpy.AnnData(X=data.X[sel,:])

for i_platform in numpy.unique(samples.Method.values[sel]):

    print("Analysing data from platform %s" %str(i_platform))

    sel_platform = samples.Method.values[sel]==i_platform
    i_data = scanpy.AnnData(X=data.X[sel_platform,:])
    i_samples  = samples.loc[sel].loc[sel_platform]
    i_platform = i_platform.replace(' ', '_')
            
    i_data = i_data.to_df().transpose()
    i_data.index=genes
    df_data_out = df_data_out.merge(i_data, how='outer', left_index=True, right_index=True)

    i_metadata = pandas.DataFrame(index=i_data.columns)
    i_metadata['Platform'] = str(i_platform)
    i_metadata['Dataset'] = 'Ding'
    i_metadata['Platform_Symbol'] = 'cross'
    i_metadata['celltype'] = i_samples.CellType.values
    df_metadata = df_metadata.merge(i_metadata.transpose(), how='outer', left_index=True, right_index=True)

    sys.stdout.flush()

del data

df_metadata = df_metadata.transpose()
df_metadata['Platform'] = df_metadata.Platform.replace(to_replace='10x_Chromium_(v3)', value='10X v3').values
df_metadata['Platform'] = df_metadata.Platform.replace(to_replace='10x_Chromium_(v2)', value='10X v2').values
df_metadata['Platform'] = df_metadata.Platform.replace(to_replace='10x_Chromium_(v2)_A', value='10X v2').values
df_metadata['Platform'] = df_metadata.Platform.replace(to_replace='10x_Chromium_(v2)_B', value='10X v2').values
df_metadata['Platform_Category'] = df_metadata.Platform.values
df_metadata['Tissue'] = 'PBMC'
df_metadata = df_metadata.transpose()

print(df_metadata.shape, df_data_out.shape, df_data_out.head())

#pbmc data from 10X website

folder = data_location+'/Single_Cell/PBMC/'

PBMC_10k_fname = folder+'pbmc_10k_v3_raw_feature_bc_matrix.h5'
#PBMC_5k_fname = folder+'5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5'
PBMC_fname = PBMC_10k_fname

#Read in and form clusters
pbmc_data    = scanpy.read_10x_h5(PBMC_fname)
print(pbmc_data.shape)
scanpy.pp.filter_cells(pbmc_data, min_counts=min_library_size, inplace=True)
scanpy.pp.filter_cells(pbmc_data, min_genes=min_genes, inplace=True)
print(pbmc_data.shape)

i_data       = pbmc_data.to_df().transpose()
i_data.index = pbmc_data.var['gene_ids'].values

i_metadata = pandas.DataFrame(index=i_data.columns)
i_metadata['Platform'] = '10X v3'
i_metadata['Platform_Category'] = '10X v3'
i_metadata['celltype'] = 'Unannotated'
i_metadata['Tissue'] = 'PBMC'
i_metadata['Dataset'] = 'Chromium 10X v3 PBMC'
i_metadata['Disease'] = 'Healthy'

df_data_out  = df_data_out.merge(i_data, how='outer', left_index=True, right_index=True)
df_metadata  = df_metadata.merge(i_metadata.transpose(), how='outer', left_index=True, right_index=True)

print(df_data_out.shape, df_metadata.shape)

del pbmc_data

# Park et al pbmc influenza and covid data

folder = data_location+'/Single_Cell/Park/'
fname  = data_location+'/Single_Cell/Park/GSE149689_matrix.mtx'

data  = scanpy.read_mtx(fname)
data = scanpy.AnnData(X=data.X.transpose())

samples = pandas.read_csv(folder+"/GSE149689_barcodes.tsv", skiprows=0, header=None, sep='\t')
donors = numpy.array([item.split("-")[1] for item in samples[0].values.astype(str)])
donor_list = ['Sample_1_nCoV_1','Sample_2_nCoV_2','Sample_3_Flu_1','Sample_4_Flu_2','Sample_5_Normal_1',
              'Sample_6_Flu_3','Sample_7_Flu_4','Sample_8_Flu_5','Sample_9_nCoV_3','Sample_10_nCoV_4',
              'Sample_11_nCoV_5','Sample_12_nCoV_6','Sample_13_Normal_2','Sample_14_Normal_3',
              'Sample_15_nCoV_7','Sample_16_nCoV_8','Sample_17_nCoV_9','Sample_18_nCoV_10','Sample_19_Normal_4',
              'Sample_20_nCoV_11']

genes = pandas.read_csv(folder+"/GSE149689_features.tsv", skiprows=0, header=None, sep='\t')
genes = genes[0].values.astype(str)

sel = (numpy.array(data.X.sum(1)).reshape(-1)>=min_library_size) & (numpy.array(data.X.astype(bool).sum(axis=1)).reshape(-1) >= min_genes)
data    = scanpy.AnnData(X=data.X[sel,:])

for i_donor in numpy.unique(donors[sel])[numpy.unique(donors[sel], return_counts=True)[1]>500]: ##Three of these donors have very few samples (100-300 or so)

    print("Analysing data from donor %s" %str(i_donor))

    sel_donor = donors[sel]==i_donor
    i_data = scanpy.AnnData(X=data.X[sel_donor,:])
    i_data = i_data.to_df().transpose()
    i_data.index=genes

    i_name = donor_list[int(i_donor)-1]
    df_data_out = df_data_out.merge(i_data, how='outer', left_index=True, right_index=True)

    i_metadata = pandas.DataFrame(index=i_data.columns)
    i_metadata['Platform'] = '10X v3'
    i_metadata['Dataset'] = 'Park'
    i_metadata['Platform_Category'] = '10X v3'
    i_metadata['Platform_Symbol'] = 'cross'
    i_metadata['display_metadata'] = i_metadata.index.values.astype(str)
    i_metadata['Disease'] = i_name.split('_')[2]
    df_metadata = df_metadata.merge(i_metadata.transpose(), how='outer', left_index=True, right_index=True)

    sys.stdout.flush()

del data

df_data_out.fillna(0,inplace=True)
df_metadata = df_metadata.transpose()
df_metadata['ID'] = df_metadata.index.values
df_metadata.index = numpy.arange(df_metadata.shape[0])
df_data_out.columns = df_metadata.index

print(df_data_out.shape, df_metadata.shape, df_data_out.head())
sys.stdout.flush()

for i_folder in ['SC_data_pack', 'SC_data_pack_pt1', 'SC_data_pack_pt5']:
    if not os.path.isdir(data_location+"/Single_Cell/"+i_folder):
        os.makedirs(data_location+"/Single_Cell/"+i_folder)

(df_data_out.index.to_series()).to_csv(data_location+"/Single_Cell/SC_data_pack/genes.tsv", sep='\t', header=False)
(df_data_out.index.to_series()).to_csv(data_location+"/Single_Cell/SC_data_pack_pt1/genes.tsv", sep='\t', header=False)
(df_data_out.index.to_series()).to_csv(data_location+"/Single_Cell/SC_data_pack_pt5/genes.tsv", sep='\t', header=False)

df_metadata.to_csv(data_location+"/Single_Cell/SC_data_pack/metadata.tsv", sep='\t')
df_metadata.to_csv(data_location+"/Single_Cell/SC_data_pack/barcodes.tsv", sep='\t', columns=[], header=False)
scipy.io.mmwrite(data_location+"/Single_Cell/SC_data_pack/matrix.mtx", scipy.sparse.csr_matrix(df_data_out.values))

sel_pt5 = (numpy.random.random(size=df_metadata.shape[0])>0.5)
df_metadata.loc[sel_pt5].to_csv(data_location+"/Single_Cell/SC_data_pack_pt5/metadata.tsv", sep='\t')
df_metadata.loc[sel_pt5].to_csv(data_location+"/Single_Cell/SC_data_pack_pt5/barcodes.tsv", sep='\t', columns=[], header=False)
scipy.io.mmwrite(data_location+"/Single_Cell/SC_data_pack_pt5/matrix.mtx", scipy.sparse.csr_matrix(df_data_out.values[:,sel_pt5]))

sel_pt1 = (numpy.random.random(size=df_metadata.shape[0])>0.9)
df_metadata.loc[sel_pt1].to_csv(data_location+"/Single_Cell/SC_data_pack_pt1/metadata.tsv", sep='\t')
scipy.io.mmwrite(data_location+"/Single_Cell/SC_data_pack_pt1/matrix.mtx", scipy.sparse.csr_matrix(df_data_out.values[:,sel_pt1]))
df_metadata.loc[sel_pt1].to_csv(data_location+"/Single_Cell/SC_data_pack_pt1/barcodes.tsv", sep='\t', columns=[], header=False)
