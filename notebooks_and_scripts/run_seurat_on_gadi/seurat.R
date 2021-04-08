library(Seurat)
library(ggplot2)

sc.data <- Read10X(data.dir = "/scratch/yf0/pa5933/Data/Single_Cell/SC_data_pack_pt1",gene.column = 2,unique.features = TRUE,strip.suffix = FALSE)

cat("Read in data")

platforms <- read.csv("/scratch/yf0/pa5933/Data/Single_Cell/SC_data_pack_pt1/metadata.tsv", sep = "\t", header = TRUE, row.names=1)

cat("Read in metadata")

sc.seurat <- CreateSeuratObject(counts = sc.data)

cat("Made seurat object")

sc.seurat <- AddMetaData(object = sc.seurat, metadata = platforms)

cat("Added metadata")

sc.list <- SplitObject(sc.seurat, split.by = "Dataset")

#sc.list <- sc.list[c("HG-U133_Plus_2", "HuGene", "Illumina V2", "Illumina V4", "RNASeq")] 
sc.list <- sc.list[unique(platforms[["Dataset"]])] 

for (i in 1:length(sc.list)) {
    sc.list[[i]] <- NormalizeData(sc.list[[i]], verbose = FALSE)
    sc.list[[i]] <- FindVariableFeatures(sc.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

#reference.list <- sc.list[c("HG-U133_Plus_2", "HuGene", "Illumina V4", "RNASeq")]
reference.list <- sc.list[unique(platforms[["Dataset"]])] 

cat("Beginning to find anchors")

sc.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, k.filter = 100)

cat("Beginning to integrate using only anchors")

sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:30)
integrated_data <- GetAssayData(object = sc.integrated, assay = "integrated", slot = "data")

nrow(sc.integrated$integrated)
ncol(sc.integrated$integrated)

write.table(integrated_data, file="/scratch/yf0/pa5933/Data/Single_Cell/SC_data_pack_pt1/integrated_anchors.tsv", quote=FALSE, sep='\t')

cat("Beginning to integrate using all features")

sc.integrated <- IntegrateData(anchorset = sc.anchors, dims = 1:30, features.to.integrate = rownames(sc.seurat))
integrated_data <- GetAssayData(object = sc.integrated, assay = "integrated", slot = "data")

nrow(sc.integrated$integrated)
ncol(sc.integrated$integrated)

#Takes up too much space...not a sparse matrix anymore...
#write.table(integrated_data, file="/scratch/yf0/pa5933/Data/Single_Cell/SC_data_pack_pt1/integrated_data.tsv", quote=FALSE, sep='\t')

