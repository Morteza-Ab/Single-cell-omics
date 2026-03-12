library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(GEOquery)

getGEOSuppFiles("GSE221363", makeDirectory = TRUE)
untar("GSE221363/GSE221363_RAW.tar", exdir = "GSE221363/raw")
list.files("GSE221363/raw", recursive = TRUE)


data_dir <- "GSE221363/raw"
sample <- "GSM6859099"

matrix_file <- file.path(data_dir, paste0(sample, "_matrix.mtx.gz"))
feature_file <- file.path(data_dir, paste0(sample, "_genes.tsv.gz"))
barcode_file <- file.path(data_dir, paste0(sample, "_barcodes.tsv.gz"))

mat <- readMM(matrix_file)

features <- read.delim(feature_file, header = FALSE)
barcodes <- read.delim(barcode_file, header = FALSE)

gene_names <- make.unique(features$V2)
rownames(mat) <- gene_names
colnames(mat) <- barcodes$V1

library(Seurat)
seurat <- CreateSeuratObject(
  counts = mat,
  project = sample,
  min.cells = 3,
  min.features = 200
)


### QC and filter


seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern="^mt-")

VlnPlot(seurat,
        features=c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol=3)

# Filter cells
seurat <- subset(
  seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)
dim(seurat)

### Normalize, find variable features, scale
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method="vst", nfeatures=2000)
VariableFeaturePlot(seurat)
seurat <- ScaleData(seurat)



### PCA, neighbors, clustering, UMAP
seurat <- RunPCA(seurat)
ElbowPlot(seurat)

seurat <- FindNeighbors(seurat, dims=1:20)
seurat <- FindClusters(seurat, resolution=0.5)
seurat <- RunUMAP(seurat, dims=1:20)

DimPlot(seurat, label=TRUE)



FeaturePlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk"),
  ncol = 3
)

DotPlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk")
) + RotatedAxis()

table(Idents(seurat))



osteoblast_clusters <- c(15,15)
osteoclast_clusters <- c(7,9)

seurat$celltype <- "Other"
seurat$celltype[Idents(seurat) %in% osteoblast_clusters] <- "Osteoblast"
seurat$celltype[Idents(seurat) %in% osteoclast_clusters] <- "Osteoclast"

DimPlot(seurat, group.by="celltype", label=TRUE)


grep("Slit3", rownames(seurat), value = TRUE) ### slit3 exist in the data set (use capital first letter in mouse dataset for genes)
FeaturePlot(seurat, features="Slit3")

VlnPlot(seurat, features="Slit3")
DotPlot(seurat, features="Slit3") + RotatedAxis()

VlnPlot(
  seurat,
  features="Slit3",
  group.by="seurat_clusters"
)









### with 6 pca based on elbow
seurat <- FindNeighbors(seurat, dims = 1:6)
seurat <- FindClusters(seurat, resolution = 0.4)
seurat <- RunUMAP(seurat, dims = 1:6)
DimPlot(seurat, label=TRUE)

FeaturePlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk"),
  ncol = 3
)

DotPlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk")
) + RotatedAxis()

table(Idents(seurat))




osteoblast_clusters <- c(15)
osteoclast_clusters <- c(9, 5)

seurat$celltype <- "Other"
seurat$celltype[Idents(seurat) %in% osteoblast_clusters] <- "Osteoblast"
seurat$celltype[Idents(seurat) %in% osteoclast_clusters] <- "Osteoclast"

DimPlot(seurat, group.by="celltype", label=TRUE)


grep("Slit3", rownames(seurat), value = TRUE) ### slit3 exist in the data set (use capital first letter in mouse dataset for genes)
FeaturePlot(seurat, features="Slit3")

VlnPlot(seurat, features="Slit3")
DotPlot(seurat, features="Slit3") + RotatedAxis()

VlnPlot(
  seurat,
  features="Slit3",
  group.by="seurat_clusters"
)











