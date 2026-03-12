install.packages("Seurat")
install.packages("Matrix")
install.packages("dplyr")
install.packages("patchwork")

if(!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("GEOquery")

library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(GEOquery)

getGEOSuppFiles("GSE289940", makeDirectory = TRUE)
untar("GSE289940/GSE289940_RAW.tar", exdir = "GSE289940/raw")
list.files("GSE289940/raw", recursive = TRUE)



samples <- c(
  "GSM8801891_SWI_L1",
  "GSM8801892_SWI_L2",
  "GSM8801893_Con_L1",
  "GSM8801894_Con_L2"
)
data_dir <- "GSE289940/raw"

seurat_list <- lapply(samples, function(sample){
  
  matrix_file <- file.path(data_dir, paste0(sample, "_matrix.mtx.gz"))
  feature_file <- file.path(data_dir, paste0(sample, "_features.tsv.gz"))
  barcode_file <- file.path(data_dir, paste0(sample, "_barcodes.tsv.gz"))
  
  mat <- readMM(matrix_file)
  
  features <- read.delim(feature_file, header = FALSE)
  barcodes <- read.delim(barcode_file, header = FALSE)
  
  gene_names <- make.unique(features$V2)
  
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- sample
  
  return(obj)
})
seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[2:4],
  add.cell.ids = samples,
  project = "GSE289940"
)
seurat$condition <- ifelse(
  grepl("SWI", seurat$sample),
  "Immobilized",
  "Control"
)
dim(seurat)


seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern="^mt-")
VlnPlot(seurat,
        features=c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol=3)

seurat <- subset(
  seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)
dim(seurat)
seurat <- NormalizeData(seurat)

seurat <- FindVariableFeatures(seurat, selection.method="vst", nfeatures=2000)
VariableFeaturePlot(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
ElbowPlot(seurat)


seurat <- FindNeighbors(seurat, dims=1:20)

seurat <- FindClusters(seurat, resolution=0.5)

seurat <- RunUMAP(seurat, dims=1:20)

DimPlot(seurat, label=TRUE)

DimPlot(seurat, group.by="condition")

seurat <- JoinLayers(seurat)
markers <- FindAllMarkers(
  seurat,
  only.pos=TRUE,
  min.pct=0.25,
  logfc.threshold=0.25
)

head(markers)

library(dplyr)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(seurat, features = top10$gene)

FeaturePlot(seurat, features = "Tnfsf11")

VlnPlot(seurat, features="Tnfsf11", group.by="condition")


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

osteoblast_clusters <- c(2,4,5)
osteoclast_clusters <- c(0,6,10)


seurat$celltype <- "Other"

seurat$celltype[Idents(seurat) %in% osteoblast_clusters] <- "Osteoblast"
seurat$celltype[Idents(seurat) %in% osteoclast_clusters] <- "Osteoclast"
DimPlot(seurat, group.by="celltype", label=TRUE)
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





VlnPlot(
  seurat,
  features="Slit3",
  group.by="celltype"
)

DotPlot(
  seurat,
  features="Slit3",
  group.by="celltype"
) + RotatedAxis()



osteoblast <- subset(seurat, idents = osteoblast_clusters)
FeaturePlot(osteoblast, features="Slit3")
VlnPlot(osteoblast, features="Slit3", group.by="condition")

osteoclast <- subset(seurat, idents = osteoclast_clusters)
FeaturePlot(osteoclast, features="Slit3")
VlnPlot(osteoclast, features="Slit3", group.by="condition")

DotPlot(
  osteoblast,
  features="Slit3",
  group.by="condition"
) + RotatedAxis()

DotPlot(
  osteoclast,
  features="Slit3",
  group.by="condition"
) + RotatedAxis()





AverageExpression(
  seurat,
  features="Slit3",
  group.by="celltype"
)

FeaturePlot(
  seurat,
  features = c("Tnfsf11","Slit3","Runx2","Bglap", "Ctsk"),
  ncol = 3
)
