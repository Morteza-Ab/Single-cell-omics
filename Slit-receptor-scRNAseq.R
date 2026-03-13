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
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk","Sost", "Dmp1"),
  ncol = 3
)

DotPlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk","Sost", "Dmp1")
) + RotatedAxis()

table(Idents(seurat))

osteoblast_clusters <- c(2,4,6)
osteoclast_clusters <- c(12)
osteocytes_cluster <- c(11)


seurat$celltype <- "Other"

seurat$celltype[Idents(seurat) %in% osteoblast_clusters] <- "Osteoblast"
seurat$celltype[Idents(seurat) %in% osteoclast_clusters] <- "Osteoclast"
seurat$celltype[Idents(seurat) %in% osteocytes_cluster] <- "Osteocytes"
DimPlot(seurat, group.by="celltype", label=TRUE)
DimPlot(seurat, group.by="celltype", label=TRUE)

grep("Slit3", rownames(seurat), value = TRUE)
grep("^Robo[1-4]$", rownames(seurat), value = TRUE)
#rownames(seurat)[rownames(seurat) %in% c("Robo1","Robo2","Robo3","Robo4")]

### slit3 exist in the data set (use capital first letter in mouse dataset for genes)
FeaturePlot(seurat, features="Slit3")
FeaturePlot(seurat, features=c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"))
VlnPlot(seurat, features="Slit3")
DotPlot(seurat, features="Slit3") + RotatedAxis()


DotPlot(seurat, features = c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"),
        group.by="seurat_clusters") + RotatedAxis()



VlnPlot(
  seurat,
  features = c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"),
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

DotPlot(
  seurat,
  features="Robo4",
  group.by="celltype"
) + RotatedAxis()



DotPlot(
  seurat,
  features = c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"),
  group.by = "celltype"
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







######## cluster 13 genes


cluster13_markers <- markers %>%
  filter(cluster == 13) %>%
  arrange(desc(avg_log2FC))

head(cluster13_markers, 20)
FeaturePlot(seurat, features = head(cluster13_markers$gene, 6))



#### results 
head(cluster13_markers, 20)
p_val avg_log2FC pct.1 pct.2     p_val_adj cluster     gene
Robo4     0.000000e+00  14.459064 0.854 0.000  0.000000e+00      13    Robo4
Zfp366    0.000000e+00  13.424998 0.562 0.000  0.000000e+00      13   Zfp366
Aplnr     0.000000e+00  13.259249 0.583 0.000  0.000000e+00      13    Aplnr
Sox7      0.000000e+00  12.735374 0.396 0.000  0.000000e+00      13     Sox7
Gm20631   0.000000e+00  11.945283 0.250 0.000  0.000000e+00      13  Gm20631
Bcl6b     0.000000e+00  11.791917 0.542 0.000  0.000000e+00      13    Bcl6b
Fam167b   0.000000e+00  11.527533 0.625 0.001  0.000000e+00      13  Fam167b
Cd34      0.000000e+00  10.811053 0.604 0.000  0.000000e+00      13     Cd34
Sema3f    0.000000e+00  10.805728 0.417 0.000  0.000000e+00      13   Sema3f
Kdr       0.000000e+00  10.773707 0.979 0.002  0.000000e+00      13      Kdr
Ptprb     0.000000e+00  10.594309 0.979 0.002  0.000000e+00      13    Ptprb
Stab2     0.000000e+00  10.456647 0.750 0.006  0.000000e+00      13    Stab2
Tek       0.000000e+00  10.188056 0.812 0.002  0.000000e+00      13      Tek
Gpihbp1   0.000000e+00  10.119093 0.750 0.001  0.000000e+00      13  Gpihbp1
Cldn5     0.000000e+00   9.965325 0.688 0.001  0.000000e+00      13    Cldn5
Arhgef15  0.000000e+00   9.912789 0.500 0.001  0.000000e+00      13 Arhgef15
Galnt15   0.000000e+00   9.882335 0.354 0.001  0.000000e+00      13  Galnt15
Gm16268  6.129661e-246   9.797359 0.292 0.001 1.144224e-241      13  Gm16268
Nova2     0.000000e+00   9.792615 0.542 0.001  0.000000e+00      13    Nova2
Prex2     0.000000e+00   9.747728 0.958 0.004  0.000000e+00      13    Prex2




### another dataset


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
Layers(seurat)

head(GetAssayData(seurat, layer = "counts"))




seurat <- NormalizeData(seurat)
head(GetAssayData(seurat, layer = "data"))
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
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk", "Sost", "Dmp1"),
  ncol = 3
)

DotPlot(
  seurat,
  features = c("Runx2","Sp7","Alpl","Bglap","Col1a1","Acp5","Ctsk")
) + RotatedAxis()

table(Idents(seurat))

grep("sost", rownames(seurat), value = TRUE, ignore.case = TRUE) ### seems the sost and dmp1 are not within the datasets
grep("Sost", rownames(seurat), value = TRUE)

grep("^Robo[1-4]$", rownames(seurat), value = TRUE)


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



FeaturePlot(seurat, features=c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"))
DotPlot(seurat, features = c("Slit1","Slit2","Slit3","Robo1","Robo2","Robo3","Robo4"),
        group.by="seurat_clusters") + RotatedAxis()






### with 6 pca based on elbow; I have tried the highest possible pca number but the best seems to be 6
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


### cluster 18 has highest slit1, robo3 and robo4, cluster 15 robo2-3 and cluster 11 slit 1-3
cluster18_markers <- FindMarkers(
  seurat,
  ident.1 = 18,
  only.pos = TRUE
)

head(cluster18_markers, 20)


# or markers <- FindAllMarkers(
markers <- FindAllMarkers(
  seurat,
  only.pos=TRUE,
  min.pct=0.25,
  logfc.threshold=0.25
)

head(markers)

cluster18_markers <- markers %>%
  filter(cluster == 18) %>%
  arrange(desc(avg_log2FC))

head(cluster18_markers, 20)
cluster18_markers$gene <- row.names(cluster18_markers)
FeaturePlot(seurat, features = head(cluster18_markers$gene, 6))


#### results 
head(cluster18_markers, 20)
head(cluster18_markers, 20)
p_val avg_log2FC pct.1 pct.2     p_val_adj
Arpp21         0.000000e+00   8.723548 0.929 0.007  0.000000e+00
Scn4b          0.000000e+00   9.722526 0.571 0.001  0.000000e+00
Gm30948        0.000000e+00   9.692579 0.536 0.001  0.000000e+00
2010300C02Rik  0.000000e+00   8.602413 0.536 0.002  0.000000e+00
Igll1         3.131152e-307   8.564443 1.000 0.015 5.029257e-303
Mmrn1         9.250365e-307   6.463341 0.679 0.005 1.485794e-302
Smtnl2        4.792562e-258   6.794571 0.821 0.011 7.697813e-254
Vpreb2        7.381909e-252   7.958890 0.536 0.004 1.185682e-247
Heyl          8.210151e-250   9.513071 0.393 0.002 1.318714e-245
Fxyd6         8.345889e-241  10.487336 0.250 0.000 1.340517e-236
Fcrl6         1.810423e-228   9.551710 0.286 0.001 2.907902e-224
Ttll11        9.642693e-218   6.825246 0.357 0.002 1.548809e-213
Dntt          3.526425e-217   7.720488 1.000 0.023 5.664144e-213
Vpreb1        1.558388e-202   8.388607 1.000 0.026 2.503082e-198
Gfra2         2.817893e-183   6.030227 0.750 0.015 4.526099e-179
1700048O20Rik 7.587524e-173   7.364736 0.536 0.007 1.218708e-168
A630023P12Rik 3.337016e-170   7.621564 0.357 0.003 5.359915e-166
Uaca          6.580529e-161   6.215528 0.429 0.005 1.056965e-156
Slit1         1.794778e-160   6.441486 0.429 0.005 2.882772e-156
Dlg2          1.010747e-157   7.837119 0.536 0.008 1.623462e-153
