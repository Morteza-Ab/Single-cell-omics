library(Seurat)
library(dplyr)
library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(harmony)
library(ggplot2)
library(hdf5r)
library(scales)
library(DESeq2)
library(SingleR)
library(celldex)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "Seurat",
  "harmony",
  "GEOquery",
  "ggplot2",
  "patchwork",
  "dplyr",
  "hdf5r",
  "glmGamPoi"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(Seurat)
library(harmony)
library(GEOquery)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(glmGamPoi)



getGEOSuppFiles("GSE171768", makeDirectory = TRUE)

untar(
  "GSE171768/GSE171768_RAW.tar",
  exdir = "GSE171768/raw"
)

list.files("GSE171768/raw", recursive = TRUE)


s1 <- Read10X_h5(
  "GSE171768/raw/GSM5233644_20pcO2_r1_filtered_feature_bc_matrix.h5"
)

s2 <- Read10X_h5(
  "GSE171768/raw/GSM5233645_20pcO2_r2_filtered_feature_bc_matrix.h5"
)

s3 <- Read10X_h5(
  "GSE171768/raw/GSM5233646_5pcO2_r1_filtered_feature_bc_matrix.h5"
)

s4 <- Read10X_h5(
  "GSE171768/raw/GSM5233647_5pcO2_r2_filtered_feature_bc_matrix.h5"
)


s1 <- CreateSeuratObject(counts = s1, project = "O20-r1", min.cells = 10, min.features = 200)
s2 <- CreateSeuratObject(counts = s2, project = "O20-r2", min.cells = 10, min.features = 200)
s3 <- CreateSeuratObject(counts = s3, project = "O5-r1",  min.cells = 10, min.features = 200)
s4 <- CreateSeuratObject(counts = s4, project = "O5-r2",  min.cells = 10, min.features = 200)


s1$Oxygen <- "20"; s1$sample <- "O20r1"
s2$Oxygen <- "20"; s2$sample <- "O20r2"
s3$Oxygen <- "5";  s3$sample <- "O5r1"
s4$Oxygen <- "5";  s4$sample <- "O5r2"


################# OR



getGEOSuppFiles("GSE171768", makeDirectory = TRUE)

untar(
  "GSE171768/GSE171768_RAW.tar",
  exdir = "GSE171768/raw"
)

files <- list.files("GSE171768/raw", pattern = "h5$", full.names = TRUE)

names(files) <- c("O20r1","O20r2","O5r1","O5r2")


## CREATE OBJECTS

objs <- lapply(names(files), function(nm) {
  
  counts <- Read10X_h5(files[nm])
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = nm,
    min.cells = 10,
    min.features = 200
  )
  
  obj$sample <- nm
  obj$Oxygen <- ifelse(grepl("O20", nm), "20", "5")
  
  obj
})

names(objs) <- names(files)

####QC METRICS

objs <- lapply(objs, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})


## FILTER

objs <- lapply(objs, function(obj) {
  subset(obj,
         subset =
           nFeature_RNA > 200 &
           nFeature_RNA < 7500 &
           percent.mt < 15)
})

### NORMALIZATION (PER SAMPLE)

objs <- lapply(objs, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj
})



anchors <- FindIntegrationAnchors(
  object.list = objs,
  dims = 1:20
)

combined <- IntegrateData(anchors)




DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30)

combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:20)

Idents(combined) <- "seurat_clusters"

DimPlot(combined, label = TRUE)
DimPlot(combined, group.by = "sample")
DimPlot(combined, group.by = "Oxygen")


combined$condition_cluster <- paste(combined$Oxygen, Idents(combined), sep = "_")

Idents(combined) <- "condition_cluster"




###########Example: cluster-wise oxygen comparison
deg_list <- list()

clusters <- levels(Idents(combined))

for (cl in clusters) {
  
  deg_list[[cl]] <- FindMarkers(
    combined,
    ident.1 = paste0("5_", str_split(cl, "_")[[1]][2]),
    ident.2 = paste0("20_", str_split(cl, "_")[[1]][2]),
    logfc.threshold = 0
  )
  
  deg_list[[cl]]$cluster <- cl
}



make_volcano <- function(df, title) {
  df$log2fc <- log2(exp(df$avg_logFC))
  df$gene <- rownames(df)
  
  ggplot(df, aes(log2fc, -log10(p_val))) +
    geom_point(alpha = 0.5) +
    ggtitle(title) +
    theme_classic()
}

volcano_plots <- mapply(make_volcano, deg_list, names(deg_list), SIMPLIFY = FALSE)



#saveRDS(combined, "GSE171768_integrated.rds")
#write.csv(markers, "cluster_markers.csv")
