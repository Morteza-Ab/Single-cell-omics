install.packages("Seurat")
library(Seurat)
library(ggplot2)
#remove.packages("ggrepel")
#install.packages("ggrepel")
library(ggrepel)
#install.packages("SingleR")
BiocManager::install("SingleR")

library(SingleR)
library(dplyr)
BiocManager::install("celldex")

library(celldex)
BiocManager::install("SingleCellExperiment")

library(RColorBrewer)
BiocManager::install("RColorBrewer")

library(SingleCellExperiment)
setwd("/Users/mortezaabyadeh/Desktop")
AST <- readRDS("AST.rds")
head(AST)
dim(AST)
colnames(AST)
class(AST)

P10 <- readRDS("13.main_analysis_seuratobject-191210.rds")
class(P10)
dim(P10)
head(P10)


### this might works ###

detach("package:Seurat", unload = TRUE)
library(SeuratObject)
P10 <- as(P10, "Seurat")
class(P10)



saveRDS(P10, file = "P10_temp.rds")
library(SeuratObject)
P11 <- readRDS("P10_temp.rds")
class(P11)

head(P11)

class(P11)
dim(P11)

packageVersion("Seurat")

remove.packages("Seurat", lib = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library")

install.packages("remotes")

remotes::install_version("Seurat", version = "4.4.0")

packageVersion("Seurat")
update.packages("Seurat")
P10 <- readRDS("13.main_analysis_seuratobject-191210.rds")
class(P10)
dim(P10)
head(P10)



head(AST)
summary(AST)
class(AST)
ElbowPlot(AST)
srat <- FindNeighbors(AST, dims = 1:10)
srat <- FindClusters(AST, resolution = 0.5)
srat <- RunUMAP(AST, dims = 1:10, verbose = F)
head(srat)
table(srat@meta.data$seurat_clusters)
DimPlot(AST,label.size = 4,repel = T,label = T, group.by = "nomatch_renamed")
DimPlot(srat,label.size = 4,repel = T, label= T, group.by = "seurat_clusters")
DimPlot(srat,label.size = 4,repel = T,label = T, group.by = "Group")
DimPlot(srat,label.size = 4,repel = T,label = T, group.by = "FINAL")



print(unique(AST$FINAL))
print(unique(AST$nomatch_renamed))
FeaturePlot(srat, features = c("Aldh1l1", "Sox9"))
VlnPlot(srat,features = "Sox9") & theme(plot.title = element_text(size=10))
levels(Idents(AST))
levels(Idents(srat))
srat1 <- srat
head(AST)
levels(Idents(srat1)) <- unique(AST$nomatch_renamed)
levels(Idents(srat1))
#DimPlot(srat,label.size = 4,repel = T,label = T)
