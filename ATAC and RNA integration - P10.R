
#--------------------------------------------------------
## Alternative methods to consider:
# Finding differentially accessible peaks between clusters: https://docs.scvi-tools.org/en/0.15.1/tutorials/notebooks/peakvi_in_R.html
# Try integration with scRNAseq using scANVI: https://docs.scvi-tools.org/en/0.15.1/tutorials/notebooks/peakvi_in_R.html
#--------------------------------------------------------

# This analysis is based on the pipeline here: https://stuartlab.org/signac/articles/mouse_brain_vignette.html#integrating-with-scrna-seq-data

#--------------------------------------------------------
# Load libraries
#--------------------------------------------------------
library(Signac)

library(Seurat)
library(ggplot2)
library(SeuratObject)
#library(GenomeInfoDb)
#library(EnsDb.Mmusculus.v79) # EnsDb.Mmusculus.v79 is mm10: from https://stuartlab.org/signac/articles/install.html
#library(patchwork)

rm(list=ls())

#--------------------------------------------------------
# load filtered and processed data
#--------------------------------------------------------
setwd("/Users/mortezaabyadeh/Desktop")
dir()
brain_atac_P10 = readRDS("brain_after_4593cells.rds")
head(brain_atac_P10)
brain_atac_P10@meta.data$Condition = ifelse(brain_atac_P10@meta.data$group %in% c("WT1", "WT2"), "WT", "KO")
brain_atac_P10@meta.data$cell.id = row.names(brain_atac_P10@meta.data)


setwd("E:/SCRNA-seq/p10-p17/relabeled/11042024 analysis")
rna_P10P17 = readRDS("Double_positive_WTKO_P10P17_11_04_2024.rds")

#--------------------------------------------------------
# Do some plotting to orient yourself
#--------------------------------------------------------

DimPlot(brain_atac_P10, group.by = 'Condition', label = TRUE, repel = TRUE)  + ggtitle('P10 scATAC-seq')
DimPlot(brain_atac_P17, group.by = 'Condition', label = TRUE, repel = TRUE)  + ggtitle('P17 scATAC-seq')

DimPlot(rna_P10P17, group.by = 'Category', label = TRUE, repel = TRUE)  + ggtitle('scRNA-seq')
DimPlot(rna_P10P17, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)  + ggtitle('scRNA-seq')
head(rna_P10P17)
#--------------------------------------------------------
# Transfer RNA labels on P10 ATAC cells
#--------------------------------------------------------
set.seed(123)

DefaultAssay(brain_atac_P10) <- "RNA"



rna_P10 <- subset(rna_P10P17, subset = Age == "P10")

rna_P10 <- FindVariableFeatures(
  object = rna_P10,
  nfeatures = 2000
)

transfer.anchors <- FindTransferAnchors(
  reference = rna_P10,
  query = brain_atac_P10,
  reduction = 'cca',
  dims = 1:40
)
# In Seurat, Canonical Correlation Analysis (CCA) is just one method used for aligning and integrating multi-modal datasets (like RNA-seq and ATAC-seq).
# Warning: npcs is smaller than the largest value requested by the dims parameter.
# Setting npcs to 40 and continuing.
# Running CCA

transfer.anchors
# AnchorSet object containing 3022 anchors between the reference and query

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = as.character(rna_P10$Category),
  weight.reduction = brain_atac_P10[['lsi']],
  dims = 2:30
)

brain_atac_P10 <- AddMetaData(object = brain_atac_P10, metadata = predicted.labels)

plot1 <- DimPlot(rna_P10, group.by = 'Category', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(brain_atac_P10, group.by = 'predicted.id', label = TRUE, repel = TRUE)  + ggtitle('scATAC-seq')

#png("atac_and_rna.png", units="in", width=10, height=6, res=600)
plot1 + plot2
#dev.off()
head(brain_atac_P10)
saveRDS(brain_atac_P10, file = "brain_scrna_atac_P10.rds")
head(Brain_p10)

colnames(brain_atac_P10@meta.data)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = as.character(rna_P10$seurat_clusters),
  weight.reduction = brain_atac_P10[['lsi']],
  dims = 2:30
)

brain_atac_P10 <- AddMetaData(object = brain_atac_P10, metadata = predicted.labels)

head(brain_atac_P10)
# Define a color palette with named colors for each cluster
cluster_colors <- c("2" = "#FF6347", "3" = "#4682B4", "5" = "#32CD32", "7" = "#9370DB")

# Plot scRNA-seq with consistent colors
plot1 <- DimPlot(rna_P10, group.by = 'seurat_clusters', label = TRUE, repel = TRUE, pt.size = 0.5) + 
  NoLegend() + 
  ggtitle('scRNA-seq') +
  scale_color_manual(values = cluster_colors)
  

# Plot scATAC-seq with the same color scheme
plot2 <- DimPlot(brain_atac_P10, group.by = 'predicted.id', label = TRUE, repel = TRUE, pt.size = 0.5) + 
  ggtitle('scATAC-seq with Transferred scRNA Clusters') +
  scale_color_manual(values = cluster_colors)

# Combine the plots
plot1 + plot2




# Define the genes or features you want to visualize in the DotPlot
features_to_plot <- c("Gdf10", "Gpr37l1", "Gria1", "Wif1", "Sept4", "Tnc", "Hsd11b1", "Cyp26b1", 
                      "Myoc", "Nfkb1", "Cxcl10", "Aqp4", "C4b", "Olig2", "Slc6a11",
                      "Slc7a10", "Fam212b", "Kcnip3", "Fgfr3", "Ppdpf", "Hist1h2bc",
                      "Rab11b", "Dtymk", "Mki67", "Top2a", "Mcm2", "Cdk1", "Gfap", "Aldh1l1", 
                      "Sox9", "Vim", "S100b", "S1pr1", "Slc1a2", "Slc1a3", "Apoe"))  # Replace with your genes of interest

# Define your groups of features
group1 <- c("Gdf10", "Gpr37l1", "Gria1", "Wif1", "Sept4", "Tnc")
group2 <- c("Hsd11b1", "Cyp26b1")
group3 <- c("Myoc", "Nfkb1", "Cxcl10", "Aqp4", "C4b")
group4 <- c("Olig2", "Slc6a11",
            "Slc7a10", "Fam212b", "Kcnip3", "Fgfr3", "Ppdpf")
group5 <- c("Hist1h2bc",
            "Rab11b", "Dtymk", "Mki67", "Top2a", "Mcm2", "Cdk1")
group6 <- c("Gfap", "Aldh1l1", 
            "Sox9", "Vim", "S100b", "S1pr1", "Slc1a2", "Slc1a3", "Apoe")



# Combine into one vector with a separator
features_ordered <- c(group1, group2, group3, 
                      group4,  group5,group6)

# Calculate positions for vertical lines based on group lengths
x_intercepts <- c(length(group1) + 0.5, 
                  length(group1) + length(group2) + 0.5, 
                  length(group1) + length(group2) + length(group3) + 0.5,
                  length(group1) + length(group2) + length(group3) + length(group4) + 0.5,
                  length(group1) + length(group2) + length(group3) + length(group4) +length(group5)+ 0.5 )
# Custom color palette with at least 8 colors
colors <- c("red", "blue", "green", "purple", "orange", "pink", "yellow", "cyan")



# Generate the DotPlot for the transferred labels in the scATAC-seq data
DotPlot(brain_atac_P10, features = features_ordered, cols = colors,  group.by = "predicted.id") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), 
        axis.text.y = element_text(size = 8))  +
  geom_vline(xintercept = x_intercepts, linetype = "dashed", color = "black", size = 1)+
  ggtitle("DotPlot for Transferred scRNA-seq Labels in scATAC-seq Data")


# Generate the DotPlot for the transferred labels in the scATAC-seq data
DotPlot(brain_atac_P10, features = features_ordered, cols = colors,  group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), 
        axis.text.y = element_text(size = 8))  +
  geom_vline(xintercept = x_intercepts, linetype = "dashed", color = "black", size = 1)+
  ggtitle("DotPlot for Transferred scRNA-seq Labels in scATAC-seq Data")

# Save the metadata from the Seurat object to a .csv file
write.csv(brain_atac_P10@meta.data, file = "brain_atac_P10_metadata.csv", row.names = TRUE)



