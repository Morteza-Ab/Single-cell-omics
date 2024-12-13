
setwd("Users/mortezaabyadeh/Desktop/BG")
BG<-readRDS("BG.rds")

install.packages("SeuratObject")
install.packages("Seurat")
library(SeuratObject)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)



# Specify the clusters to keep/BG_prog
clusters_to_keep <- c("P10_BG", "P10_PROG","P10_BG_KO",  "P10_PROG_KO",
                      "P17_BG", "P17_BG_KO",
                      "P17_PROG_KO")
levels(Idents(AST))
# Set identities to the column you are interested in
Idents(AST) <- AST$FINAL

BG_PROG <- subset(AST, idents = clusters_to_keep)

umap_plot<- DimPlot(BG_PROG, 
                    reduction = "umap",group.by = "FINAL", 
                    label = F, 
                    label.size = 3,  # Adjust label size here
                    pt.size = 0.5,   # Adjust point size
                    shuffle = F) # Shuffle the labels to avoid overlap

# Add labels with repelling
LabelClusters(umap_plot, id = "FINAL", repel = TRUE, size = 3)
BG_PROG <- saveRDS(BG_PROG, file = "BG_PROG.rds")
head(BG_PROG)
####
# Specify the clusters to keep/BG
clusters_to_keep <- c("P10_BG", "P10_BG_KO",  
                      "P17_BG", "P17_BG_KO")
levels(Idents(AST))
# Set identities to the column you are interested in
Idents(AST) <- AST$FINAL

BG <- subset(AST, idents = clusters_to_keep)

umap_plot<- DimPlot(BG, 
                    reduction = "umap",group.by = "FINAL", 
                    label = F, 
                    label.size = 3,  # Adjust label size here
                    pt.size = 0.5,   # Adjust point size
                    shuffle = F) # Shuffle the labels to avoid overlap
BG <- saveRDS(BG, file = "BG.rds")
# Add labels with repelling
LabelClusters(umap_plot, id = "FINAL", repel = TRUE, size = 3)
setwd("C:/Users/zareikheirm/Desktop/relabeled/BG")
dir()
BG<- readRDS("C:/Users/zareikheirm/Desktop/relabeled/BG/BG.rds")

# View the cluster identities
levels(Idents(BG))


head(BG)

unique(BGT@meta.data$Age)




# Plot UMAP
DimPlot(BG, reduction = "umap", label = TRUE,  pt.size = 0.3, 
        group.by= "FINAL" )
umap_plot<- DimPlot(BG, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.3, 
        group.by = "nomatch_renamed", repel = TRUE)
# Remove the title by setting it to NULL
umap_plot <- umap_plot + ggtitle(NULL)
# Display the plot
print(umap_plot)

saveRDS(BG, file = "BG.rds")
BG = readRDS("C:/Users/zareikheirm/Documents/p10-p17/relabeled/BG/BG.rds")
head(BG)

#### BG Top genes
top_genes_all_clusters <- FindAllMarkers(BG)
head(top_genes_all_clusters)

library(EnhancedVolcano)
library(ggrepel)


# View the cluster identities
levels(Idents(BG))

# Set Idents to RNA_snn_res.0.5
BG <- SetIdent(BG, value = "nomatch_renamed")

# Check the identity levels
levels(Idents(BG))

# Perform differential expression analysis
de_results <- FindMarkers(BG, ident.1 = "P10_BG" , ident.2 = "P17_BG")

# Add columns for volcano plot
de_results$log2FoldChange <- de_results$avg_log2FC
de_results$pvalue <- de_results$p_val_adj

# Create the volcano plot
EnhancedVolcano(de_results,
                lab = rownames(de_results),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = "Mean Log Fold Change",
                ylab = "-log10(p-value)",
                title = "wt BG: P17 vs P10",
                pCutoff = 0.05,  # Adjust the p-value cutoff
                FCcutoff = 1.5,  # Adjust the log2 fold change cutoff
                pointSize = 1.0,  # Size of the points
                labSize = 3,    # Size of the gene labels
                col = c('black', 'red'), # Colors for non-significant and significant points
                legendLabels = c('NS', 'Significant'),
                legendPosition = 'top',
                legendLabSize = 10)

##################P10
# Perform differential expression analysis
de_results_P10 <- FindMarkers(BG, ident.1 = "P10_BG", ident.2 = "P10_BG_KO")

# Add columns for volcano plot
de_results_P10$log2FoldChange <- de_results_P10$avg_log2FC
de_results_P10$pvalue <- de_results_P10$p_val_adj


p<-EnhancedVolcano(de_results_P10,
                lab = rownames(de_results_P10),
                x = "avg_log2FC",  # Column with log fold change
                y = "p_val_adj",  # Column with adjusted p-values
                xlab = "Log2 Fold Change",
                ylab = "-log10(p-value)",
                title = "BG-P10: WT vs KO",
                
                pCutoff = 0.05,  # Significance threshold for p-value
                FCcutoff = 1.0,  # Significance threshold for fold change
                pointSize = 1.5,
                labSize = 3.0,
                col = c('gray', 'red', 'blue', 'purple'),  # Color for non-sig, upregulated, downregulated, and both sig
                colAlpha = 0.6,
                legendPosition = "top",
                legendLabSize = 10,
                legendIconSize = 5,
                gridlines.major = TRUE,
                gridlines.minor = FALSE, max.overlaps = 60 )

# Modify the font size for axis labels and title
p <- p + theme(axis.title.x = element_text(size = 10),  # X-axis label font size
               axis.title.y = element_text(size = 10),  # Y-axis label font size
               plot.title = element_text(size = 12),    # Plot title font size
               plot.subtitle = element_text(size = 5), # Subtitle font size
               axis.text.x = element_text(size = 10),    # Font size for x-axis tick labels
               axis.text.y = element_text(size = 10))     # Axis tick labels font size

# Display the plot
print(p)
head(de_results_P10, n = 20)

# Ppia, Them4, Ddx5, Vapa, Csnk1a1, Bex3, Cirbp, H3f3a, Gdf10, Tsix, Ndufb11, 
# Bex2, Tmsb4x, Hnrnpa2b1, Srsf2, Cox8a, Rpl41, Rps13 
 

# Sort by p-value
top_genes_pvalue <- de_results %>%
  arrange(p_val_adj) %>%          # Arrange by adjusted p-value
  head(n = 10)                     # Select top 10 genes

# Sort by fold change
top_genes_fc <- de_results %>%
  arrange(desc(abs(avg_log2FC))) %>% # Arrange by absolute fold change
  head(n = 10)     

top_genes_combined <- de_results %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>% # Arrange by p-value and then fold change
  head(n = 10)  # Select top 10 genes
# Set cutoffs
pvalue_cutoff <- 0.05
fc_cutoff <- 1.0

significant_genes <- de_results %>%
  filter(p_val_adj < pvalue_cutoff & abs(avg_log2FC) > fc_cutoff)

# Create volcano plot with highlighted top genes
P<- EnhancedVolcano(de_results,
                    lab = rownames(de_results),
                    x = "avg_log2FC",
                    y = "p_val_adj",
                    xlab = "Log2 Fold Change",
                    ylab = "-log10(p-value)",
                    title = "BG-p10: KO vs WT",
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    col = c('gray', 'red', 'blue', 'purple'),
                    colAlpha = 0.6,
                    legendPosition = "top",
                    legendLabSize = 10,
                    
                    gridlines.major = TRUE,
                    gridlines.minor = FALSE, labSize = 3.0)

# Modify the font size for axis labels and title
p <- p + theme(axis.title.x = element_text(size = 10),  # X-axis label font size
               axis.title.y = element_text(size = 10),  # Y-axis label font size
               plot.title = element_text(size = 12),    # Plot title font size
               plot.subtitle = element_text(size = 5), # Subtitle font size
               axis.text.x = element_text(size = 10),    # Font size for x-axis tick labels
               axis.text.y = element_text(size = 10))     # Axis tick labels font size

# Display the plot
print(p)



########################## P17-BG

# Perform differential expression analysis
de_results <- FindMarkers(BG, ident.1 = "P17_BG", ident.2 = "P17_BG_KO")

# Add columns for volcano plot
de_results$log2FoldChange <- de_results$avg_log2FC
de_results$pvalue <- de_results$p_val_adj


p<-EnhancedVolcano(de_results,
                   lab = rownames(de_results),
                   x = "avg_log2FC",  # Column with log fold change
                   y = "p_val_adj",  # Column with adjusted p-values
                   xlab = "Log2 Fold Change",
                   ylab = "-log10(p-value)",
                   title = "BG-p17: WT vs KO",
                   pCutoff = 0.05,  # Significance threshold for p-value
                   FCcutoff = 1.0,  # Significance threshold for fold change
                   pointSize = 1.5,
                   labSize = 3.0,
                   col = c('gray', 'red', 'blue', 'purple'),  # Color for non-sig, upregulated, downregulated, and both sig
                   colAlpha = 0.6,
                   legendPosition = "top",
                   legendLabSize = 10,
                   legendIconSize = 5,
                   drawConnectors = F,
                   widthConnectors = 0.1,
                   gridlines.major = TRUE,
                   gridlines.minor = FALSE, max.overlaps = 20)

# Modify the font size for axis labels and title
p <- p + theme(axis.title.x = element_text(size = 10),  # X-axis label font size
               axis.title.y = element_text(size = 10),  # Y-axis label font size
               plot.title = element_text(size = 12),    # Plot title font size
               plot.subtitle = element_text(size = 5), # Subtitle font size
               axis.text.x = element_text(size = 10),    # Font size for x-axis tick labels
               axis.text.y = element_text(size = 10))     # Axis tick labels font size

# Display the plot
print(p)

print(de_results, 10)

# Assuming de_results is your data frame with differential expression results
library(dplyr)

# Sort by p-value
top_genes_pvalue <- de_results %>%
  arrange(p_val_adj) %>%          # Arrange by adjusted p-value
  head(n = 10)                     # Select top 10 genes

# Sort by fold change
top_genes_fc <- de_results %>%
  arrange(desc(abs(avg_log2FC))) %>% # Arrange by absolute fold change
  head(n = 10)     

top_genes_combined <- de_results %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>% # Arrange by p-value and then fold change
  head(n = 10)  # Select top 10 genes
# Set cutoffs
pvalue_cutoff <- 0.05
fc_cutoff <- 1.0

significant_genes <- de_results %>%
  filter(p_val_adj < pvalue_cutoff & abs(avg_log2FC) > fc_cutoff)

# Create volcano plot with highlighted top genes
P<- EnhancedVolcano(de_results,
                lab = rownames(de_results),
                x = "avg_log2FC",
                y = "p_val_adj",
                xlab = "Log2 Fold Change",
                ylab = "-log10(p-value)",
                title = "BG-p17: KO vs WT",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                col = c('gray', 'red', 'blue', 'purple'),
                colAlpha = 0.6,
                legendPosition = "top",
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = F,
                widthConnectors = 0.5,
                gridlines.major = TRUE,
                gridlines.minor = FALSE, labSize = 3.0)

# Modify the font size for axis labels and title
p <- p + theme(axis.title.x = element_text(size = 10),  # X-axis label font size
               axis.title.y = element_text(size = 10),  # Y-axis label font size
               plot.title = element_text(size = 12),    # Plot title font size
               plot.subtitle = element_text(size = 5), # Subtitle font size
               axis.text.x = element_text(size = 10),    # Font size for x-axis tick labels
               axis.text.y = element_text(size = 10))     # Axis tick labels font size

# Display the plot
print(p)


head(de_results, n = 20)   # Display the first 20 rows

# Cldn10, Paqr6, Ntsr2, Scrg1, Nkain4, Dao, Cd9, Tuba1a, Gdf10,
# Aldoc, Pla2g16, Cd81, Ramp1, Csrp1, Slc7a10, Ech1, Ednrb, Psmb5, Cst3, Gria4 


### dotplot with new genes
markers_of_interest <- c( "Ephx1", "Ccnd2", "Tagin2", "Tnfrsf12a", "Rps27l", "Asns",
                          "Iqgap2", "Dohh",
                          "Them4", "Bex3", "Ech1", "Gpx8","Gdf10", "Gpr37l1",
                          "Ndufb11", "S100b", "S1pr1", "Slc1a2", "Slc1a3", "Aqp4", 
                          "Ldhb", "Apoe", "Glul", "Gja1", "Gjb6", "Gpc5")
)
# Custom color palette with at least 8 colors
colors <- c("red", "blue", "green", "purple", "orange", "pink", "yellow", "cyan")

DotPlot(BG, features = markers_of_interest, cols = colors, group.by = "nomatch_renamed", scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.title.x = element_text(size = 12),  # X-axis label font size
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10))  # Rotate x-axis labels


