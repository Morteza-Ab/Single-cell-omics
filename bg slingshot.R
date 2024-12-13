#---------------
setwd("/Users/mortezaabyadeh/Desktop")
dir()
AST = readRDS("AST.rds")
levels(Idents(AST))
head(AST)

table(AST$Category)


# Specify the clusters to keep/BG_prog wt
clusters_to_keep <- c("BG_P10_WT", "BGPROG_P10_WT", 
                      "BG_P17_WT", "BGPROG_P17_WT")

# Subset the Seurat object
wtBG_Prog <- subset(BG_PROG, idents = clusters_to_keep)

umap_plot<- DimPlot(wtBG_Prog, 
                    reduction = "umap",group.by = "Category", 
                    label = F, 
                    label.size = 3,  # Adjust label size here
                    pt.size = 0.5,   # Adjust point size
                    shuffle = F) # Shuffle the labels to avoid overlap


# Add labels with repelling
LabelClusters(umap_plot, id = "Category", repel = TRUE, size = 3)

# Assuming 'BG' is your Seurat object
sce <- as.SingleCellExperiment(wtBG_Prog)
# Apply slingshot on the SingleCellExperiment object
sce <- slingshot(sce, clusterLabels = 'Category', reducedDim = 'PCA')



unique(sce$Category)
colors <- c("BGPROG_P10_WT" = "blue","BG_P10_WT" = "red",  "BGPROG_P17_WT" = "green",
            "BG_P17_WT"="yellow")

head(reducedDims(sce)$UMAP)




# Create the data frame with UMAP coordinates and clusters
df <- as.data.frame(reducedDims(sce)$UMAP)
df$cluster <- sce$Category

# Ensure column names are correct
colnames(df) <- c("umap_1", "umap_2", "cluster")


library(slingshot)
sce <- as.SingleCellExperiment(wtBG_Prog)
sce <- slingshot(sce, clusterLabels = 'Category', reducedDim = 'UMAP')


slingshot_results <- slingCurves(sce)

# Extract UMAP coordinates and trajectory data
umap_coords <- as.data.frame(reducedDims(sce)$UMAP)
umap_coords$cluster <- sce$Category

# Extract trajectories
trajectories <- slingshot_results$lineages

# Create ggplot
p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size=0.5) +
  scale_color_manual(values = colors) +  # Apply custom colors
  theme_minimal() +
  labs(title = "UMAP Projection with Slingshot Trajectory",
       x = "UMAP 1",
       y = "UMAP 2")

# Add trajectory lines
for (trajectory in trajectories) {
  p <- p + geom_path(data = trajectory, aes(x = umap_1, y = umap_2), color = "black", size = 1)
}

print(p)
str(slingshot_results)

# Plot the trajectory with UMAP
plot(reducedDims(sce)$UMAP, col = sce$seurat_clusters, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')

setwd("E:/SCRNA-seq/p10-p17/relabeled/11052024/BG")
BG_PROG <- readRDS("BG_PROG.rds")
head(BG_PROG)

table(BG_PROG$Category)

# ############################Specify the clusters to keep/BG_prog ko
clusters_to_keep <- c("BG_P10_cKO", "BGPROG_P10_cKO", "BG_P17_cKO", "BGPROG_P17_cKO")

# Subset the Seurat object
koBG_Prog <- subset(BG_PROG, idents = clusters_to_keep)
table(koBG_Prog$Category)
# BG_P10_cKO     BG_P17_cKO BGPROG_P10_cKO BGPROG_P17_cKO 
#  118            498             20             32 
umap_plot<- DimPlot(koBG_Prog, 
                    reduction = "umap",group.by = "Category", 
                    label = F, 
                    label.size = 3,  # Adjust label size here
                    pt.size = 0.5,   # Adjust point size
                    shuffle = F) # Shuffle the labels to avoid overlap


# Add labels with repelling
LabelClusters(umap_plot, id = "Category", repel = TRUE, size = 3)

# Assuming 'BG' is your Seurat object
sce <- as.SingleCellExperiment(koBG_Prog)
# Apply slingshot on the SingleCellExperiment object
sce <- slingshot(sce, clusterLabels = 'Category', reducedDim = 'PCA')



unique(sce$Category)
colors <- c("BGPROG_P10_cKO" = "blue","BG_P10_cKO" = "red", "BGPROG_P17_cKO" = "green", "BG_P17_cKO"="pink")

head(reducedDims(sce)$UMAP)




# Create the data frame with UMAP coordinates and clusters
df <- as.data.frame(reducedDims(sce)$UMAP)
df$cluster <- sce$Category

# Ensure column names are correct
colnames(df) <- c("umap_1", "umap_2", "cluster")


library(slingshot)
sce <- as.SingleCellExperiment(koBG_Prog)
sce <- slingshot(sce, clusterLabels = 'Category', reducedDim = 'UMAP')


slingshot_results <- slingCurves(sce)

# Extract UMAP coordinates and trajectory data
umap_coords <- as.data.frame(reducedDims(sce)$UMAP)
umap_coords$cluster <- sce$Category

# Extract trajectories
trajectories <- slingshot_results$lineages

# Create ggplot
p <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size=0.5) +
  scale_color_manual(values = colors) +  # Apply custom colors
  theme_minimal() +
  labs(title = "UMAP Projection with Slingshot Trajectory",
       x = "UMAP 1",
       y = "UMAP 2")

# Add trajectory lines
for (trajectory in trajectories) {
  p <- p + geom_path(data = trajectory, aes(x = umap_1, y = umap_2), color = "black", size = 0.5)
}

print(p)
str(slingshot_results)

# Plot the trajectory with UMAP
plot(reducedDims(sce)$UMAP, col = sce$seurat_clusters, pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')


