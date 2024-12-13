if (!requireNamespace("remotes", quietly = TRUE))
  install.packages('remotes')
remotes::install_github('immunogenomics/presto')
library(Rcpp)
library(data.table)
library(presto)
library(Seurat)
library(SeuratObject)
library(Signac)

### run your atatc seq data
setwd("/Users/mortezaabyadeh/Desktop")
# dir()
# brain <- readRDS('brain_after_4593cells.rds')
# head(brain)
# table(brain$group)
# # KO2  KO1  WT2  WT1 
# # 375 1164  979 2075
brain_P10 <- readRDS("brain_scrna_atac_P10.rds")
head(brain_P10)
table(brain_P10$group)

# KO2  KO1  WT2  WT1 
# 375 1164  979 2075

brain_P10$group_binary <- ifelse(brain_P10@meta.data$group%in%c("KO2","KO1"),"KO","WT")
head(brain_P10)
table(brain_P10$group_binary)
# KO   WT 
# 1539 3054
table(brain_P10$predicted.id)
# BG_P10_cKO      BG_P10_WT BGPROG_P10_cKO  BGPROG_P10_WT   MYOC_P10_cKO 
# 1290           1236            156            268             57 
# MYOC_P10_WT  OLIG2_P10_cKO   OLIG2_P10_WT   VAFA_P10_cKO    VAFA_P10_WT 
# 161             69            581            320            455 

# change back to working with peaks instead of gene activities
DefaultAssay(brain_P10) <- 'peaks'
DefaultAssay(brain_P10)
#[1] "peaks"

# levels(Idents(brain_P10))
Idents(brain_P10)
Idents(brain_P10) <- "predicted.id"

# wilcox is the default option for test.use
da_peaks <- FindMarkers(
  object = brain_P10,
  ident.1 = "BG_P10_cKO",
  ident.2 = "BG_P10_WT",
  test.use = 'wilcox',
  min.pct = 0.1
)
# min.pct stands for the minimum percentage of cells 
# that should express a gene in at least one of the two groups being compared (e.g., "BG_P10_cKO" vs. "BG_P10_WT").

head(da_peaks)
dim(da_peaks)
# [1] 85031     5
da_peaks <- da_peaks[da_peaks$p_val_adj<0.05,]
head(da_peaks)
dim(da_peaks)
# [1] 39234 5
summary(da_peaks)
# p_val             avg_log2FC          pct.1            pct.2          p_val_adj        
# Min.   :0.000e+00   Min.   :-3.1149   Min.   :0.0230   Min.   :0.0280   Min.   :0.0000000  
# 1st Qu.:0.000e+00   1st Qu.: 0.5614   1st Qu.:0.1810   1st Qu.:0.0950   1st Qu.:0.0000000  
# Median :2.490e-12   Median : 0.8047   Median :0.2710   Median :0.1510   Median :0.0000005  
# Mean   :1.440e-08   Mean   : 0.7112   Mean   :0.3145   Mean   :0.1889   Mean   :0.0028970  
# 3rd Qu.:2.423e-09   3rd Qu.: 0.9969   3rd Qu.:0.4070   3rd Qu.:0.2390   3rd Qu.:0.0004874  
# Max.   :2.485e-07   Max.   : 1.9906   Max.   :0.9850   Max.   :0.9470   Max.   :0.0499883  

plot1 <- VlnPlot(
  object = brain_P10,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("BG_P10_cKO", "BG_P10_WT")
)
plot2 <- FeaturePlot(
  object = brain_P10,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

nrow(brain_P10)


# Peak coordinates can be difficult to interpret alone. We can find the closest gene to each of these peaks using the ClosestFeature() function.
######################## up
open_up <- rownames(da_peaks[da_peaks$avg_log2FC > 1, ])


# summary(da_peaks$avg_log2FC)
summary(open_up)
# Length     Class      Mode 
# 9669 character character
head(open_up)
head(da_peaks)
dim(open_up)
length(open_up)
# [1] 9669

closest_genes_up <- ClosestFeature(brain_P10, regions = open_up)


df <- as.data.frame(da_peaks)
rownames(df)
head(da_peaks)
da_peaks$query_region <- row.names(df)

head(da_peaks)

dim(closest_genes_up)
# [1] 9669    8
dim(da_peaks)
#39234     6
################
library(dplyr)

# Perform a left join based on the "query_region" column
overal_up <- left_join(closest_genes_up, da_peaks, by = "query_region")

# View the combined dataset
head(overal_up,2)
# tx_id    gene_name            gene_id   gene_biotype type           closest_region
# 1 ENSMUST00000110844        Kcnh1 ENSMUSG00000058248 protein_coding  gap chr1-192337830-192413016
# 2 ENSMUST00000195148 RP24-411K2.1 ENSMUSG00000103451        lincRNA exon chr1-183563512-183564045
# query_region distance        p_val avg_log2FC pct.1 pct.2    p_val_adj
# 1 chr1-192380397-192381287        0 1.043413e-85   1.139048 0.799 0.392 2.099117e-80
# 2 chr1-183504047-183504970    58541 3.622211e-82   1.069741 0.809 0.393 7.287091e-77
head(closest_genes_up,2)
 # tx_id    gene_name            gene_id   gene_biotype type
 # ENSMUST00000110844 ENSMUST00000110844        Kcnh1 ENSMUSG00000058248 protein_coding  gap
 # ENSMUSE00001345758 ENSMUST00000195148 RP24-411K2.1 ENSMUSG00000103451        lincRNA exon
 # closest_region             query_region distance
 # ENSMUST00000110844 chr1-192337830-192413016 chr1-192380397-192381287        0
 # ENSMUSE00001345758 chr1-183563512-183564045 chr1-183504047-183504970    58541
 
dim(overal_up)
dim(overal_up)
# [1] 9669   13
dim(closest_genes_up)
#[1] 9669    8
#####################

# install.packages("writexl")
library(writexl)
# Assuming 'your_data' is the data frame you want to save
write_xlsx(overal_up, "closest_genes_up.xlsx")

####################### down
open_down <- rownames(da_peaks[da_peaks$avg_log2FC < -1, ])
closest_genes_down <- ClosestFeature(brain_P10, regions = open_down)
head(closest_genes_down)
# tx_id     gene_name            gene_id   gene_biotype type
# ENSMUST00000168338 ENSMUST00000168338       Trmt61a ENSMUSG00000060950 protein_coding  utr
# ENSMUSE00000315744 ENSMUST00000022124         Cd180 ENSMUSG00000021624 protein_coding exon
# ENSMUST00000165114 ENSMUST00000165114       Zfp36l1 ENSMUSG00000021127 protein_coding  cds
# ENSMUST00000032462 ENSMUST00000032462         Timp4 ENSMUSG00000030317 protein_coding  utr
# ENSMUST00000005953 ENSMUST00000005953         Sqrdl ENSMUSG00000005803 protein_coding  utr
# ENSMUSE00000809893 ENSMUST00000153921 4933427I22Rik ENSMUSG00000085928        lincRNA exon
# closest_region              query_region distance
# ENSMUST00000168338 chr12-111682870-111683902 chr12-111693685-111694603     9782
# ENSMUSE00000315744 chr13-102693558-102693789 chr13-102573167-102574085   119472
# ENSMUST00000165114   chr12-80112819-80112875   chr12-80112747-80113635        0
# ENSMUST00000032462  chr6-115251779-115251849  chr6-115251831-115252803        0
# ENSMUST00000005953  chr2-122809339-122809569  chr2-123386920-123387843   577350
# ENSMUSE00000809893  chr4-139943408-139943519  chr4-139948262-139949151     4742
# 
df <- as.data.frame(da_peaks)
rownames(df)
head(da_peaks)
da_peaks$query_region <- row.names(df)

head(da_peaks)

dim(closest_genes_down)
# [1] 1032    8
dim(da_peaks)
# [1] 39234     6
################
library(dplyr)

# Perform a left join based on the "query_region" column
overal_down <- left_join(closest_genes_down, da_peaks, by = "query_region")

# View the combined dataset
head(overal_down,2)
# tx_id gene_name            gene_id   gene_biotype type            closest_region
# 1 ENSMUST00000168338   Trmt61a ENSMUSG00000060950 protein_coding  utr chr12-111682870-111683902
# 2 ENSMUST00000022124     Cd180 ENSMUSG00000021624 protein_coding exon chr13-102693558-102693789
# query_region distance        p_val avg_log2FC pct.1 pct.2    p_val_adj
# 1 chr12-111693685-111694603     9782 8.753562e-49  -1.531584 0.530 0.678 1.761024e-43
# 2 chr13-102573167-102574085   119472 2.257749e-48  -1.987180 0.306 0.523 4.542094e-43
head(closest_genes_down,2)
# tx_id gene_name            gene_id   gene_biotype type
# ENSMUST00000168338 ENSMUST00000168338   Trmt61a ENSMUSG00000060950 protein_coding  utr
# ENSMUSE00000315744 ENSMUST00000022124     Cd180 ENSMUSG00000021624 protein_coding exon
# closest_region              query_region distance
# ENSMUST00000168338 chr12-111682870-111683902 chr12-111693685-111694603     9782
# ENSMUSE00000315744 chr13-102693558-102693789 chr13-102573167-102574085   119472
dim(closest_genes_down)
# [1] 1032    8
dim(overal_down)
# [1] 1032   13 
write_xlsx(overal_down, "closest_genes_down.xlsx")
#####################

# GO analysis
library(clusterProfiler)
# BiocManager::install("org.Hs.eg.db")
library(org.Mm.eg.db)

library(enrichplot)
head(overal_up)

up_genes_ego <- enrichGO(
  gene = overal_up$gene_id,
  keyType = "ENSEMBL",
  OrgDb = org.Mm.eg.db,
  ont = "BP", # Biological Process ontology
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


barplot(up_genes_ego,showCategory = 20)




down_gene_ego <- enrichGO(
  gene = closest_genes_down$gene_id,
  keyType = "ENSEMBL",
  OrgDb = org.Mm.eg.db,
  ont = "BP", # Biological Process ontology
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


barplot(down_gene_ego,showCategory = 20)
################## UP$DOWN
install.packages("org.Mm.eg.db")
library(clusterProfiler)
library(org.Mm.eg.db)

# Get unique DEGs from both comparisons
deg_genes <- unique(c(rownames(degs_P10), rownames(degs_P17)))

go_results <- enrichGO(gene = deg_genes,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",  # Use "ENTREZID" if genes are in Entrez ID format
                       ont = "BP",          # Biological Process
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

library(enrichplot)
barplot(go_results, 
        showCategory = 20, 
        title = "GO Enrichment - Combined P10 and P17 WT/KO", 
        x = "GeneRatio")  # You can change to "Count" or "GeneRatio"


library(ggplot2)


# Convert results to a data frame
go_df <- as.data.frame(go_results)

# Add a column for NES (here, assigning random positive/negative values for demonstration)
# Replace this with your actual NES or custom grouping
set.seed(123)
go_df$NES <- runif(nrow(go_df), -2, 2)

# Categorize terms based on NES or another logical grouping
go_df$Category <- ifelse(go_df$NES > 0, "Positive", "Negative")

library(ggplot2)
# Select top 10 positive and top 10 negative terms based on NES
top_terms <- rbind(
  go_df[go_df$NES > 0, ][order(go_df$NES, decreasing = TRUE), ][1:10, ],
  go_df[go_df$NES < 0, ][order(go_df$NES), ][1:10, ]
)

library(ggplot2)

ggplot(top_terms, aes(x = reorder(Description, NES), y = NES, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("blue", "green")) +
  labs(x = "", y = "Normalized Enrichment Score (NES)", title = "GO Enrichment - P10 and P17 WT/KO") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),          # Reduce label size
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

############ select neuro-related pathways 
library(dplyr)
colnames(go_df)

neuro_pathways <- go_df %>%
  filter(grepl("neuro|synapse|axon|astrocyte|gliogenesis|nervous system", Description, ignore.case = TRUE))

neuro_pathways <- go_df %>%
  filter(grepl("neuro|synapse|axon", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust) %>%
  slice_head(n =30)

print(neuro_pathways)

library(ggplot2)

ggplot(neuro_pathways, aes(x = reorder(Description, NES), y = NES, fill = NES > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Neuro Pathway", y = "Normalized Enrichment Score (NES)", 
       title = "Neuro-Related GO Enrichment") +
  scale_fill_manual(values = c("blue", "green")) +
  theme_minimal()






###########################
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("stuart-lab/signac", ref = "develop")

# # install.packages("Signac")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("Signac")]
remotes::install_github("timoast/signac")

library(Signac)
packageVersion("Signac")

brain_P10 <- SortIdents(brain_P10)
# Creating pseudobulk profiles for 4 variables across 201178 features
# Computing euclidean distance between pseudobulk profiles
# Clustering distance matrix


# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_up), LookupGeneCoords(brain_P10, "Gfap"))

CoveragePlot(
  object = brain_P10,
  region = "Gfap",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)
