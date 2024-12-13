# Load necessary libraries
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  # For mouse data; replace with org.Hs.eg.db for human data
library(fgsea)
library(ggplot2)
library(dplyr)

setwd("Users/mortezaabyadeh/Desktop")
dir()
BG_PROG<-readRDS("BG_PROG.rds")
# Define identities based on sample types, e.g., P10_WT, P10_KO
levels(Idents(BG_PROG))

# Find DEGs between P10_WT and P10_KO
degs_P10 <- FindMarkers(BG_PROG, ident.1 = "BG_P10_WT", ident.2 = "BG_P10_cKO", min.pct = 0.25, logfc.threshold = 0.25)
degs_P17 <- FindMarkers(BG_PROG, ident.1 = "BG_P17_WT" , ident.2 = "BG_P17_cKO", min.pct = 0.25, logfc.threshold = 0.25)


BiocManager::install("clusterProfiler")

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




