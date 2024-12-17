#memory.limit(size = 50000) # Adjust value depending on your system (e.g., 50 GB)
# Example

#rm(list = ls())  # Clear environment
#gc()             # Force garbage collection
library(dplyr)
#sapply(ls(), function(x) object.size(get(x))) %>% sort(decreasing = TRUE)

R.version.string
library(installr)

BiocManager :: install("installr")

updateR()

setRepositories(ind=1:3)
install.packages("Signac")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")




library(Rsamtools)
library(Signac)
library(Seurat)
# Install JASPAR2020 if you haven't already
install.packages("BiocManager")
BiocManager::install("JASPAR2020")

# Load the JASPAR2020 package
library(JASPAR2020)

install.packages("TFBSTools")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("patchwork")
BiocManager::install("MotifDb")
BiocManager::install("motifmatchr")
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(enrichplot)
library(MotifDb)


system.file()

getwd()
setwd("/Users/mortezaabyadeh/Desktop")
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
#brain_P10<- subset(brain_P10, cells = Cells(brain_P10)[1:3000])
#table(brain_P10$predicted.id)
# BG_P10_cKO      BG_P10_WT BGPROG_P10_cKO  BGPROG_P10_WT   MYOC_P10_cKO 
# 1290           1236            156            268             57 
# MYOC_P10_WT  OLIG2_P10_cKO   OLIG2_P10_WT   VAFA_P10_cKO    VAFA_P10_WT 
# 161             69            581            320               455

Assays(brain_P10)
DefaultAssay(brain_P10) <- "peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
# The getMatrixSet() function retrieves a set 
# of Position Frequency Matrices (PFMs) from the JASPAR database.
# Position Frequency Matrix (PFM)
# A PFM represents the frequency of each nucleotide (A, C, G, T) at each position in the transcription factor's binding motif. 
# These matrices are used in motif analysis to identify potential binding sites in DNA sequences.

brain_P10 <- AddMotifs(
  object = brain_P10,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

class(brain_P10[["peaks"]])
levels(Idents(brain_P10))
Idents(brain_P10) <- brain_P10$predicted.id


da_peaks <- FindMarkers(
  object = brain_P10,
  ident.1 = "BG_P10_WT",
  ident.2 =  "BG_P10_cKO",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
# 'LR' stands for Likelihood Ratio Test. It is often used in single-cell 
# RNA-seq and ATAC-seq analyses to test for differences in feature accessibility or expression.
# min.pct = 0.05, A threshold for filtering features (peaks)


summary(da_peaks)

dim(da_peaks)
# [1] 31345     5
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$pct.1 > 0.2, ])
length(top.da.peak)
# [1]  4285

## Matching GC.percent distribution
# peaks.matched can then be used as the background peak set by setting background=peaks.matched in FindMotifs().
# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(brain_P10, idents = c("BG_P10_WT",  "BG_P10_cKO"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(brain_P10, assay = "peaks", layer = "meta.features")

summary(meta.feature[top.da.peak, ])
# we have NAs and I removed them
meta.feature.clean <- na.omit(meta.feature[top.da.peak, ])

# you can impute it; here i used simply mean
# meta.feature[top.da.peak, "GC.percent"][is.na(meta.feature[top.da.peak, "GC.percent"])] <- mean(meta.feature[top.da.peak, "GC.percent"], na.rm = TRUE)


peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature.clean,
  n = 50000
)



#enriched.motifs <- FindMotifs(object = brain_P10,features = top.da.peak)



# test enrichment
enriched.motifs <- FindMotifs(
  object = brain_P10,
  features = meta.feature
)
