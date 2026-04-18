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
library(tidyr)
library(pheatmap)
library(MAST)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(writexl)
library(enrichplot)
library(enrichR)
# BiocManager::install("rmsigdbr")
# install.packages("msigdbr")
library(clusterProfiler)
library(org.Mm.eg.db)
library(rWikiPathways)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(msigdbr)

setwd("/Users/mortezaabyadeh/Desktop/npc scrna")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(Matrix)
library(hdf5r)
library(pheatmap)

setwd("/Users/mortezaabyadeh/Desktop/npc scrna")

s1 <- Read10X_h5("filtered1_feature_bc_matrix.h5")
s2 <- Read10X_h5("filtered2_feature_bc_matrix.h5")
s3 <- Read10X_h5("filtered3_feature_bc_matrix.h5")
s4 <- Read10X_h5("filtered4_feature_bc_matrix.h5")

s1 <- CreateSeuratObject(s1, project="F_CTRL", min.cells=3, min.features=200)
s2 <- CreateSeuratObject(s2, project="M_CTRL", min.cells=3, min.features=200)
s3 <- CreateSeuratObject(s3, project="F_TREAT", min.cells=3, min.features=200)
s4 <- CreateSeuratObject(s4, project="M_TREAT", min.cells=3, min.features=200)

s1$sex <- "Female"; s1$treatment <- "Control"
s2$sex <- "Male";   s2$treatment <- "Control"
s3$sex <- "Female"; s3$treatment <- "Treated"
s4$sex <- "Male";   s4$treatment <- "Treated"

combined <- merge(s1, y=c(s2,s3,s4),
                  add.cell.ids=c("FCTRL","MCTRL","FTREAT","MTREAT"))

combined$sex_treat <- interaction(combined$sex, combined$treatment)

### QC and normalize
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern="^mt-")

combined <- SCTransform(combined, vars.to.regress="percent.mt", verbose=FALSE)

combined <- RunPCA(combined, verbose=FALSE)
combined <- RunHarmony(combined, group.by.vars="orig.ident")

combined <- FindNeighbors(combined, reduction="harmony", dims=1:30)
combined <- FindClusters(combined, resolution=0.5)
combined <- RunUMAP(combined, reduction="harmony", dims=1:30)

DimPlot(combined, label=TRUE)

### MARKER SETS

markers <- list(
  
  ## --- Neurons ---
  Glutamatergic = c("Slc17a7","Camk2a","Neurod1","Satb2","Tbr1","Cux1","Cux2","Rorb","Bcl11b","Fezf2","Foxp2"),
  
  GABAergic = c("Slc32a1","Gad1","Gad2","Lhx6","Pvalb","Sst","Vip","Lamp5","Nkx2-1","Reln"),
  
  PV = c("Pvalb","Kcnc1","Gad1","Gad2","Slc32a1","Nkx2-1","Sox6","Sp9"),
  
  SST = c("Sst","Npy","Gad1","Nos1","Cort","Htr3a","Calb1","Pcdh7"),
  
  VIP = c("Vip","Htr3a","Cck","Crh","Calb2","Reln","Slc32a1"),
  
  CGE = c("Lamp5","Sncg","Gad2","Id2","Ndnf","Nos1","Reln","Kit"),
  
  Dopaminergic = c("Th","Slc6a3","Nr4a2","Ddc","Lmx1b","Foxa2","Aldh1a1","Slc18a2","En1"),
  
  Serotonergic = c("Tph2","Slc6a4","Fev","Ddc","Pet1","Gata2","Gata3","Slc18a2"),
  
  Cholinergic = c("Chat","Slc5a7","Ache","Chrnb4","Isl1","Nkx2-1","Gbx1"),
  
  Noradrenergic = c("Th","Dbh","Slc6a2","Ddc","Phox2a","Phox2b","Hand2"),
  
  ## --- Cerebellum ---
  Purkinje = c("Calb1","Car8","Pcp2","Itpr1","Rora","Grid2","Nfatc4","Aldoc"),
  
  Cerebellar_Granule = c("Gabra6","Neurod1","Rbfox3","Pax6","Atoh1","Meis2","Gap43","Zic1"),
  
  ## --- Hippocampus ---
  CA1 = c("Fibcd1","Wfs1","Cdh13","Slc17a7","Camk2a","Satb2","Cck","Cux2"),
  
  CA3 = c("Pvrl3","Stxbp6","Chst9","Slc17a7","Meis2","Spp1","Pfkl"),
  
  Dentate = c("Prox1","Dsp","Lct","Neurod1","C1ql2","Gabra2","Slc17a7"),
  
  ## --- Hypothalamus ---
  Hypothalamic = c("Agrp","Pomc","Hcrt","Npy","Oxt","Avp","Pmch","Nms","Trh","Ghrh","Crh"),
  
  ## --- Striatum ---
  MSN_D1 = c("Drd1","Pdyn","Chrm4","Gad1","Bcl11b","Foxp1","Ppp1r1b"),
  
  MSN_D2 = c("Drd2","Penk","Gpr52","Gad1","Bcl11b","Foxp1","Ppp1r1b","Adora2a"),
  
  ## --- Glia ---
  Astrocytes = c("Aldh1l1","Gfap","Aqp4","S100b","Sox9","Slc1a3","Glul","Fgfr3","Acsl6","Gjb6"),
  
  OPC = c("Pdgfra","Cspg4","Olig2","Sox10","Nkx2-2","Ascl1","Gpr17","Smc4"),
  
  Oligodendrocytes = c("Mbp","Plp1","Mog","Mag","Cnp","Cldn11","Ugt8a","Sox10","Nkx6-2"),
  
  NewlyFormed_Oligo = c("Tcf7l2","Enpp6","Olig1","Olig2","Sox10","Tnr","Gjb6"),
  
  ## --- Microglia ---
  Microglia = c("Cx3cr1","P2ry12","Tmem119","Csf1r","Hexb","Fcrls","Trem2","Sall1","Adgrg1"),
  
  Activated_Microglia = c("Trem2","Apoe","Lpl","Axl","Clec7a","Itgax","Cst7","Cd9","Lgals3"),
  
  ## --- Other brain cells ---
  Endothelial = c("Cldn5","Pecam1","Flt1","Tie1","Tek","Emcn","Sox17","Sox18","Esam","Cd34"),
  
  Pericytes = c("Pdgfrb","Rgs5","Kcnj8","Notch3","Cd248","Anpep","Cspg4","Abcc9"),
  
  SmoothMuscle = c("Acta2","Tagln","Myh11","Cnn1","Des","Mustn1","Myl9","Notch3"),
  
  BAM = c("Mrc1","Lyve1","F13a1","Ms4a7","Pf4","Cx3cr1","Cbr2","Folr2"),
  
  Mast = c("Cpa3","Hdc","Kit","Mcpt4","Fcer1a","Ms4a2","Slc45a3"),
  
  Ependymal = c("Foxj1","Ccdc153","Cfap54","Sox9","Vim","S100b","Ezr","Enkur","Dydc1"),
  
  ## --- Progenitors ---
  RadialGlia = c("Pax6","Sox2","Nes","Vim","Slc1a3","Fabp7","Hes1","Hes5"),
  
  IPC = c("Eomes","Svet1","Ngn2","Insm1","Tgif2","Neurog1","Gadd45g"),
  
  NSC = c("Nestin","Gfap","Sox2","Ascl1","Mcm2","Ki67","Dcx","Egfr","Id4"),
  
  Neuroblasts = c("Dcx","Psancam","Stmn1","Tubb3","Neurod2","Nrxn3","Pak3")
)


### MODULE SCORE

combined <- AddModuleScore(
  combined,
  features = markers,
  name = names(markers)
)

### EXTRACT SCORES
colnames(combined@meta.data)

celltype_names <- c(
  "Glutamatergic","GABAergic","PV","SST","VIP","CGE",
  "Dopaminergic","Serotonergic","Cholinergic","Noradrenergic",
  "Purkinje","Cerebellar_Granule",
  "CA1","CA3","Dentate","Hypothalamic",
  "MSN_D1","MSN_D2",
  "Astrocytes","OPC","Oligodendrocytes","NewlyFormed_Oligo",
  "Microglia","Activated_Microglia",
  "Endothelial","Pericytes","SmoothMuscle",
  "BAM","Mast","Ependymal",
  "RadialGlia","IPC","NSC","Neuroblasts"
)
score_cols <- unlist(lapply(celltype_names, function(ct) {
  grep(paste0("^", ct), colnames(combined@meta.data), value = TRUE)
}))



meta <- combined@meta.data

score_mat <- sapply(celltype_names, function(ct) {
  cols <- grep(paste0("^", ct), colnames(meta), value = TRUE)
  
  if (length(cols) == 0) {
    return(rep(NA, nrow(meta)))  # important
  } else if (length(cols) == 1) {
    return(meta[, cols])
  } else {
    return(apply(meta[, cols], 1, max, na.rm = TRUE))
  }
})


combined$manual_celltype <- apply(score_mat, 1, function(x) {
  
  if (all(is.na(x))) {
    return("Unknown")
  }
  
  max_val <- max(x, na.rm = TRUE)
  
  if (max_val < 0.05) {
    return("Unknown")
  } else {
    return(names(x)[which.max(x)])
  }
})


library(ggplot2)
library(RColorBrewer)

n <- length(unique(combined$manual_celltype))

colors <- colorRampPalette(brewer.pal(12, "Paired"))(n)

DimPlot(
  combined,
  group.by = "manual_celltype",
  split.by="sex_treat",
  #label = TRUE,
  repel = TRUE,
  label.size = 5
)






### ASSIGN CELL TYPE (WINNER TAKES ALL)

combined$manual_celltype <- colnames(score_mat)[apply(score_mat, 1, which.max)]
combined$manual_celltype <- gsub("\\d+", "", combined$manual_celltype)


### VISUALIZATION

DimPlot(combined, group.by="manual_celltype", label=TRUE)

DimPlot(combined, group.by="manual_celltype", split.by="sex_treat")


# QC CHECK MARKERS




# CELL TYPE COMPOSITION

df <- as.data.frame(table(
  CellType = combined$manual_celltype,
  Group = combined$sex_treat
))

colnames(df) <- c("CellType","Group","Freq")

df <- df %>%
  group_by(Group) %>%
  mutate(Percent = Freq / sum(Freq) * 100)

ggplot(df, aes(x=Group, y=Percent, fill=CellType)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))


df <- as.data.frame(table(
  CellType = combined$manual_celltype,
  Group = combined$sex_treat
))

colnames(df) <- c("CellType","Group","Freq")

df <- df %>%
  group_by(Group) %>%
  mutate(Percent = Freq / sum(Freq) * 100)

# Convert to wide format
df_wide <- df %>%
  select(CellType, Group, Percent) %>%
  pivot_wider(names_from = Group, values_from = Percent)

df_wide

library(writexl)

write_xlsx(df_wide, "celltype_percentages.xlsx")


# HEATMAP (optional)


df <- as.data.frame(table(
  CellType = combined$manual_celltype,
  Group = combined$sex_treat
))

colnames(df) <- c("CellType","Group","Freq")

df <- df %>%
  group_by(Group) %>%
  mutate(Percent = Freq / sum(Freq) * 100)
df_wide <- df %>%
  select(Group, CellType, Percent) %>%
  tidyr::pivot_wider(names_from = CellType,
                     values_from = Percent,
                     values_fill = 0)

mat <- as.data.frame(df_wide)
rownames(mat) <- mat$Group
mat$Group <- NULL
mat <- as.matrix(mat)

mat_scaled <- scale(mat)

library(pheatmap)

pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 9,
  main = "Cell Composition (Z-score scaled)"
)

############################################# so here we may loose some cells due to low number of identified; tehrefore we wont
######## have them in degs analysis

library(tidyverse)

df <- read.table(text = "
CellType Female_Control Female_Treated Male_Control Male_Treated
Activated_Microglia 9 3 4 3
Astrocytes 228 200 200 272
BAM 5 9 4 7
CA 370 266 334 260
Cerebellar_Granule 830 872 947 921
CGE 20 30 22 26
Cholinergic 2 3 4 5
Dentate 19 30 19 29
Dopaminergic 3 6 10 15
Endothelial 50 52 63 89
Ependymal 9 9 8 8
GABAergic 1 0 1 1
Glutamatergic 35 43 40 46
Hypothalamic 9 0 0 0
IPC 1 0 1 0
Mast 8 10 10 12
Microglia 80 89 48 102
MSN_D 249 214 205 248
Neuroblasts 3322 3666 3006 4051
NewlyFormed_Oligo 52 69 65 81
Noradrenergic 1 0 0 0
NSC 14 10 12 9
Oligodendrocytes 964 1067 632 971
OPC 78 39 55 43
Pericytes 31 22 22 17
Purkinje 2761 2945 2398 3310
PV 25 23 34 31
RadialGlia 6 1 3 4
Serotonergic 14 7 3 6
SmoothMuscle 2 5 3 3
SST 41 29 29 38
VIP 3 27 22 27
", header = TRUE)

df_long <- df %>%
  pivot_longer(-CellType, names_to = "Group", values_to = "Count")


ggplot(df_long, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~CellType, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )



############### DEGs ###############

combined <- PrepSCTFindMarkers(combined)


combined$group <- paste(combined$sex, combined$treatment, sep = "_")
Idents(combined) <- "group"
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)



#### cell number in each cells ##### if number of cells are low the deg os not statistically significant and not have the power
table(combined$manual_celltype, combined$group)

cell_count <- as.data.frame.matrix(
  table(combined$manual_celltype, combined$group)
)

cell_count$CellType <- rownames(cell_count)


cell_count <- cell_count[, c("CellType", setdiff(colnames(cell_count), "CellType"))]

write_xlsx(cell_count, "CellType_CellCounts.xlsx")




run_degs <- function(obj, ct, g1, g2) {
  
  # check if cell type exists
  if (!(ct %in% obj$manual_celltype)) {
    return(NULL)
  }
  
  tmp <- subset(obj, subset = manual_celltype == ct)
  
  # skip empty or too small
  if (is.null(tmp) || ncol(tmp) < 30) return(NULL)
  
  Idents(tmp) <- "group"
  
  # check groups exist
  if (!all(c(g1, g2) %in% unique(tmp$group))) return(NULL)
  
  deg <- FindMarkers(
    tmp,
    ident.1 = g1,
    ident.2 = g2,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  
  if (nrow(deg) == 0) return(NULL)
  
  deg$gene <- rownames(deg)
  deg$log2FC <- deg$avg_log2FC
  deg$FC <- 2^(deg$avg_log2FC)
  deg$pct_expr <- (deg$pct.1 + deg$pct.2) / 2
  
  deg <- deg %>%
    dplyr::filter(p_val_adj < 0.05, abs(log2FC) > 0.25)
  
  return(deg)
}


male_deg_list <- list()

for (ct in celltype_names) {
  cat("Male:", ct, "\n")
  
  res <- run_degs(
    combined,
    ct = ct,
    g1 = "Male_Treated",
    g2 = "Male_Control"
  )
  
  if (!is.null(res)) {
    male_deg_list[[ct]] <- res
  }
}

female_deg_list <- list()

for (ct in celltype_names) {
  cat("Female:", ct, "\n")
  
  res <- run_degs(
    combined,
    ct = ct,
    g1 = "Female_Treated",
    g2 = "Female_Control"
  )
  
  if (!is.null(res)) {
    female_deg_list[[ct]] <- res
  }
}



male_deg_all <- bind_rows(male_deg_list, .id = "celltype")
female_deg_all <- bind_rows(female_deg_list, .id = "celltype")


write_xlsx(
  list(
    Male_Treated_vs_Control = male_deg_all,
    Female_Treated_vs_Control = female_deg_all
  ),
  path = "DEGs_by_celltype_sex_treatment.xlsx"
)


########### all genes ###########
run_all_genes <- function(obj, ct, g1, g2) {
  
  if (!(ct %in% obj$manual_celltype)) return(NULL)
  
  tmp <- subset(obj, subset = manual_celltype == ct)
  
  if (is.null(tmp) || ncol(tmp) < 30) return(NULL)
  
  Idents(tmp) <- "group"
  
  if (!all(c(g1, g2) %in% unique(tmp$group))) return(NULL)
  
  res <- FindMarkers(
    tmp,
    ident.1 = g1,
    ident.2 = g2,
    test.use = "wilcox",
    logfc.threshold = 0,   # KEEP ALL GENES
    min.pct = 0            # KEEP ALL GENES
  )
  
  if (nrow(res) == 0) return(NULL)
  
  res$gene <- rownames(res)
  res$log2FC <- res$avg_log2FC
  res$FC <- 2^(res$avg_log2FC)
  res$pct_expr <- (res$pct.1 + res$pct.2) / 2
  
  # reorder columns nicely
  res <- res[, c("gene", "log2FC", "FC", "p_val", "p_val_adj", "pct.1", "pct.2", "pct_expr")]
  
  return(res)
}


male_sheets <- list()

for (ct in celltype_names) {
  cat("Male:", ct, "\n")
  
  res <- run_all_genes(
    combined,
    ct = ct,
    g1 = "Male_Treated",
    g2 = "Male_Control"
  )
  
  if (!is.null(res)) {
    male_sheets[[ct]] <- res
  }
}


female_sheets <- list()

for (ct in celltype_names) {
  cat("Female:", ct, "\n")
  
  res <- run_all_genes(
    combined,
    ct = ct,
    g1 = "Female_Treated",
    g2 = "Female_Control"
  )
  
  if (!is.null(res)) {
    female_sheets[[ct]] <- res
  }
}




write_xlsx(male_sheets, "Male_AllGenes_by_CellType.xlsx")
write_xlsx(female_sheets, "Female_AllGenes_by_CellType.xlsx")


####### volcano 

plot_volcano <- function(df, ct, sex_label) {
  
  df <- df %>%
    mutate(
      significance = case_when(
        p_val_adj < 0.05 & log2FC > 0.25 ~ "Up",
        p_val_adj < 0.05 & log2FC < -0.25 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(df, aes(x = log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    theme_classic() +
    ggtitle(paste(sex_label, "-", ct)) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey"))
  
  ggsave(paste0("Volcano_", sex_label, "_", ct, ".png"), p, width = 6, height = 5)
}

for (ct in names(male_sheets)) {
  plot_volcano(male_sheets[[ct]], ct, "Male")
}

for (ct in names(female_sheets)) {
  plot_volcano(female_sheets[[ct]], ct, "Female")
}




rank_celltypes <- function(sheet_list) {
  
  data.frame(
    CellType = names(sheet_list),
    DEG_Count = sapply(sheet_list, function(df) {
      sum(df$p_val_adj < 0.05 & abs(df$log2FC) > 0.25)
    })
  ) %>%
    arrange(desc(DEG_Count))
}

male_rank <- rank_celltypes(male_sheets)
female_rank <- rank_celltypes(female_sheets)

write_xlsx(
  list(Male = male_rank, Female = female_rank),
  "CellType_Response_Ranking.xlsx"
)


######### patwhaya analysis


write_xlsx(male_sheets, "Male_AllGenes_by_CellType.xlsx")
write_xlsx(female_sheets, "Female_AllGenes_by_CellType.xlsx")
setwd("/Users/mortezaabyadeh/Desktop/npc scrna/new")
library(readxl)
male_sheets <- lapply(excel_sheets("Male_AllGenes_by_CellType.xlsx"), function(sh) {read_excel("Male_AllGenes_by_CellType.xlsx", sheet = sh)
})

names(male_sheets) <- excel_sheets("Male_AllGenes_by_CellType.xlsx")

female_sheets <- lapply(excel_sheets("Female_AllGenes_by_CellType.xlsx"), function(sh) {
  read_excel("Female_AllGenes_by_CellType.xlsx", sheet = sh)
})

names(female_sheets) <- excel_sheets("Female_AllGenes_by_CellType.xlsx")

count_degs <- function(sheet_list, sex_label) {
  
  df <- data.frame(
    CellType = names(sheet_list),
    Up = sapply(sheet_list, function(x) {
      sum(x$p_val_adj < 0.05 & x$log2FC > 0)
    }),
    Down = sapply(sheet_list, function(x) {
      sum(x$p_val_adj < 0.05 & x$log2FC < 0)
    })
  )
  
  df$Total <- df$Up + df$Down
  df$Sex <- sex_label
  
  return(df)
}


male_deg_counts <- count_degs(male_sheets, "Male")
female_deg_counts <- count_degs(female_sheets, "Female")

deg_counts <- rbind(male_deg_counts, female_deg_counts)


ggplot(deg_counts, aes(x = reorder(CellType, Total), y = Total, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_classic() +
  ylab("Number of DEGs (adj p < 0.05)") +
  xlab("Cell Type") +
  ggtitle("DEG Count per Cell Type") +
  scale_fill_manual(values = c("Male" = "steelblue", "Female" = "tomato"))




deg_long <- deg_counts %>%
  tidyr::pivot_longer(cols = c("Up", "Down"),
                      names_to = "Direction",
                      values_to = "Count")

ggplot(deg_long, aes(x = CellType, y = Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Sex) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  ggtitle("Up vs Down regulated genes")



split_genes <- function(df) {
  
  up_genes <- df %>%
    filter(p_val_adj < 0.05 & log2FC > 0) %>%
    pull(gene)
  
  down_genes <- df %>%
    filter(p_val_adj < 0.05 & log2FC < 0) %>%
    pull(gene)
  
  return(list(up = up_genes, down = down_genes))
}



run_enrich <- function(genes) {
  
  if (length(genes) < 10) return(NULL)
  
  gene_df <- bitr(genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  
  enrich <- enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pvalueCutoff = 0.05
  )
  
  return(enrich@result)
}



ora_results <- list()

for (sex in c("Male", "Female")) {
  
  sheets <- if (sex == "Male") male_sheets else female_sheets
  
  for (ct in names(sheets)) {
    
    df <- sheets[[ct]]
    
    genes <- split_genes(df)
    
    up_res <- run_enrich(genes$up)
    down_res <- run_enrich(genes$down)
    
    if (!is.null(up_res)) {
      ora_results[[paste(sex, ct, "UP", sep = "_")]] <- up_res
    }
    
    if (!is.null(down_res)) {
      ora_results[[paste(sex, ct, "DOWN", sep = "_")]] <- down_res
    }
  }
}


write_xlsx(ora_results, "ORA_Up_Down_Pathways.xlsx")








gsea_results <- list()

for (sex in c("Male", "Female")) {
  
  sheets <- if (sex == "Male") male_sheets else female_sheets
  
  for (ct in names(sheets)) {
    
    df <- sheets[[ct]]
    
    gene_list <- df$log2FC
    names(gene_list) <- df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    gene_df <- bitr(names(gene_list),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
    
    gene_list <- gene_list[gene_df$SYMBOL]
    names(gene_list) <- gene_df$ENTREZID
    
    gsea <- tryCatch({
      gseGO(
        geneList = gene_list,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pvalueCutoff = 0.05,
        verbose = FALSE
      )
    }, error = function(e) NULL)
    
    if (!is.null(gsea) && nrow(gsea@result) > 0) {
      gsea_results[[paste(sex, ct, sep = "_")]] <- gsea@result
    }
  }
}

write_xlsx(gsea_results, "GSEA_All_CellTypes.xlsx")








####### wikipath mouse 2024

run_wp_enrich <- function(genes) {
  
  if (length(genes) < 10) return(NULL)
  
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  enrichWP(
    gene = gene_df$ENTREZID,
    organism = "Mus musculus",
    pvalueCutoff = 0.05
  )
}

get_up_down <- function(df) {
  
  up <- df %>%
    filter(p_val_adj < 0.05 & log2FC > 0) %>%
    pull(gene)
  
  down <- df %>%
    filter(p_val_adj < 0.05 & log2FC < 0) %>%
    pull(gene)
  
  list(up = up, down = down)
}

plot_wp_dotplot <- function(enrich_obj, ct, sex, direction) {
  
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) return(NULL)
  
  p <- dotplot(enrich_obj, showCategory = 5) +
    ggtitle(paste(sex, "-", ct, "-", direction, "WikiPathways"))
  
  ggsave(
    paste0("WP_", sex, "_", ct, "_", direction, ".png"),
    p,
    width = 6,
    height = 5
  )
}


wp_results <- list()

for (sex in c("Male", "Female")) {
  
  sheets <- if (sex == "Male") male_sheets else female_sheets
  
  for (ct in names(sheets)) {
    
    cat(sex, "-", ct, "\n")
    
    df <- sheets[[ct]]
    
    genes <- get_up_down(df)
    
    # ---- UP ----
    wp_up <- run_wp_enrich(genes$up)
    
    if (!is.null(wp_up)) {
      name <- paste(sex, ct, "UP", sep = "_")
      wp_results[[name]] <- wp_up@result
      
      plot_wp_dotplot(wp_up, ct, sex, "UP")
    }
    
    # ---- DOWN ----
    wp_down <- run_wp_enrich(genes$down)
    
    if (!is.null(wp_down)) {
      name <- paste(sex, ct, "DOWN", sep = "_")
      wp_results[[name]] <- wp_down@result
      
      plot_wp_dotplot(wp_down, ct, sex, "DOWN")
    }
  }
}


library(writexl)

write_xlsx(wp_results, "WikiPathways_Up_Down_All_CellTypes.xlsx")


###### KEGG 2026

run_kegg <- function(genes) {
  
  if (length(genes) < 10) return(NULL)
  
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  ekegg <- enrichKEGG(
    gene         = gene_df$ENTREZID,
    organism     = "mmu",
    pvalueCutoff = 0.05
  )
  
  if (is.null(ekegg) || nrow(ekegg@result) == 0) return(NULL)
  
  return(ekegg)
}


get_up_down <- function(df) {
  
  up <- df %>%
    filter(p_val_adj < 0.05 & log2FC > 0) %>%
    pull(gene)
  
  down <- df %>%
    filter(p_val_adj < 0.05 & log2FC < 0) %>%
    pull(gene)
  
  list(up = up, down = down)
}




plot_kegg_dotplot <- function(kegg_obj, ct, sex, direction) {
  
  if (is.null(kegg_obj)) return(NULL)
  
  p <- dotplot(kegg_obj, showCategory = 5) +
    ggtitle(paste(sex, ct, direction, "KEGG 2026"))
  
  ggsave(
    filename = paste0("KEGG_", sex, "_", ct, "_", direction, ".png"),
    plot = p,
    width = 6,
    height = 5
  )
}





kegg_results <- list()

for (sex in c("Male", "Female")) {
  
  sheets <- if (sex == "Male") male_sheets else female_sheets
  
  for (ct in names(sheets)) {
    
    cat("KEGG:", sex, "-", ct, "\n")
    
    df <- sheets[[ct]]
    
    genes <- get_up_down(df)
    
    # ---------------- UP ----------------
    kegg_up <- run_kegg(genes$up)
    
    if (!is.null(kegg_up)) {
      name <- paste(sex, ct, "UP", sep = "_")
      kegg_results[[name]] <- kegg_up@result
      
      plot_kegg_dotplot(kegg_up, ct, sex, "UP")
    }
    
    # ---------------- DOWN ----------------
    kegg_down <- run_kegg(genes$down)
    
    if (!is.null(kegg_down)) {
      name <- paste(sex, ct, "DOWN", sep = "_")
      kegg_results[[name]] <- kegg_down@result
      
      plot_kegg_dotplot(kegg_down, ct, sex, "DOWN")
    }
  }
}



write_xlsx(kegg_results, "KEGG_2026_Up_Down_By_CellType.xlsx")



######### wiki pathyaw mouse 2024

run_wikipathways <- function(genes) {
  
  if (length(genes) < 10) return(NULL)
  
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  if (is.null(gene_df) || nrow(gene_df) < 10) return(NULL)
  
  wp <- enrichWP(
    gene         = gene_df$ENTREZID,
    organism     = "Mus musculus",
    pvalueCutoff = 0.05
  )
  
  if (is.null(wp) || nrow(wp@result) == 0) return(NULL)
  
  return(wp)
}


get_up_down <- function(df) {
  
  df$log2FC <- df$log2FC
  
  up <- df %>%
    filter(p_val_adj < 0.05 & log2FC > 0) %>%
    pull(gene)
  
  down <- df %>%
    filter(p_val_adj < 0.05 & log2FC < 0) %>%
    pull(gene)
  
  list(up = up, down = down)
}

plot_wp_combined <- function(wp_up, wp_down, ct, sex) {
  
  up_df <- if (!is.null(wp_up)) {
    wp_up@result %>%
      mutate(Direction = "Up")
  } else NULL
  
  down_df <- if (!is.null(wp_down)) {
    wp_down@result %>%
      mutate(Direction = "Down")
  } else NULL
  
  df <- bind_rows(up_df, down_df)
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    arrange(p.adjust) %>%
    group_by(Direction) %>%
    slice_head(n = 5) %>%   # top 5 per direction
    ungroup() %>%
    mutate(logP = -log10(p.adjust))
  
  p <- ggplot(df, aes(x = GeneRatio, y = logP)) +
    
    geom_point(aes(size = Count, color = Direction), alpha = 0.85) +
    
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    
    geom_text(aes(label = Description),
              size = 3,
              hjust = 0,
              nudge_x = 0.02,
              show.legend = FALSE) +
    
    theme_classic() +
    
    labs(
      title = paste(sex, ct, "WikiPathways 2024 Mouse"),
      x = "Gene Ratio (Fold Enrichment)",
      y = "-log10 adjusted p-value",
      size = "Gene Count"
    )
  
  ggsave(
    filename = paste0("WP_", sex, "_", ct, "_COMBINED.png"),
    plot = p,
    width = 8,
    height = 5
  )
  
  return(p)
}


for (sex in c("Male", "Female")) {
  
  sheets <- if (sex == "Male") male_sheets else female_sheets
  
  for (ct in names(sheets)) {
    
    cat("WP combined:", sex, "-", ct, "\n")
    
    df <- sheets[[ct]]
    if (is.null(df)) next
    
    genes <- get_up_down(df)
    
    wp_up <- run_wikipathways(genes$up)
    wp_down <- run_wikipathways(genes$down)
    
    plot_wp_combined(wp_up, wp_down, ct, sex)
  }
}


######## facet for patwhays below 0.05 #####

collect_wikipathways <- function(sheets, sex) {
  
  all_results <- list()
  
  for (ct in names(sheets)) {
    
    df <- sheets[[ct]]
    if (is.null(df)) next
    
    genes <- get_up_down(df)
    
    wp_up <- run_wikipathways(genes$up)
    wp_down <- run_wikipathways(genes$down)
    
    # UP
    if (!is.null(wp_up)) {
      up_df <- wp_up@result %>%
        mutate(
          CellType = ct,
          Direction = "Up"
        )
      all_results[[length(all_results)+1]] <- up_df
    }
    
    # DOWN
    if (!is.null(wp_down)) {
      down_df <- wp_down@result %>%
        mutate(
          CellType = ct,
          Direction = "Down"
        )
      all_results[[length(all_results)+1]] <- down_df
    }
  }
  
  df <- bind_rows(all_results)
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # filter ONLY significant
  df <- df %>%
    filter(pvalue < 0.05) %>%
    mutate(
      Pathway = Description,
      logFC = GeneRatio
    )
  
  return(df)
}
male_wp <- collect_wikipathways(male_sheets, "Male")
female_wp <- collect_wikipathways(female_sheets, "Female")

head(male_wp)
plot_megawp <- function(df, sex) {
  
  df$Pathway <- reorder(df$Pathway, df$FoldEnrichment)
  
  p <- ggplot(df, aes(x = FoldEnrichment, y = Pathway)) +
    
    geom_point(aes(size = Count, color = Direction), alpha = 0.8) +
    
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    
    facet_wrap(~CellType, scales = "free_y") +
    
    theme_classic() +
    
    labs(
      title = paste(sex, "WikiPathways 2024 Mouse (P.value < 0.05)"),
      x = "Fold Enrichment",
      y = "Pathway",
      size = "Gene Count"
    )
  
  ggsave(paste0("MEGA_WP_", sex, ".png"), p, width = 14, height = 10)
  
  return(p)
}

plot_megawp(male_wp, "Male")
plot_megawp(female_wp, "Female")

colnames(wp@result)


##################################################### msig 2020



run_hallmark <- function(genes, db = "MSigDB_Hallmark_2020") {
  
  genes <- unique(na.omit(genes))
  
  if (length(genes) < 10) return(NULL)
  
  res <- enrichr(genes, c(db))
  
  df <- res[[db]]
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    filter(P.value < 0.05)
  
  return(df)
}



get_up_down <- function(df) {
  
  up <- df %>%
    filter(p_val_adj < 0.05 & log2FC > 0) %>%
    pull(gene)
  
  down <- df %>%
    filter(p_val_adj < 0.05 & log2FC < 0) %>%
    pull(gene)
  
  list(up = up, down = down)
}

collect_hallmark <- function(sheets, sex) {
  
  all_results <- list()
  
  for (ct in names(sheets)) {
    
    df <- sheets[[ct]]
    if (is.null(df)) next
    
    genes <- get_up_down(df)
    
    up <- run_hallmark(genes$up)
    down <- run_hallmark(genes$down)
    
    # ---------------- UP ----------------
    if (!is.null(up)) {
      up <- up %>%
        mutate(
          CellType = ct,
          Direction = "Up"
        )
      all_results[[length(all_results) + 1]] <- up
    }
    
    # ---------------- DOWN ----------------
    if (!is.null(down)) {
      down <- down %>%
        mutate(
          CellType = ct,
          Direction = "Down"
        )
      all_results[[length(all_results) + 1]] <- down
    }
  }
  
  df <- bind_rows(all_results)
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- df %>%
    filter(P.value < 0.05)
  
  # convert Enrichr score → fold enrichment proxy
  df <- df %>%
    mutate(
      FoldEnrichment = Combined.Score / max(Combined.Score)
    )
  
  return(df)
}


plot_megahallmark <- function(df, sex) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df$Term <- reorder(df$Term, df$FoldEnrichment)
  
  # keep top pathways per cell type (optional but recommended)
  df <- df %>%
    group_by(CellType, Direction) %>%
    slice_max(FoldEnrichment, n = 5) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = FoldEnrichment, y = Term)) +
    
    geom_point(
      aes(
        size = as.numeric(sub("/.*", "", Overlap)),  # ✅ extract count
        color = Direction
      ),
      alpha = 0.8
    ) +
    
    scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
    
    facet_wrap(~CellType, scales = "free_y") +
    
    theme_classic() +
    
    labs(
      title = paste(sex, "MSigDB Hallmark (P.value < 0.05)"),
      x = "Fold Enrichment",
      y = "Hallmark Pathway",
      size = "Gene Count"
    )
  
  ggsave(
    paste0("MEGA_HALLMARK_", sex, ".png"),
    p,
    width = 14,
    height = 10
  )
  
  return(p)
}


male_hallmark <- collect_hallmark(male_sheets, "Male")
female_hallmark <- collect_hallmark(female_sheets, "Female")

head(male_hallmark)

plot_megahallmark(male_hallmark, "Male")
plot_megahallmark(female_hallmark, "Female")
