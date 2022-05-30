# 4.Subset

table(scc_integrated$cell_type)

# Hepo     Cho     HSC    Endo Kupffer     LCM    pDCs  Neutro B lymph T lymph      NK 
# 8197     457     273   35179    9018    6822    1314    9734   12873    9082    3006

## 1.Hepo
Hepo <- subset(scc_integrated,idents = "Hepo")
saveRDS(Hepo,file="0.Data/Hepo.rds")
Hepo@active.assay <- "SCT"
Hepo@active.ident <- as.factor(Hepo$type)
DimPlot(Hepo)
Hepo_4vC_DEG <- FindMarkers(Hepo,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo_4vC_DEG$change = ifelse(Hepo_4vC_DEG$p_val_adj < 0.05 & abs(Hepo_4vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hepo_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo_4vC_DEG$change)
# Down Stable     Up 
# 63     18    196

Hepo_8vC_DEG <- FindMarkers(Hepo,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo_8vC_DEG$change = ifelse(Hepo_8vC_DEG$p_val_adj < 0.05 & abs(Hepo_8vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hepo_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo_8vC_DEG$change)
# Down Stable     Up 
# 160     39   1140

## 2.Cho
Cho <- subset(scc_integrated,idents = "Cho")
saveRDS(Cho,file="0.Data/Cho.rds")
Cho@active.assay <- "SCT"
Cho@active.ident <- as.factor(Cho$type)
DimPlot(Cho)
Cho_4vC_DEG <- FindMarkers(Cho,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Cho_4vC_DEG$change = ifelse(Cho_4vC_DEG$p_val_adj < 0.05 & abs(Cho_4vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(Cho_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Cho_4vC_DEG$change)
# Down Stable     Up 
# 4    171     19

Cho_8vC_DEG <- FindMarkers(Cho,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Cho_8vC_DEG$change = ifelse(Cho_8vC_DEG$p_val_adj < 0.05 & abs(Cho_8vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(Cho_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Cho_8vC_DEG$change)
# Down Stable     Up 
# 105    651    409

## 3.HSC
HSC <- subset(scc_integrated,idents = "HSC")
saveRDS(HSC,file="0.Data/HSC.rds")
HSC@active.assay <- "SCT"
HSC@active.ident <- as.factor(HSC$type)
DimPlot(HSC)
HSC_4vC_DEG <- FindMarkers(HSC,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
HSC_4vC_DEG$change = ifelse(HSC_4vC_DEG$p_val_adj < 0.05 & abs(HSC_4vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(HSC_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HSC_4vC_DEG$change)
# Down Stable     Up 
# 1    419      1

HSC_8vC_DEG <- FindMarkers(HSC,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
HSC_8vC_DEG$change = ifelse(HSC_8vC_DEG$p_val_adj < 0.05 & abs(HSC_8vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(HSC_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HSC_8vC_DEG$change)
# Down Stable     Up 
# 15   1028    257 
Cho_8vC_DEG <- FindMarkers(Cho,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Cho_8vC_DEG$change = ifelse(Cho_8vC_DEG$p_val_adj < 0.05 & abs(Cho_8vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(Cho_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Cho_8vC_DEG$change)
# Down Stable     Up 
# 105    651    409 

## 4.Endo
Endo <- subset(scc_integrated,idents = "Endo")
saveRDS(Endo,file="0.Data/Endo.rds")
Endo@active.assay <- "SCT"
Endo@active.ident <- as.factor(Endo$type)
DimPlot(Endo)
Endo_4vC_DEG <- FindMarkers(Endo,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Endo_4vC_DEG$change = ifelse(Endo_4vC_DEG$p_val_adj < 0.05 & abs(Endo_4vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Endo_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Endo_4vC_DEG$change)
# Down   Up 
# 35  154

Endo_8vC_DEG <- FindMarkers(Endo,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Endo_8vC_DEG$change = ifelse(Endo_8vC_DEG$p_val_adj < 0.05 & abs(Endo_8vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Endo_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Endo_8vC_DEG$change)
# Down   Up 
# 35  877

## 5.Kupffer
Kupffer <- subset(scc_integrated,idents = "Kupffer")
saveRDS(Kupffer,file="0.Data/Kupffer.rds")
Kupffer@active.assay <- "SCT"
Kupffer@active.ident <- as.factor(Kupffer$type)
DimPlot(Kupffer)
Kupffer_4vC_DEG <- FindMarkers(Kupffer,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Kupffer_4vC_DEG$change = ifelse(Kupffer_4vC_DEG$p_val_adj < 0.05 & abs(Kupffer_4vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(Kupffer_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Kupffer_4vC_DEG$change)
# Down Stable     Up 
# 57      4    170 

Kupffer_8vC_DEG <- FindMarkers(Kupffer,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Kupffer_8vC_DEG$change = ifelse(Kupffer_8vC_DEG$p_val_adj < 0.05 & abs(Kupffer_8vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(Kupffer_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Kupffer_8vC_DEG$change)
# Down Stable     Up 
# 59      3    1243


## 6.LCM
LCM <- subset(scc_integrated,idents = "LCM")
saveRDS(LCM,file="0.Data/LCM.rds")
LCM@active.assay <- "SCT"
LCM@active.ident <- as.factor(LCM$type)
DimPlot(LCM)
LCM_4vC_DEG <- FindMarkers(LCM,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LCM_4vC_DEG$change = ifelse(LCM_4vC_DEG$p_val_adj < 0.05 & abs(LCM_4vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(LCM_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LCM_4vC_DEG$change)
# Down Stable     Up 
# 14      6     67  

LCM_8vC_DEG <- FindMarkers(LCM,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LCM_8vC_DEG$change = ifelse(LCM_8vC_DEG$p_val_adj < 0.05 & abs(LCM_8vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(LCM_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LCM_8vC_DEG$change)
# Down Stable     Up 
# 18     24   1139

## 7.pDCs
pDCs <- subset(scc_integrated,idents = "pDCs")
saveRDS(pDCs,file="0.Data/pDCs.rds")
pDCs@active.assay <- "SCT"
pDCs@active.ident <- as.factor(pDCs$type)
DimPlot(pDCs)
pDCs_4vC_DEG <- FindMarkers(pDCs,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
pDCs_4vC_DEG$change = ifelse(pDCs_4vC_DEG$p_val_adj < 0.05 & abs(pDCs_4vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(pDCs_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(pDCs_4vC_DEG$change)
# Down Stable     Up 
# 12     23     30 

pDCs_8vC_DEG <- FindMarkers(pDCs,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
pDCs_8vC_DEG$change = ifelse(pDCs_8vC_DEG$p_val_adj < 0.05 & abs(pDCs_8vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(pDCs_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(pDCs_8vC_DEG$change)
# Down Stable     Up 
# 7    320    774

## 8.Neutro
Neutro <- subset(scc_integrated,idents = "Neutro")
saveRDS(Neutro,file="0.Data/Neutro.rds")
Neutro@active.assay <- "SCT"
Neutro@active.ident <- as.factor(Neutro$type)
DimPlot(Neutro)
Neutro_4vC_DEG <- FindMarkers(Neutro,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Neutro_4vC_DEG$change = ifelse(Neutro_4vC_DEG$p_val_adj < 0.05 & abs(Neutro_4vC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Neutro_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Neutro_4vC_DEG$change)
# Down   Up 
# 24   99 

Neutro_8vC_DEG <- FindMarkers(Neutro,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Neutro_8vC_DEG$change = ifelse(Neutro_8vC_DEG$p_val_adj < 0.05 & abs(Neutro_8vC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(Neutro_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Neutro_8vC_DEG$change)
# Down Stable     Up 
# 24     13   1025

## 9.B_lymph
B_lymph <- subset(scc_integrated,idents = "B lymph")
saveRDS(B_lymph,file="0.Data/B_lymph.rds")
B_lymph@active.assay <- "SCT"
B_lymph@active.ident <- as.factor(B_lymph$type)
DimPlot(B_lymph)
B_lymph_4vC_DEG <- FindMarkers(B_lymph,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
B_lymph_4vC_DEG$change = ifelse(B_lymph_4vC_DEG$p_val_adj < 0.05 & abs(B_lymph_4vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(B_lymph_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(B_lymph_4vC_DEG$change)
# Down   Up 
# 18   53 

B_lymph_8vC_DEG <- FindMarkers(B_lymph,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
B_lymph_8vC_DEG$change = ifelse(B_lymph_8vC_DEG$p_val_adj < 0.05 & abs(B_lymph_8vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(B_lymph_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(B_lymph_8vC_DEG$change)
# Down   Up 
# 13  746 


## 10.T_lymph
T_lymph <- subset(scc_integrated,idents = "T lymph")
saveRDS(T_lymph,file="0.Data/T_lymph.rds")
T_lymph@active.assay <- "SCT"
T_lymph@active.ident <- as.factor(T_lymph$type)
DimPlot(T_lymph)
T_lymph_4vC_DEG <- FindMarkers(T_lymph,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
T_lymph_4vC_DEG$change = ifelse(T_lymph_4vC_DEG$p_val_adj < 0.05 & abs(T_lymph_4vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(T_lymph_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_lymph_4vC_DEG$change)
# Down Stable     Up 
# 13      2     54

T_lymph_8vC_DEG <- FindMarkers(T_lymph,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                               min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
T_lymph_8vC_DEG$change = ifelse(T_lymph_8vC_DEG$p_val_adj < 0.05 & abs(T_lymph_8vC_DEG$avg_log2FC) >= 0.25, 
                                ifelse(T_lymph_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(T_lymph_8vC_DEG$change)
# Down Stable     Up 
# 30      4   1098  

## 11.NK
NK <- subset(scc_integrated,idents = "NK")
saveRDS(NK,file="0.Data/NK.rds")
NK@active.assay <- "SCT"
NK@active.ident <- as.factor(NK$type)
DimPlot(NK)
NK_4vC_DEG <- FindMarkers(NK,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                          min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
NK_4vC_DEG$change = ifelse(NK_4vC_DEG$p_val_adj < 0.05 & abs(NK_4vC_DEG$avg_log2FC) >= 0.25, 
                           ifelse(NK_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(NK_4vC_DEG$change)
# Down   Up 
# 18   53 

NK_8vC_DEG <- FindMarkers(NK,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                          min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
NK_8vC_DEG$change = ifelse(NK_8vC_DEG$p_val_adj < 0.05 & abs(NK_8vC_DEG$avg_log2FC) >= 0.25, 
                           ifelse(NK_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(NK_8vC_DEG$change)
# Down   Up 
# 13  746 



##############
library(pheatmap)
library(patchwork)
library(viridis)


Cell_type_DEGs <- read.csv("Cell_type_DEGs.csv",row.names = 1)
Cell_type_DEGs
Cell_type_DEGs_up <- Cell_type_DEGs[1:2,]
Cell_type_DEGs_down <- Cell_type_DEGs[3:4,]

Cell_type_DEGs_heatmap_up <- pheatmap(Cell_type_DEGs_up, 
                                      show_rownames = T, show_colnames = T,
                                      cluster_cols = F, cluster_rows = F,
                                      # scale = "column",
                                      cellwidth = 20,cellheight = 20,
                                      # annotation_col = annotation_col,
                                      # annotation_colors = ann_colors,
                                      treeheight_row = 0,
                                      treeheight_col = 0,
                                      border_color = "black",viridis(100))
Cell_type_DEGs_heatmap_down <- pheatmap(Cell_type_DEGs_down, 
                                        show_rownames = T, show_colnames = T,
                                        cluster_cols = F, cluster_rows = F,
                                        # scale = "column",
                                        cellwidth = 20,cellheight = 20,
                                        # annotation_col = annotation_col,
                                        # annotation_colors = ann_colors,
                                        treeheight_row = 0,
                                        treeheight_col = 0,
                                        border_color = "black",viridis(100))

















