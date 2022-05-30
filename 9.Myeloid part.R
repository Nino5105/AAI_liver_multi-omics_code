
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(monocle) # BiocManager::install("monocle")
library(ggpubr)
library(ggsci)
library(cowplot)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(VennDiagram)

Data <- readRDS("0.Data/SCT_combine.anno.rds")
DimPlot(Data)
table(Data@active.ident)

Myeloid <-  subset(Data,idents = c("Neutro","Kupffer","LCM","pDCs"))
Myeloid

DimPlot(Myeloid)
DimPlot(Myeloid,split.by = "type")

Myeloid <- RunPCA(object = Myeloid,verbose = FALSE)
Myeloid <- FindNeighbors(Myeloid,dim=1:15)
Myeloid <- FindClusters(Myeloid,resolution = 0.8)
Myeloid <- RunUMAP (Myeloid,reduction="pca", dims = 1:15)

DimPlot(Myeloid,pt.size = 1,label.size = 5,label = T)

Myeloid2 <-  Myeloid
Myeloid2$id <- Myeloid2@active.ident

Myeloid_all <- DimPlot(Myeloid2,pt.size = 1,label.size = 5,
                 label = T)
Neutro

Myeloid@active.assay <- "SCT"
DoHeatmap(Myeloid,c(
  "S100a8", # NeuTRO
  
  "Siglech", # pDC
  
  "Itgax", # LCM
  "S100a4",
  # "Ccr2",
  
  "C1qa",
  # "Csf1r",
  "Clec4f", # Kuffer
  
  # "Cd80",
  "Cd86",  # M1
  "Cd68",  
  
  "Tgfb2", # M2
  "Wnt5a",
  "Mrc1",
  
  "Mki67", # Pro
  "Cdca3"
  # "StmN1", 
  # "Cenpf",
  # "Cnnb2"
)) 

subset(Myeloid,idents = c(17,19,14),invert =T)

# 
Myeloid2 <-  Myeloid
# Idents(Myeloid2) <- Myeloid2$cell_type
# DimPlot(Myeloid2)

new.cluster.ids <- c("Neutro","Neutro","Neutro","LCM_M1","Kuffer_M1","Kuffer_M1","Kuffer_M2","Kuffer_M1","Kuffer_M1","LCM_M1","pDC","LCM_M1","Neutro","Kuffer_M1","LCM_M1","LCM_M1","LCM_M1","Kuffer_Pro","LCM_M1","Neutro","Neutro","Neutro","pDC","LCM_Pro","Kuffer_M1")
table(new.cluster.ids)
names(new.cluster.ids) <- levels(Myeloid2)
Myeloid2 <- RenameIdents(Myeloid2, new.cluster.ids)

levels(Myeloid2) <- c("Kuffer_M1","Kuffer_M2","Kuffer_Pro","LCM_M1","LCM_Pro","Neutro","pDC")

Myeloid2@active.assay <- "SCT"

DoHeatmap(Myeloid2,c(
  "S100a8", # Neutro
  
  "Siglech", # pDC
  
  "Itgax", # LCM
  "S100a4",
  # "Ccr2",
  
  "C1qa",
  # "Csf1r",
  "Clec4f", # Kuffer
  
  # "Cd80",
  "Cd86",  # M1
  "Cd68",  
  
  "Tgfb2", # M2
  "Wnt5a",
  "Mrc1",
  
  "Mki67", # Pro
  "Cdca3"))

DimPlot(Myeloid2,label = TRUE,reduction = "umap",label.size = 5) + scale_color_tron()
DimPlot(Myeloid2,label = TRUE,reduction = "umap",label.size = 5,split.by = "type") + scale_color_tron()

FeaturePlot(Myeloid2,c("S100a8","Siglech",
                       "C1qa","Itgax", 
                       "Cd86","Cd68",
                       "Mrc1","Mki67"),cols = c("white","red"),ncol = 4)

## heatmap

AverageExpression_value <-  AverageExpression(Myeloid2, assays = "SCT", 
                                              features = c(
                                                "C1qa",
                                                # "Csf1r",
                                                "Clec4f", # Kuffer
                                                
                                                "Itgax", # LCM
                                                "S100a4",
                                                "Ccr2",
                                                
                                                "Cd86",  # M1
                                                "Cd68",  
                                                
                                                "Tgfb2", # M2
                                                # "Cd163",
                                                "Mrc1",
                                                
                                                "Mki67", # Pro
                                                "Cdca3"),
                                              return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 45,
                   gaps_col = c(3),
                   gaps_row = c(2,5,7,9),
                   cellwidth = 15,cellheight = 10,cluster_cols = F,scale = "row")

# sankey plot

data <- reshape2::melt(table(Myeloid2@active.ident,Myeloid2$type))
head(data)

ggplot(data,aes(x = Var2, stratum = Var1, alluvium = Var1,
                y = value, fill = Var1)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  # geom_text(stat = "stratum", size = 3) +
  labs(fill="Cell subtype") + ylab("Cell number") +
  theme_bw() +
  # theme(legend.position = c(0.99,0.99),legend.justification = c(1,1)) +
  theme(axis.ticks = element_blank()) + scale_fill_tron()

write.csv(table(Myeloid2@active.ident,Myeloid2$type),"3.4 Macro/proportion.csv")

# save RDS

Kuffer_M1 <- subset(Myeloid2,ident="Kuffer_M1")
LCM_M1 <- subset(Myeloid2,ident="LCM_M1")
Neutro <- subset(Myeloid2,ident="Neutro")

saveRDS(Myeloid2,"3.4 Macro/Myeloid_anno.rds")
saveRDS(Kuffer_M1,"3.4 Macro/Kuffer_M1.rds")
saveRDS(LCM_M1,"3.4 Macro/LCM_M1.rds")
saveRDS(Neutro,"3.4 Macro/Neutro.rds")


# Kuffer part
Kuffer_M1 <- readRDS("3.4 Macro/Kuffer_M1.rds")
Idents(Kuffer_M1) <- Kuffer_M1$type
DimPlot(Kuffer_M1)

Kuffer_M1_4vC_DEG <- FindMarkers(Kuffer_M1,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                                min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Kuffer_M1_4vC_DEG$change = ifelse(Kuffer_M1_4vC_DEG$p_val_adj < 0.05 & abs(Kuffer_M1_4vC_DEG$avg_log2FC) >= 0.25, 
                                 ifelse(Kuffer_M1_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Kuffer_M1_4vC_DEG$change)
# Down Stable     Up 
# 53      5    192

Kuffer_M1_4vC_DEG$type <- "Kuffer_M1_4vC"

# Kuffer_M1_4vC_DEG_up <- Kuffer_M1_4vC_DEG[Kuffer_M1_4vC_DEG$change == "Up",]
# Kuffer_M1_4vC_DEG_down <- Kuffer_M1_4vC_DEG[Kuffer_M1_4vC_DEG$change == "Down",]

Kuffer_M1_8vC_DEG <- FindMarkers(Kuffer_M1,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                                min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Kuffer_M1_8vC_DEG$change = ifelse(Kuffer_M1_8vC_DEG$p_val_adj < 0.05 & abs(Kuffer_M1_8vC_DEG$avg_log2FC) >= 0.25, 
                                 ifelse(Kuffer_M1_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Kuffer_M1_8vC_DEG$change)
# Down Stable     Up 
# 65      2   1283

Kuffer_M1_8vC_DEG$type <- "Kuffer_M1_8vC"


# LCM part
LCM_M1 <- readRDS("3.4 Macro/LCM_M1.rds")
Idents(LCM_M1) <- LCM_M1$type
DimPlot(LCM_M1)

LCM_M1_4vC_DEG <- FindMarkers(LCM_M1,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LCM_M1_4vC_DEG$change = ifelse(LCM_M1_4vC_DEG$p_val_adj < 0.05 & abs(LCM_M1_4vC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(LCM_M1_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LCM_M1_4vC_DEG$change)
# Down Stable     Up 
# 13      7     64

LCM_M1_4vC_DEG$type <- "LCM_M1_4vC"

# LCM_M1_4vC_DEG_up <- LCM_M1_4vC_DEG[LCM_M1_4vC_DEG$change == "Up",]
# LCM_M1_4vC_DEG_down <- LCM_M1_4vC_DEG[LCM_M1_4vC_DEG$change == "Down",]

LCM_M1_8vC_DEG <- FindMarkers(LCM_M1,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                              min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LCM_M1_8vC_DEG$change = ifelse(LCM_M1_8vC_DEG$p_val_adj < 0.05 & abs(LCM_M1_8vC_DEG$avg_log2FC) >= 0.25, 
                               ifelse(LCM_M1_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LCM_M1_8vC_DEG$change)
# Down Stable     Up 
# 20     29   1126

LCM_M1_8vC_DEG$type <- "LCM_M1_8vC"

# Merge

ALL_DEG <- rbind(Kuffer_M1_Mid_4vC_DEG,Kuffer_M1_Mid_8vC_DEG,LCM_M1_4vC_DEG,LCM_M1_8vC_DEG)
ALL_DEG$type  <- factor(ALL_DEG$type,levels = c("LCM_M1_8vC","LCM_M1_4vC",
                                                "Kuffer_M1_8vC","Kuffer_M1_4vC"      ))

ALL_DEG_sub <- ALL_DEG[ALL_DEG$change != "Stable",]
table(ALL_DEG_sub$type)

ggplot() + geom_point() +
  geom_jitter(ALL_DEG_sub, mapping=aes(x=type, y= avg_log2FC, color=change),width=0.4,size=1) +
  theme_test()    +
  theme(axis.text.x=element_text(size=10,angle=0,face ='bold'),
        axis.text.y=element_text(size=10,face ='bold'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 14)) + 
  xlab(NULL) + ylab('avg_log2FC')  + 
  ylim(-1.5,1.5)+
  geom_hline(yintercept = 0,lty=1,lwd=2,alpha=0.5)+
  scale_colour_manual(values = c("#6B9AC7","#E24D36")) +
  guides(color=guide_legend(override.aes = list(size=6))) + coord_flip() # 6X10

LCM_M1_4vC_DEG_up <- rownames(LCM_M1_4vC_DEG)[LCM_M1_4vC_DEG$change == "Up"]
LCM_M1_8vC_DEG_up <- rownames(LCM_M1_8vC_DEG)[LCM_M1_8vC_DEG$change == "Up"]
Kuffer_M1_4vC_DEG_up <- rownames(Kuffer_M1_Mid_4vC_DEG)[Kuffer_M1_Mid_4vC_DEG$change == "Up"]
Kuffer_M1_8vC_DEG_up <- rownames(Kuffer_M1_Mid_8vC_DEG)[Kuffer_M1_Mid_8vC_DEG$change == "Up"]

# venn.plot

venn.plot <- 
  venn.diagram(
    x = list(
      'LCM_M1_4vC \n (n = 64)' = LCM_M1_4vC_DEG_up,
      'LCM_M1_8vC \n (n = 1126)' = LCM_M1_8vC_DEG_up,
      'Kuffer_M1_4vC \n (n = 192)' = Kuffer_M1_4vC_DEG_up,
      'Kuffer_M1_8vC \n (n = 1283)' = Kuffer_M1_8vC_DEG_up),
    filename = NULL,
    col = "black",
    fill = c("#00AFBB", "#E7B800", "#FC4E07","skyblue"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap up DEG in different methods",
    main.cex = 1.2)

pdf(file="Overlap up DEP&DEG in different subgroup.pdf")
grid.draw(venn.plot)
dev.off()

