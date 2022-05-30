library(Seurat)
library(ggsci)
library(ggplot2) 
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)

# loading the Endo dataset
Endo <- readRDS("0.Data/Endo.rds")
DimPlot(Endo)
head(Endo@meta.data)
FeaturePlot(Endo,c("Cdh5","Alb"),cols = c("grey","red")) # VE-cadherin(all ECs)


# re-clustering the Endo dataset
Endo <- RunPCA(object = Endo,verbose = FALSE)
# ElbowPlot(Endo,ndims = 50) # 选择主成分为15
Endo <- FindNeighbors(Endo,dim=1:20)
Endo <- FindClusters(Endo,resolution = 0.8)
Endo <- RunUMAP(Endo,reduction="pca", dims = 1:15)
DimPlot(Endo,label = T,label.size = 5,pt.size = 1)

# annotating cellular type of the Endo dataset

DoHeatmap(Endo,c( "Vwf","Fcgr2b","Gpr182","Clec4g", #EC
                  "Rspo3","Wnt9b",# Central
                  "Dll4","Efnb2","Gja5", # Portal
                  "Fabp4","Kit","Thbd", # Pericentral
                  "Ctsl","Hamp","Igfbp2", # Middle
                  "Ltbp4","Ntn4", # Periportal
                  "Lyve1","Flt4","Pdpn","Prox1"))

sub <- subset(Endo,idents = c(3,8,12,13))
DoHeatmap(sub,c(
  "Lyve1","Flt4","Pdpn","Prox1"
))

# head(Endo@active.ident)
Endo@active.ident <- Endo$integrated_snn_res.0.8
new.cluster.ids <- c("LSEC_Mid",
                     "LSEC_Mid",
                     "LSEC_PP",
                     "LSEC_Mid2",
                     "LSEC_Mid",
                     "LSEC_Mid",
                     "LSEC_Mid",
                     "LSEC_PC",
                     "LSEC_Mid2",
                     "LSEC_PP",
                     "LVEC_Cent",
                     "LVEC_Port",
                     "LSEC_Mid2",
                     "LYEC",
                     "LSEC_Mid")
table(new.cluster.ids)

names(new.cluster.ids) <- levels(Endo)
Endo <- RenameIdents(Endo, new.cluster.ids)
levels(Endo) <- c("LVEC_Port",
                  "LSEC_PP",
                  "LSEC_Mid",
                  "LSEC_Mid2",
                  "LSEC_PC",
                  "LVEC_Cent",
                  "LYEC")


AverageExpression_value <-  AverageExpression(Endo, assays = "SCT", 
                                              features = c("Cdh5","Vwf","Fcgr2b","Gpr182","Clec4g", #EC
                                                           "Dll4","Efnb2","Gja5", # Portal
                                                           "Ltbp4","Ntn4", # Periportal
                                                            "Ctsl","Hamp","Igfbp2", # Middle
                                                          "Fabp4","Kit","Thbd", # Pericentral
                                                          "Rspo3","Wnt9b",# Central
                                                          "Lyve1","Flt4"),# lymph
                                              return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 45,
                   gaps_col = c(1,5,6,7),
                   gaps_row = c(4,6,9,11,14,17),
                   cellwidth = 20,cellheight = 10,cluster_cols = F,scale = "row")


table(Endo@active.ident,Endo$type)
            Control AA (4 weeks) AA (8 weeks)
LVEC_Port     402          228          234
LSEC_PP      2492         2033         1034
LSEC_Mid    10392         7716         4408
LSEC_Mid2    3132         1787         1731
LSEC_PC      1084          670          419
LVEC_Cent     443          269          185
LYEC          263          226           42

DimPlot(Endo,label = T,label.size = 5,pt.size = 1) + scale_color_locuszoom()

Endo$cell_subtype <- Endo@active.ident
table(Endo$cell_subtype)
# LVEC_Port   LSEC_PP  LSEC_Mid LSEC_Mid2   LSEC_PC LVEC_Cent      LYEC 
#    864      5559     22516      6650      2173       897       531
saveRDS(Endo,"3.2.Endo/Endo_anno.rds")

# 4.stack plot of cell_subtype proportion 

library(plyr)

dfsam <- as.data.frame(table(Endo$type,Endo$cell_subtype))
dfsam$Var1 <- as.numeric(dfsam$Var1)
dfsam_prop <- ddply(dfsam,"Var1",transform,Percent = Freq / sum(Freq)*100)
ggplot(dfsam_prop,aes(x=Var1,y=Percent,fill=Var2,position = "fill")) +
  geom_area(colour = "black",size =0.5)+theme_classic() + 
  scale_fill_locuszoom() + 
  theme(axis.title.x = element_blank())   #5*4
write.csv(dfsam_prop,"dfsam_prop.csv")


table(Endo$cell_subtype)
# LVEC_Port   LSEC_PP  LSEC_Mid LSEC_Mid2   LSEC_PC LVEC_Cent      LYEC 
    # 864      5559     22516      6650      2173       897       531

LVEC_Port <- subset(Endo,idents = "LVEC_Port")
LSEC_PP <- subset(Endo,idents = "LSEC_PP")
LSEC_Mid <- subset(Endo,idents = "LSEC_Mid")
LSEC_Mid2 <- subset(Endo,idents = "LSEC_Mid2")
LSEC_PC <- subset(Endo,idents = "LSEC_PC")
LVEC_Cent <- subset(Endo,idents = "LVEC_Cent")
LYEC <- subset(Endo,idents = "LYEC")

# Idents(Endo) <- as.factor(Endo$type)
# Idents(LVEC_Port) <- as.factor(LVEC_Port$type)
# Idents(LSEC_PP) <- as.factor(LSEC_PP$type)
Idents(LSEC_Mid) <- as.factor(LSEC_Mid$type)
Idents(LSEC_Mid2) <- as.factor(LSEC_Mid2$type)
# Idents(LSEC_PC) <- as.factor(LSEC_PC$type)
# Idents(LVEC_Cent) <- as.factor(LVEC_Cent$type)
# Idents(LYEC) <- as.factor(LYEC$type)

VlnPlot(LSEC_Mid,c("Col1a1","Col1a2","Tgfb2","Tgfb3","Tgfbr3","Tgfbi"),
        pt.size = 0,slot = "scale.data") + scale_fill_d3()+ 
  theme(axis.title.x = element_blank())
 

AverageExpression_value_Mid <-  AverageExpression(LSEC_Mid, assays = "SCT", 
                                                   features = c(
                                                     "Col1a1","Col1a2","Tgfb2","Tgfb3","Tgfbr3","Tgfbi",
                                                     "Twist1","Zeb1","Zeb2","Vim", # EMT markers
                                                     "Acta2","Sm22","Fn1","S100a4", # Mesenchymal markers
                                                     "Klf2","Fos","Junb",
                                                     # "Klf4","Fosb","Junb", # Vascular Tone regulation
                                                     "Mrc1","Stab1","Stab2",#"Scarb1","Scarb2",
                                                     "Lamp2"  # Endocytic receptors
                                                   ), 
                                                   return.seurat = FALSE, group.by = "ident")


pheatmap::pheatmap(AverageExpression_value_Mid$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 90,
                   # gaps_col = c(1,5,6,7),
                   gaps_row = c(6,10,13,16),
                   cellwidth = 15,cellheight = 12,cluster_cols = F,scale = "none",viridis(100))

AverageExpression_value_Mid2 <-  AverageExpression(LSEC_Mid2, assays = "SCT", 
                                                  features = c(
                                                    "Col1a1","Col1a2","Tgfb2","Tgfb3","Tgfbr3","Tgfbi",
                                                    "Twist1","Zeb1","Zeb2","Vim", # EMT markers
                                                    "Acta2","Sm22","Fn1","S100a4", # Mesenchymal markers
                                                    "Klf2","Fos","Junb",
                                                    # "Klf4","Fosb","Junb", # Vascular Tone regulation
                                                    "Mrc1","Stab1","Stab2",#"Scarb1","Scarb2",
                                                               "Lamp2"  # Endocytic receprots
                                                               ), 
                                                  return.seurat = FALSE, group.by = "ident")


pheatmap::pheatmap(AverageExpression_value_Mid2$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 90,
                   # gaps_col = c(1,5,6,7),
                   gaps_row = c(6,10,13,16),
                   cellwidth = 15,cellheight = 12,cluster_cols = F,scale = "none",viridis(100))

VlnPlot(LSEC_Mid,c("Junb"),
        pt.size = 0,slot = "scale.data") + scale_fill_d3()+ 
  theme(axis.title.x = element_blank())


FeaturePlot(LSEC_Mid,c("Junb"),split.by = "type")


LSEC_Mid@active.assay = "SCT"
DimPlot(LSEC_Mid)

LSEC_Mid_4vC_DEG <- FindMarkers(LSEC_Mid,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LSEC_Mid_4vC_DEG$change = ifelse(LSEC_Mid_4vC_DEG$p_val_adj < 0.05 & abs(LSEC_Mid_4vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(LSEC_Mid_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LSEC_Mid_4vC_DEG$change)
# Down   Up 
# 41  168

LSEC_Mid_4vC_DEG_up <- LSEC_Mid_4vC_DEG[LSEC_Mid_4vC_DEG$change == "Up",]
LSEC_Mid_4vC_DEG_down <- LSEC_Mid_4vC_DEG[LSEC_Mid_4vC_DEG$change == "Down",]

LSEC_Mid_8vC_DEG <- FindMarkers(LSEC_Mid,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
LSEC_Mid_8vC_DEG$change = ifelse(LSEC_Mid_8vC_DEG$p_val_adj < 0.05 & abs(LSEC_Mid_8vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(LSEC_Mid_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(LSEC_Mid_8vC_DEG$change)
# Down   Up 
# 32  957
LSEC_Mid_8vC_DEG_up <- LSEC_Mid_8vC_DEG[LSEC_Mid_8vC_DEG$change == "Up",]
LSEC_Mid_8vC_DEG_down <- LSEC_Mid_8vC_DEG[LSEC_Mid_8vC_DEG$change == "Down",]


LSEC_Mid_4vC_DEG_up_GO <- enrichGO(gene = row.names(LSEC_Mid_4vC_DEG_up), OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                  ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
 
LSEC_Mid_8vC_DEG_up_GO <- enrichGO(gene = row.names(LSEC_Mid_8vC_DEG_up), OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

LSEC_Mid_4vC_DEG_down_GO <- enrichGO(gene = row.names(LSEC_Mid_4vC_DEG_down), OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

LSEC_Mid_8vC_DEG_down_GO <- enrichGO(gene = row.names(LSEC_Mid_8vC_DEG_down), OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                                   ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)


dotplot(LSEC_Mid_4vC_DEG_up_GO, showCategory = 10)
dotplot(LSEC_Mid_8vC_DEG_up_GO, showCategory = 10)
dotplot(LSEC_Mid_4vC_DEG_down_GO, showCategory = 10)
dotplot(LSEC_Mid_8vC_DEG_down_GO, showCategory = 10)

LSEC_Mid_4vC_DEG_up_GO_5 <- LSEC_Mid_4vC_DEG_up_GO@result[c(1,2,5,6,7),c("Description","pvalue","Count")]
LSEC_Mid_4vC_DEG_up_GO_5$log10Pvalue <- -log10(LSEC_Mid_4vC_DEG_up_GO_5$pvalue)
LSEC_Mid_4vC_DEG_up_GO_5$subtype <- "LSEC_Mid_4vC_DEG_up"
LSEC_Mid_4vC_DEG_up_GO_5

LSEC_Mid_8vC_DEG_up_GO_5 <- LSEC_Mid_8vC_DEG_up_GO@result[c(1,2,3,5,9),c("Description","pvalue","Count")]
LSEC_Mid_8vC_DEG_up_GO_5$log10Pvalue <- -log10(LSEC_Mid_8vC_DEG_up_GO_5$pvalue)
LSEC_Mid_8vC_DEG_up_GO_5$subtype <- "LSEC_Mid_8vC_DEG_up"
LSEC_Mid_8vC_DEG_up_GO_5

Go_LSEC_Mid_DEG_up <- rbind(LSEC_Mid_4vC_DEG_up_GO_5,LSEC_Mid_8vC_DEG_up_GO_5)
Go_LSEC_Mid_DEG_up$item <- row.names(Go_LSEC_Mid_DEG_up)
Go_LSEC_Mid_DEG_up$subtype <- factor(Go_LSEC_Mid_DEG_up$subtype,levels = c("LSEC_Mid_4vC_DEG_up","LSEC_Mid_8vC_DEG_up"))


p1 <- ggplot(Go_LSEC_Mid_DEG_up,aes(x=subtype,y=Description)) +
  geom_point(aes(size=Count,color=log10Pvalue)) + 
  scale_color_gradient(low="black",high = "red") + 
  theme_bw() + ylab("")+ xlab("")+
  theme(axis.text.y=element_text(size=12,colour="black"), 
        axis.title.y=element_text(size = 12,face="bold",colour="black"), 
        axis.text.x=element_text(size=12,face="bold",colour="black",vjust = 0.5, hjust = 0.4,angle = 45), 
        axis.title.x=element_text(size = 12,face="bold",colour="black"),
        legend.title=element_text(size=12,face="bold",colour="black"),
        # legend.text=element_text(size=12,face="bold",colour="black"), 
        title = element_text(size = 13,face="bold",colour="black"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p1



dotplot(LSEC_Mid_4vC_DEG_down_GO, showCategory = 10)
dotplot(LSEC_Mid_8vC_DEG_up_GO, showCategory = 10)

LSEC_Mid_4vC_DEG_down_GO_5 <- LSEC_Mid_4vC_DEG_down_GO@result[c(1,3,4,5,7),c("Description","pvalue","Count")]
LSEC_Mid_4vC_DEG_down_GO_5$log10Pvalue <- -log10(LSEC_Mid_4vC_DEG_down_GO_5$pvalue)
LSEC_Mid_4vC_DEG_down_GO_5$subtype <- "LSEC_Mid_4vC_DEG_down"
LSEC_Mid_4vC_DEG_down_GO_5

LSEC_Mid_8vC_DEG_down_GO_5 <- LSEC_Mid_8vC_DEG_down_GO@result[c(1,2,5,7,8),c("Description","pvalue","Count")]
LSEC_Mid_8vC_DEG_down_GO_5$log10Pvalue <- -log10(LSEC_Mid_8vC_DEG_down_GO_5$pvalue)
LSEC_Mid_8vC_DEG_down_GO_5$subtype <- "LSEC_Mid_8vC_DEG_down"
LSEC_Mid_8vC_DEG_down_GO_5

Go_LSEC_Mid_DEG_down <- rbind(LSEC_Mid_4vC_DEG_down_GO_5,LSEC_Mid_8vC_DEG_down_GO_5)
Go_LSEC_Mid_DEG_down$item <- row.names(Go_LSEC_Mid_DEG_down)
Go_LSEC_Mid_DEG_down$subtype <- factor(Go_LSEC_Mid_DEG_down$subtype,levels = c("LSEC_Mid_4vC_DEG_down","LSEC_Mid_8vC_DEG_down"))

p2 <- ggplot(Go_LSEC_Mid_DEG_down,aes(x=subtype,y=Description)) +
  geom_point(aes(size=Count,color=log10Pvalue)) + 
  scale_color_gradient(low="black",high = "#00B1FF") + 
  theme_bw() + ylab("")+ xlab("")+
  theme(axis.text.y=element_text(size=12,colour="black"), 
        axis.title.y=element_text(size = 12,face="bold",colour="black"), 
        axis.text.x=element_text(size=12,face="bold",colour="black",vjust = 0.5, hjust = 0.4,angle = 45), 
        axis.title.x=element_text(size = 12,face="bold",colour="black"),
        legend.title=element_text(size=12,face="bold",colour="black"),
        # legend.text=element_text(size=12,face="bold",colour="black"), 
        title = element_text(size = 13,face="bold",colour="black"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
p2


p1

