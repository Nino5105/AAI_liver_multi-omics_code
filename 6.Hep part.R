library(Seurat)
library(ggsci)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

Hep <- readRDS("0.Data/Hep.rds")
Hep
DimPlot(Hep,split.by = "type")

Hep <- RunPCA(object = Hep,verbose = FALSE)
Hep <- FindNeighbors(Hep,dim=1:15)
Hep <- FindClusters(Hep,resolution = 0.8)
Hep <- RunUMAP (Hep,reduction="pca", dims = 1:15)

Hep_all <- DimPlot(Hep,pt.size = 1,label.size = 5,label = T)

DoHeatmap(Hep,c(
  "Cdh1",#("Ecad"),
  "Alb","Ass1","Asl","Cyp2f2", # PV
  "Hamp","Hamp2","Igfbp2", # Mid
  "Glul","Cyp2e1" # PC 
)) 

# 
Hep2 <-  Hep

new.cluster.ids <- c("Hep1",#0
                     "Hep1",
                     "Hep2",
                     "Hep2",
                     "Hep1",
                     "Hep2",#5
                     "Hep1",
                     "Hep1",
                     "Hep2",
                     "Hep2",
                     "Hep1",#10
                     "Hep3",
                     "Hep1",
                     "Hep1",
                     "Hep2")

names(new.cluster.ids) <- levels(Hep2)
Hep2 <- RenameIdents(Hep2, new.cluster.ids)

levels(Hep2) <- c("Hep1","Hep2","Hep3")

Hep@active.assay <- "SCT"

DimPlot(Hep2,label = TRUE,reduction = "umap",label.size = 5) + scale_color_npg()


modify_vlnplot <- function(obj,
                           feature,
                           pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "mm"),
                           ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

my36colors <-c(
  '#E64933', '#4DBBD5','#009C81')

feature = c(
  "Alb", # PV
  "Ass1", 
  "Hamp", # Mid
  "Cyp2e1",
  "Kdr",
  "Nrp1",
  "Flna",
  "Actg1"
)

StackedVlnPlot(Hep2, c(feature), pt.size=0, cols=my36colors)


# Find DEGS

markers <- FindAllMarkers(object = Hep2, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)


cluster_go<-function(para_cluster) {
  markers_a <- markers[which(markers$cluster == para_cluster),]
  genenames <- as.character(markers_a$gene)
  Go_BP <- enrichGO(gene = genenames, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                    ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  write.csv(Go_BP@result,file=paste0('GO_analysis_in_',para_cluster,'.csv'))
  p <- dotplot(Go_BP, showCategory = 10,title = paste0("The GO enrichment(BP) analysis of cluster",para_cluster))
  return(p)
}

cluster<-unique(Idents(Hep2))

table(markers$cluster)
# p0<-cluster_go(cluster[0])

i = NULL

pdf("Cluster_BP.pdf",width = 8,height = 8)

# cluster_go("Hep2")

for (i in c("Hep1","Hep2", "Hep3")){
  BP_plot <- cluster_go(i)
  print(BP_plot)
}
dev.off()

write.csv(markers,"Hep_DEGs.csv")

saveRDS(Hep2,"Hep_anno.rds")


#######

GO <- read.csv("Hep/GO-subtype.csv")
GO$log10Pvalue <- -log10(GO$p.adjust)
GO
levels(GO$subtype) <- c("Hep1","Hep2","Hep3")

library(ggpubr)

ggbarplot(GO, x="Description", y="log10Pvalue", fill = "subtype", color = "white",
          palette =  c("Hep1" = '#E64933', "Hep2" = '#4DBBD5',"Hep3" = '#009C81'),
          sort.val = "desc", 
          sort.by.grodowns=TRUE, 
          x.text.angle=0, 
          xlab = NULL) + coord_flip() # 5x10


Hep <- readRDS("Hep/Hep_anno.rds")
Hep$cell_subtype <- Hep@active.ident

DimPlot(Hep,label = T,label.size = 5,cols = c("Hep1" = '#E64933',
                                               "Hep2" = '#4DBBD5',
                                               "Hep3" = '#009C81'),
        split.by = "type",reduction = "umap")
 
Hep1 <- subset(Hep,ident = "Hep1")
DimPlot(Hep1)

table(Hep1$type)
# Control AA (4 weeks) AA (8 weeks) 
# 886          980          529

Hep1@active.assay = "SCT"
Hep1@active.ident <- as.factor(Hep$type)
DimPlot(Hep1)

Hep1_4vC_DEG <- FindMarkers(Hep1,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hep1_4vC_DEG$change = ifelse(Hep1_4vC_DEG$p_val_adj < 0.05 & abs(Hep1_4vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hep1_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hep1_4vC_DEG$change)
# Down Stable     Up 
# 111     56    215

Hep1_4vC_DEG_up <- Hep1_4vC_DEG[Hep1_4vC_DEG$change == "Up",]
Hep1_4vC_DEG_down <- Hep1_4vC_DEG[Hep1_4vC_DEG$change == "Down",]

Hep1_8vC_DEG <- FindMarkers(Hep1,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hep1_8vC_DEG$change = ifelse(Hep1_8vC_DEG$p_val_adj < 0.05 & abs(Hep1_8vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hep1_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hep1_8vC_DEG$change)
# Down Stable     Up 
# 110     82   1193 

Hep1_8vC_DEG_up <- Hep1_8vC_DEG[Hep1_8vC_DEG$change == "Up",]
Hep1_8vC_DEG_down <- Hep1_8vC_DEG[Hep1_8vC_DEG$change == "Down",]

Hep1_4vC_uni_up <- setdiff(rownames(Hep1_4vC_DEG_up),rownames(Hep1_8vC_DEG_up))
Hep1_8vC_uni_up <- setdiff(rownames(Hep1_8vC_DEG_up),rownames(Hep1_4vC_DEG_up))
Hep1_4_8_overlap_up <- intersect(rownames(Hep1_4vC_DEG_up),rownames(Hep1_8vC_DEG_up))
# Hep1_4vC_uni %in% rownames(Hep1_4vC_DEG_up)

Hep1_4vC_uni_down <- setdiff(rownames(Hep1_4vC_DEG_down),rownames(Hep1_8vC_DEG_down))
Hep1_8vC_uni_down <- setdiff(rownames(Hep1_8vC_DEG_down),rownames(Hep1_4vC_DEG_down))
Hep1_4_8_overlap_down <- intersect(rownames(Hep1_4vC_DEG_down),rownames(Hep1_8vC_DEG_down))

AverageExpression_value <-  AverageExpression(Hep1, assays = "SCT", 
                                              features = c(Hep1_4vC_uni_up,Hep1_4_8_overlap_up,Hep1_8vC_uni_up,
                                                           Hep1_4vC_uni_down,Hep1_4_8_overlap_down,Hep1_8vC_uni_down), 
                                              return.seurat = FALSE, 
                                              group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,
                   cluster_rows = F,show_rownames = F,
                   border_color = "black",
                   angle_col = 45,
                   # gaps_col = c(4,7,8),
                   gaps_row = c(1262),
                   cellwidth = 20,cellheight = 0.1,
                   cluster_cols = F,
                   scale = "row") #5*5

Hep2 <- subset(Hep,ident = "Hep2")
DimPlot(Hep2)

table(Hep2$type)
# Control AA (4 weeks) AA (8 weeks) 
# 689          725          255

Hep2@active.assay = "SCT"
Hep2@active.ident <- as.factor(Hep2$type)
DimPlot(Hep2)

Hep2_4vC_DEG <- FindMarkers(Hep2,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hep2_4vC_DEG$change = ifelse(Hep2_4vC_DEG$p_val_adj < 0.05 & abs(Hep2_4vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(Hep2_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hep2_4vC_DEG$change)
# Down Stable     Up 
# 46     27    273

Hep2_4vC_DEG_up <- Hep2_4vC_DEG[Hep2_4vC_DEG$change == "Up",]
Hep2_4vC_DEG_down <- Hep2_4vC_DEG[Hep2_4vC_DEG$change == "Down",]

Hep2_8vC_DEG <- FindMarkers(Hep2,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hep2_8vC_DEG$change = ifelse(Hep2_8vC_DEG$p_val_adj < 0.05 & abs(Hep2_8vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(Hep2_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hep2_8vC_DEG$change)
# Down Stable     Up 
# 106    227    745

Hep2_8vC_DEG_up <- Hep2_8vC_DEG[Hep2_8vC_DEG$change == "Up",]
Hep2_8vC_DEG_down <- Hep2_8vC_DEG[Hep2_8vC_DEG$change == "Down",]

Hep2_4vC_uni_up <- setdiff(rownames(Hep2_4vC_DEG_up),rownames(Hep2_8vC_DEG_up))
Hep2_8vC_uni_up <- setdiff(rownames(Hep2_8vC_DEG_up),rownames(Hep2_4vC_DEG_up))
Hep2_4_8_overlap_up <- intersect(rownames(Hep2_4vC_DEG_up),rownames(Hep2_8vC_DEG_up))
# Hep2_4vC_uni %in% rownames(Hep2_4vC_DEG_up)

Hep2_4vC_uni_down <- setdiff(rownames(Hep2_4vC_DEG_down),rownames(Hep2_8vC_DEG_down))
Hep2_8vC_uni_down <- setdiff(rownames(Hep2_8vC_DEG_down),rownames(Hep2_4vC_DEG_down))
Hep2_4_8_overlap_down <- intersect(rownames(Hep2_4vC_DEG_down),rownames(Hep2_8vC_DEG_down))

AverageExpression_value <-  AverageExpression(Hep2, assays = "SCT", 
                                              features = c(Hep2_4vC_uni_up,Hep2_4_8_overlap_up,Hep2_8vC_uni_up,
                                                           Hep2_4vC_uni_down,Hep2_4_8_overlap_down,Hep2_8vC_uni_down), 
                                              return.seurat = FALSE, 
                                              group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,
                   cluster_rows = F,show_rownames = F,
                   border_color = "black",
                   angle_col = 45,
                   # gaps_col = c(4,7,8),
                   gaps_row = c(833),
                   cellwidth = 20,cellheight = 0.1,
                   cluster_cols = F,
                   scale = "row") #5*5

 Hep <- readRDS("3.1. Hep/ Hep_anno.rds")
DimPlot( Hep)

 Hep$cell_subtype <-  Hep@active.ident

dfsam <- as.data.frame(table( Hep$type, Hep$cell_subtype))

control <- dfsam[dfsam$Var1 == "Control",]
myLabel = as.vector(control$Var2)
myLabel = paste(myLabel, "(", round(control$Freq / sum(control$Freq) * 100, 2), "%)"
                , sep = "")  

p1 <- ggplot(control, aes(x = "", y = Freq,fill = Var2),colour = cols ) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "Control")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +     
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_npg(labels = myLabel) +
  labs(title =  paste0("Control"),fill = "Subtype") 

p1


treatment_E <- dfsam[dfsam$Var1 == "AA (4 weeks)",]
myLabel = as.vector(treatment_E$Var2)
myLabel = paste(myLabel, "(", round(treatment_E$Freq / sum(treatment_E$Freq) * 100, 2), "%)"
                , sep = "") 
p2 <- ggplot(treatment_E, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "AA (4 weeks)")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +     
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_npg(labels = myLabel) +
  labs(title =  paste0("AA (4 weeks)"),fill = "Subtype") 
p2

treatment_L <- dfsam[dfsam$Var1 == "AA (8 weeks)",]
myLabel = as.vector(treatment_L$Var2)
myLabel = paste(myLabel, "(", round(treatment_L$Freq / sum(treatment_L$Freq) * 100, 2), "%)"
                , sep = "") 
p3 <- ggplot(treatment_L, aes(x = "", y = Freq,fill = Var2)) + 
  geom_bar(stat = "identity")+
  coord_polar(theta = "y") + 
  theme_bw() + 
  labs(x = "", y = "", title = "AA (8 weeks)")+
  theme(axis.ticks = element_blank())+
  theme(axis.text.x = element_blank()) + 
  theme(panel.grid=element_blank()) +     
  theme(panel.border=element_blank()) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + 
  scale_fill_npg(labels = myLabel) +
  labs(title =  paste0("AA (8 weeks)"),fill = "Subtype") 
p3

p1+p2+p3 #5*15
