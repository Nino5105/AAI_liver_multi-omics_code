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

Lymph <- readRDS("0.Data/Lymph.rds")
Lymph
DimPlot(Lymph)
DimPlot(Lymph,split.by = "type")

Lymph <- RunPCA(object = Lymph,verbose = FALSE)
Lymph <- FindNeighbors(Lymph,dim=1:15)
Lymph <- FindClusters(Lymph,resolution = 0.8)
Lymph <- RunUMAP (Lymph,reduction="pca", dims = 1:15)

Lymph_all <- DimPlot(Lymph,pt.size = 1,label.size = 5,label = T)
Lymph_all

Lymph2 <-  Lymph
Lymph2$id <- Lymph2@active.ident

Lymph_all <- DimPlot(Lymph,pt.size = 1,label.size = 5,
                 label = T)
# Lymph_class <- DimPlot(Lymph,pt.size = 1,label.size = 5,
#                    label = T,split.by = "class")
# Lymph_all 
# Lymph_class

# 分布位置
Lymph@active.assay <- "SCT"
DoHeatmap(Lymph,c(
  "Cd3d", # T lymph
  "Cd3e", 
  "Cd3g",
  "Cd4",
  "Cd8a", 
  "Cd8b1",
  "Ncr1",# NK
  "Nkg7", 
  "Cd79a", # B lymph
  "Cd79b"
)) 


# B subtype
DoHeatmap(subset(Lymph,idents = c(0,1,2,6,11,19)),c(
  "Ighd",# B lymph naive
  "Fcmr",
  "Ighm",
  "Tcl1a",
  "Fam129c",
  "Igha",
  
  # "Ebi3",# B lymph memory
  # "Fcrl4",
  # "Dusp4",
  # "Gpr183",
  # "Ccr6",
  # "Clecl1",
  
  "Prdm1",# B lymph plasma
  "Xbp1",
  "Ccr10",
  "Txndc5",
  "Ssr4",
  "Lman1"
)) 


# 
Lymph2 <-  Lymph

new.cluster.ids <- c("B naïve",
                     "B naïve",
                     "B naïve",
                     "T Memory",
                     "CD8+ naïve",
                     "NK cyto",
                     "B naïve",
                     "NK inflam",
                     "CD8+ CTL",
                     "CD4+ naïve",
                     "CD4+ effector",
                     "B naïve",
                     "CD8+ CTL",
                     "CD4+ Treg",
                     "CD8+ naïve",
                     "CD8+ CTL",
                     "CD8+ naïve",
                     "T Memory",
                     "NK inflam",
                     "B plasma")
table(new.cluster.ids)
names(new.cluster.ids) <- levels(Lymph2)
Lymph2 <- RenameIdents(Lymph2, new.cluster.ids)

levels(Lymph2) <- c("CD8+ naïve","CD8+ CTL","CD4+ naïve","CD4+ effector","CD4+ Treg",
                    "T Memory","B naïve","B plasma","NK cyto","NK inflam")

Lymph2@active.assay <- "SCT"

DoHeatmap(Lymph2,c(
  "Cd3d",  "Cd4",  "Cd8a",   "Cd8b1",
  
  "Foxp3",# T lymph regulatory
  "Ctla4",
  
  "Lef1",# T lymph naive
  "Sell",
  "Ccr7",
  "Tcf7",
  
  "Cxcr3",# T lymph memory
  "Cxcr6",
  "Cd40lg",
  "S100a4",
  
  "Fasl",# T lymph effector
  "Ifng",
  "Gzmk",

  "Ncr1",# NK
  "Nkg7",
  "Prf1",
  "Xcl1",
  
  "Ighd",# B lymph naive
  "Fcmr",
  "Ighd",
  "Jchain"))


DimPlot(Lymph2,label = TRUE,reduction = "umap",label.size = 5)

## heatmap

AverageExpression_value <-  AverageExpression(Lymph2, assays = "SCT", 
                                              features = c("Cd3d","Cd4","Cd8a","Cd8b1",
                                                           "Cd79a","Cd79b", "Ncr1","Nkg7",
                                                           "Sell","Ccr7",
                                                           "Fasl","Ifng",
                                                           "Foxp3", "Ctla4",# T lymph regulatory
                                                           "Cxcr3","Cxcr6",
                                                           "Ighd","Fcmr",# B lymph naive
                                                           "Igha","Jchain","Prf1","Xcl1"                                                           ),
                                              return.seurat = FALSE, group.by = "ident")

pheatmap::pheatmap(AverageExpression_value$SCT,cluster_rows = F,border_color = "black",
                   angle_col = 45,
                   gaps_col = c(6,8,10),
                   gaps_row = c(4,6,8,10,12,14,16,20),
                   cellwidth = 15,cellheight = 10,cluster_cols = F,scale = "row")

# sankey plot

data <- reshape2::melt(table(Lymph2@active.ident,Lymph2$type))
head(data)

ggplot(data,aes(x = Var2, stratum = Var1, alluvium = Var1,
                y = value, fill = Var1, label = Var1)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  labs(fill="Cell subtype") + ylab("Cell number") +
  theme_bw() +
  theme(legend.position = c(0.99,0.99),legend.justification = c(1,1)) +
  theme(axis.ticks = element_blank())

# CTL Subset
CTL <- subset(Lymph2,idents = "CD8+ CTL")
CTL@active.ident <- CTL$type

DimPlot(CTL,pt.size = 1,ncol = 1) + scale_color_d3() + NoLegend()  
DimPlot(CTL,split.by = "type",pt.size = 1,ncol = 1) + scale_color_d3() 
#### Find DEGS

table(CTL$type)
# Control AA (4 weeks) AA (8 weeks) 
# 439          705         1265

CTL@active.assay = "SCT"
CTL@active.ident <- as.factor(CTL$type)
DimPlot(CTL)

CTL_4vC_DEG <- FindMarkers(CTL,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
CTL_4vC_DEG$change = ifelse(CTL_4vC_DEG$p_val_adj < 0.05 & abs(CTL_4vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(CTL_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTL_4vC_DEG$change)
# Down Stable     Up 
# 18     37     98

CTL_4vC_DEG_up <- CTL_4vC_DEG[CTL_4vC_DEG$change == "Up",]
CTL_4vC_DEG_down <- CTL_4vC_DEG[CTL_4vC_DEG$change == "Down",]

CTL_8vC_DEG <- FindMarkers(CTL,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                           min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
CTL_8vC_DEG$change = ifelse(CTL_8vC_DEG$p_val_adj < 0.05 & abs(CTL_8vC_DEG$avg_log2FC) >= 0.25, 
                            ifelse(CTL_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTL_8vC_DEG$change)
# Down Stable     Up 
# 52     26    908 

CTL_4vC_DEG$type <- "4w_vs_Con"
CTL_8vC_DEG$type <- "8w_vs_Con"

ALL_DEG <- rbind(CTL_4vC_DEG,CTL_8vC_DEG)
head(ALL_DEG)
ALL_DEG$type  <- factor(ALL_DEG$type,levels = rev(c("4w_vs_Con","8w_vs_Con")))

ALL_DEG_sub <- ALL_DEG[ALL_DEG$change != "Stable",]
table(ALL_DEG_sub$change)
head(ALL_DEG_sub)
ggplot() + geom_point() +
  geom_jitter(ALL_DEG_sub, mapping=aes(x=type, y= avg_log2FC, color=change),width=0.4,size=1) +
  theme_test()    +
  theme(axis.text.x=element_text(size=10,angle=0,face ='bold'),
        axis.text.y=element_text(size=10,angle=90,face ='bold'),
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 14),
        legend.position = c(0.95,0.95),legend.justification = c(1,1)) + 
  xlab(NULL) + ylab('avg_log2FC')  + 
  geom_hline(yintercept = 0,lty=1,lwd=1.2,alpha=0.5)+
  scale_colour_manual(values = c("#6B9AC7","#E24D36")) +
  guides(color=guide_legend(override.aes = list(size=6))) + coord_flip() # 6X10 +

table(row.names(ALL_DEG_sub)[ALL_DEG_sub$type == "4w_vs_Con"] %in% row.names(ALL_DEG_sub)[ALL_DEG_sub$type == "8w_vs_Con"])
# FALSE 
# 116


#######

rownames(CTL_4vC_DEG_up)

rownames(CTL_8vC_DEG_up)

CTL_4vC_DEG_up_Go <- enrichGO(gene = rownames(CTL_4vC_DEG_up), 
                              OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = 'BH',
                           pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CTL_8vC_DEG_up_Go <- enrichGO(gene = rownames(CTL_8vC_DEG_up), 
                              OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL",
                              ont = "BP", pAdjustMethod = 'BH',
                              pvalueCutoff = 0.05,qvalueCutoff = 0.05)

CTL_4vC_DEG_up_Go@result[c(1:20),c("Description")]
dotplot(CTL_4vC_DEG_up_Go)
CTL_4vC_DEG_up_Go_5 <- CTL_4vC_DEG_up_Go@result[c(1,2,9,10,16),c("Description","pvalue")]
CTL_4vC_DEG_up_Go_5$log10Pvalue <- -log10(CTL_4vC_DEG_up_Go_5$pvalue)
CTL_4vC_DEG_up_Go_5$subtype <- "4w_vs_Con"

CTL_8vC_DEG_up_Go@result[c(1:30),c("Description")]
dotplot(CTL_8vC_DEG_up_Go)
CTL_8vC_DEG_up_Go_5 <- CTL_8vC_DEG_up_Go@result[c(13,16,17,20,29),c("Description","pvalue")]
CTL_8vC_DEG_up_Go_5$log10Pvalue <- -log10(CTL_8vC_DEG_up_Go_5$pvalue)
CTL_8vC_DEG_up_Go_5$subtype <- "8w_vs_Con"

CTL_DEG_up_Go <- rbind(CTL_8vC_DEG_up_Go_5,CTL_4vC_DEG_up_Go_5)
CTL_DEG_up_Go$item <- rownames(CTL_DEG_up_Go)
colnames(CTL_DEG_up_Go)[3] <- "-log10Pvalue"

CTL_DEG_up_Go$subtype <- factor(CTL_DEG_up_Go$subtype,levels = c("8w_vs_Con","4w_vs_Con"))

library(ggpubr)

ggbarplot(CTL_DEG_up_Go, x="Description", y="-log10Pvalue", 
          fill = "subtype", color = "white",width = 0.8,
          # palette =  c("Lymph1" = '#E64933', "Lymph2" = '#4DBBD5',"Lymph3" = '#009C81'),
          sort.val = "desc", 
          sort.by.grodowns=TRUE, 
          x.text.angle=0, 
          xlab = NULL) + coord_flip() # 5x10

library(viridis)

Hepo <- readRDS("Hepo/Hepo_anno.rds")
Hepo$cell_subtype <- Hepo@active.ident

DimPlot(Hepo,label = T,label.size = 5,cols = c("Hepo1" = '#E64933',
                                               "Hepo2" = '#4DBBD5',
                                               "Hepo3" = '#009C81'),
        split.by = "type",reduction = "umap")
 
Hepo1 <- subset(Hepo,ident = "Hepo1")
DimPlot(Hepo1)

table(Hepo1$type)
# Control AA (4 weeks) AA (8 weeks) 
# 886          980          529

Hepo1@active.assay = "SCT"
Hepo1@active.ident <- as.factor(Hepo$type)
DimPlot(Hepo1)

Hepo1_4vC_DEG <- FindMarkers(Hepo1,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo1_4vC_DEG$change = ifelse(Hepo1_4vC_DEG$p_val_adj < 0.05 & abs(Hepo1_4vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hepo1_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo1_4vC_DEG$change)
# Down Stable     Up 
# 111     56    215

Hepo1_4vC_DEG_up <- Hepo1_4vC_DEG[Hepo1_4vC_DEG$change == "Up",]
Hepo1_4vC_DEG_down <- Hepo1_4vC_DEG[Hepo1_4vC_DEG$change == "Down",]

Hepo1_8vC_DEG <- FindMarkers(Hepo1,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                            min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo1_8vC_DEG$change = ifelse(Hepo1_8vC_DEG$p_val_adj < 0.05 & abs(Hepo1_8vC_DEG$avg_log2FC) >= 0.25, 
                             ifelse(Hepo1_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo1_8vC_DEG$change)
# Down Stable     Up 
# 110     82   1193 

Hepo1_8vC_DEG_up <- Hepo1_8vC_DEG[Hepo1_8vC_DEG$change == "Up",]
Hepo1_8vC_DEG_down <- Hepo1_8vC_DEG[Hepo1_8vC_DEG$change == "Down",]

Hepo1_4vC_uni_up <- setdiff(rownames(Hepo1_4vC_DEG_up),rownames(Hepo1_8vC_DEG_up))
Hepo1_8vC_uni_up <- setdiff(rownames(Hepo1_8vC_DEG_up),rownames(Hepo1_4vC_DEG_up))
Hepo1_4_8_overlap_up <- intersect(rownames(Hepo1_4vC_DEG_up),rownames(Hepo1_8vC_DEG_up))
# Hepo1_4vC_uni %in% rownames(Hepo1_4vC_DEG_up)

Hepo1_4vC_uni_down <- setdiff(rownames(Hepo1_4vC_DEG_down),rownames(Hepo1_8vC_DEG_down))
Hepo1_8vC_uni_down <- setdiff(rownames(Hepo1_8vC_DEG_down),rownames(Hepo1_4vC_DEG_down))
Hepo1_4_8_overlap_down <- intersect(rownames(Hepo1_4vC_DEG_down),rownames(Hepo1_8vC_DEG_down))

AverageExpression_value <-  AverageExpression(Hepo1, assays = "SCT", 
                                              features = c(Hepo1_4vC_uni_up,Hepo1_4_8_overlap_up,Hepo1_8vC_uni_up,
                                                           Hepo1_4vC_uni_down,Hepo1_4_8_overlap_down,Hepo1_8vC_uni_down), 
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

Hepo2 <- subset(Hepo,ident = "Hepo2")
DimPlot(Hepo2)

table(Hepo2$type)
# Control AA (4 weeks) AA (8 weeks) 
# 689          725          255

Hepo2@active.assay = "SCT"
Hepo2@active.ident <- as.factor(Hepo2$type)
DimPlot(Hepo2)

Hepo2_4vC_DEG <- FindMarkers(Hepo2,ident.1 = "AA (4 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo2_4vC_DEG$change = ifelse(Hepo2_4vC_DEG$p_val_adj < 0.05 & abs(Hepo2_4vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(Hepo2_4vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo2_4vC_DEG$change)
# Down Stable     Up 
# 46     27    273

Hepo2_4vC_DEG_up <- Hepo2_4vC_DEG[Hepo2_4vC_DEG$change == "Up",]
Hepo2_4vC_DEG_down <- Hepo2_4vC_DEG[Hepo2_4vC_DEG$change == "Down",]

Hepo2_8vC_DEG <- FindMarkers(Hepo2,ident.1 = "AA (8 weeks)",ident.2 = "Control",
                             min.pct = 0.1, logfc.threshold = 0.25) # ident.1 vs ident.2
Hepo2_8vC_DEG$change = ifelse(Hepo2_8vC_DEG$p_val_adj < 0.05 & abs(Hepo2_8vC_DEG$avg_log2FC) >= 0.25, 
                              ifelse(Hepo2_8vC_DEG$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(Hepo2_8vC_DEG$change)
# Down Stable     Up 
# 106    227    745

Hepo2_8vC_DEG_up <- Hepo2_8vC_DEG[Hepo2_8vC_DEG$change == "Up",]
Hepo2_8vC_DEG_down <- Hepo2_8vC_DEG[Hepo2_8vC_DEG$change == "Down",]

Hepo2_4vC_uni_up <- setdiff(rownames(Hepo2_4vC_DEG_up),rownames(Hepo2_8vC_DEG_up))
Hepo2_8vC_uni_up <- setdiff(rownames(Hepo2_8vC_DEG_up),rownames(Hepo2_4vC_DEG_up))
Hepo2_4_8_overlap_up <- intersect(rownames(Hepo2_4vC_DEG_up),rownames(Hepo2_8vC_DEG_up))
# Hepo2_4vC_uni %in% rownames(Hepo2_4vC_DEG_up)

Hepo2_4vC_uni_down <- setdiff(rownames(Hepo2_4vC_DEG_down),rownames(Hepo2_8vC_DEG_down))
Hepo2_8vC_uni_down <- setdiff(rownames(Hepo2_8vC_DEG_down),rownames(Hepo2_4vC_DEG_down))
Hepo2_4_8_overlap_down <- intersect(rownames(Hepo2_4vC_DEG_down),rownames(Hepo2_8vC_DEG_down))

AverageExpression_value <-  AverageExpression(Hepo2, assays = "SCT", 
                                              features = c(Hepo2_4vC_uni_up,Hepo2_4_8_overlap_up,Hepo2_8vC_uni_up,
                                                           Hepo2_4vC_uni_down,Hepo2_4_8_overlap_down,Hepo2_8vC_uni_down), 
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

data <- readRDS("3.3 Lymph/Lympo_anno.rds")
data <- subset(data,ident = "CD8+ CTL")
table(data@active.ident)

DefaultAssay(data) <- "RNA"
DimPlot(data)


expr_matrix <- as(as.matrix(data@assays$RNA@counts), 'sparseMatrix')


colnames(data@meta.data)
p_data <- data@meta.data[,c(4,5,13)]


f_data <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))


pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()) # 表达矩阵是counts值

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#Removing 177 outliers

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 100))
length(expressed_genes) # m = 8805

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~class",cores=10)



# 差异表达基因作为轨迹构建的基因,阈值为qval<0.001,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.1) #选出2431个基因
table(diff$qval <= 0.1)

deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

# 差异基因的结果文件保存
# write.csv(deg,file="../8.Cell trajectory/monocle_data_DEG.csv",quote=F)

# ①轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)   #将“轨迹构建基因嵌入cds对象
table(cds@featureData@data[["use_for_ordering"]]) # 查看“轨迹构建基因”
plot_ordering_genes(cds) # 轨迹构建基因

# ②对数据进行降维
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')

# ③根据order gene的表达趋势，将细胞排序并完成轨迹构建

cds <- orderCells(cds)
# cds2 <- orderCells(cds,root_state = 2) ##⚠️使用root_state参数可以设置拟时间轴的根

cds$State[cds$State == 2] = "3" 
cds$State[cds$State == 4] = "3" 
cds$State[cds$State == 5] = "2" 

p0 <- plot_cell_trajectory(cds,color_by="Pseudotime",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=F) 
p1 <- plot_cell_trajectory(cds,color_by="State",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=F)
# scale_color_manual(breaks = c("1", "2", "3"), values=c("#979797","#EF4F5B","#7990C8"))
p2 <- plot_cell_trajectory(cds,color_by="type",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=F) + scale_color_d3()
# p3 <- plot_cell_trajectory(cds,color_by="cell_type_sub",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=F) + 
#   scale_color_manual(breaks = c("data M1", "data M2", "data Pro"), values=c("#1F77B4","#FF7F0E","#2CA02C"))
# # p4 <- plot_cell_trajectory(cds2,color_by="sample",show_tree = T,cell_size = 1,cell_link_size = 2,show_backbone=TRUE)+ scale_color_nejm()

p0 + p1 + p2

 
# plot_cell_trajectory(cds, color_by = "cell_type_sub", cell_size = 1,cell_link_size = 2,show_backbone=TRUE) + facet_wrap("~type", nrow = 1)+ scale_color_lancet()

# 检测拟时的差异基因
# 
# Time_diff <- differentialGeneTest(cds[ordergene,], cores = 5,ullModelFormulaStr = "~sm.ns(Pseudotime)")
# head(Time_diff)
# Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]
# 
# Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
# # head(Time_genes)
# 
# p <- plot_pseudotime_heatmap(cds[Time_diff,], num_clusters = 4, show_rownames=T, return_heatmap=T)
# p$tree_row

# 通过BEAM进行统计分析

plot_cell_trajectory(cds, color_by = "State")

BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 10)
head(BEAM_res)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

table(BEAM_res$qval < 0.001)
BEAM_genes <- row.names(subset(BEAM_res,qval < 0.001)) # 选择qval< 0.001,n =3361

plot <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                    branch_point = 1,num_clusters = 5,
                                    cores = 10, 
                                    use_gene_short_name = T,
                                    show_rownames = F,
                                    return_heatmap = T) #

# ggsave(plot$ph_res, "BEAM_heatmap.pdf", width = 4, height = 5)



plot <- plot_genes_branched_heatmap(cds[BEAM_genes,],
                                    branch_point = 1,num_clusters = 5,
                                    cores = 10, 
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T) #

plot$ph_res


gene_group <- plot$annotation_row
gene_group$gene <- rownames(gene_group)
allcluster_go=data.frame()

C1_gene <- gene_group[gene_group$Cluster == 1,]$gene #n=945
C1_BP <- enrichGO(gene = C1_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C1_BP@result[c(1:20),c("Description","pvalue")]

# Description       pvalue
# GO:0034341                            response to interferon-gamma 9.413039e-17
# GO:0022407                        regulation of cell-cell adhesion 7.135159e-16
# GO:0007159                            leukocyte cell-cell adhesion 2.056747e-15
# GO:1903037              regulation of leukocyte cell-cell adhesion 1.155498e-14
# GO:0071346                   cellular response to interferon-gamma 1.649048e-14
# GO:2000377 regulation of reactive oxygen species metabolic process 2.232410e-14
# GO:0001819              positive regulation of cytokine production 2.452758e-14
# GO:0045862                      positive regulation of proteolysis 3.423095e-14
# GO:1901652                                     response to peptide 4.500434e-14
# GO:0050900                                     leukocyte migration 4.661031e-14


C2_gene <- gene_group[gene_group$Cluster == 2,]$gene #n=791
C2_BP <- enrichGO(gene = C2_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C2_BP@result[c(1:10),c("Description","pvalue")]
# Description       pvalue
# GO:0006914                                      autophagy 1.044263e-22
# GO:0061919         process utilizing autophagic mechanism 1.044263e-22
# GO:0016236                                 dataautophagy 2.918399e-19
# GO:0007033                           vacuole organization 2.472665e-17
# GO:0002221 pattern recognition receptor signaling pathway 6.199727e-15
# GO:0010506                        regulation of autophagy 1.311984e-14
# GO:0002274                   data leukocyte activation 1.618353e-13
# GO:1905037                     autophagosome organization 3.956518e-12
# GO:0006898                  receptor-mediated endocytosis 1.262085e-11
# GO:0002224           toll-like receptor signaling pathway 3.485216e-11

C3_gene <- gene_group[gene_group$Cluster == 3,]$gene #n=464
C3_BP <- enrichGO(gene = C3_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C3_BP@result[c(1:10),c("Description","pvalue")]
# Description       pvalue
# GO:0044282             small molecule catabolic process 8.269818e-27
# GO:0016054               organic acid catabolic process 9.865204e-24
# GO:0046395            carboxylic acid catabolic process 9.865204e-24
# GO:0015711                      organic anion transport 1.239687e-22
# GO:1901605           alpha-amino acid metabolic process 1.876707e-21
# GO:0006520        cellular amino acid metabolic process 5.082649e-20
# GO:0006631                 fatty acid metabolic process 8.547026e-20
# GO:1901615   organic hydroxy compound metabolic process 6.738344e-19
# GO:0072521 purine-containing compound metabolic process 1.324269e-18
# GO:0006790            sulfur compound metabolic process 1.699003e-18

C4_gene <- gene_group[gene_group$Cluster == 4,]$gene #n=635
C4_BP <- enrichGO(gene = C4_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C4_BP@result[c(1:10),c("Description","pvalue")]
# Description       pvalue
# GO:0030099                         data cell differentiation 3.514756e-14
# GO:0002683         negative regulation of immune system process 3.410562e-13
# GO:1903706                            regulation of hemopoiesis 3.527257e-13
# GO:0030098                           lymphocyte differentiation 4.792726e-13
# GO:0032103 positive regulation of response to external stimulus 4.107262e-12
# GO:1902105              regulation of leukocyte differentiation 5.108378e-12
# GO:0009615                                    response to virus 1.552339e-11
# GO:0051607                            defense response to virus 2.392963e-11
# GO:0022407                     regulation of cell-cell adhesion 3.194645e-11
# GO:0001819           positive regulation of cytokine production 3.632018e-11

C5_gene <- gene_group[gene_group$Cluster == 5,]$gene #n=526
C5_BP <- enrichGO(gene = C5_gene,OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = "BP",pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C5_BP@result[c(1:10),c("Description","pvalue")]
# 
# Description       pvalue
# GO:0006119                              oxidative phosphorylation 1.525688e-40
# GO:0006091         generation of precursor metabolites and energy 3.333301e-40
# GO:0046034                                  ATP metabolic process 4.779036e-37
# GO:0045333                                   cellular respiration 5.636872e-33
# GO:0042775 mitochondrial ATP synthesis coupled electron transport 5.929860e-32
# GO:0022900                               electron transport chain 6.333431e-32
# GO:0022904                   respiratory electron transport chain 4.892449e-31
# GO:0042773               ATP synthesis coupled electron transport 6.992298e-31
# GO:0015980    energy derivation by oxidation of organic compounds 1.188526e-26
# GO:0010257                    NADH dehydrogenase complex assembly 2.076714e-26



head(allcluster_go[,c("ID","Description","qvalue","cluster")])