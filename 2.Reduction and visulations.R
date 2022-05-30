# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)

# （1）数据降维

scc_integrated <-readRDS("../3.Seurat_workflow/0.Data/SCT_combine.sample.pca.rds")
scc_integrated <- RunPCA(object = scc_integrated,verbose = FALSE)
pdf("PCA_Elbowplot.pdf",height=10,width=10)
ElbowPlot(scc_integrated,ndims = 50)
dev.off()

# （2）选择PCA个数

# pdf("cluster.pdf",,height=10,width=10)
# for (i in 20:40){
#   scc_integrated <- FindNeighbors(scc_integrated,dim=1:i)
#   scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
#   scc_integrated <- RunUMAP (scc_integrated,reduction="pca", dims = 1:i)
#   p1 <- DimPlot(scc_integrated,label = TRUE,reduction = "umap") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 20)) + labs(title =  paste0("PCA n = ",i))
#   scc_integrated <- RunTSNE(scc_integrated,dims = 1:i)
#   p2 <- TSNEPlot (scc_integrated,label = TRUE) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 20)) + labs(title =  paste0("PCA n = ",i))
#   print(p1/p2)
#   # }
#   # dev.off()
#   saveRDS(scc_integrated,file="SCT_combine.sample.pca.rds")
  
  scc_integrated <- FindNeighbors(scc_integrated,dim=1:40)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
  scc_integrated <- RunUMAP (scc_integrated,reduction="pca", dims = 1:40)
  scc_integrated <- RunTSNE(scc_integrated,dims = 1:40)
  scc_integrated$id <- scc_integrated@active.ident
  
  # （3）选择resolution
  
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.1)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.2)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.4)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.6)
  scc_integrated <- FindClusters(scc_integrated,resolution = 0.8)
  scc_integrated <- FindClusters(scc_integrated,resolution = 1)
  clustree(scc_integrated)
  
  scc_integrated$integrated_snn_res.0.8 #选择0.8为resolution 值
  
  # （4）绘制聚类图
  
  # DimPlot(scc_integrated,label = T,pt.size = 1,label.size = 5,reduction = "tsne")
  DimPlot(scc_integrated,label = T,label.size = 5,reduction = "umap")
  
  label.size = 5
  
  p1 <- DimPlot(scc_integrated,label = T,group.by = "id") + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cluster id"))
  p2 <- DimPlot(scc_integrated,label = F,group.by = "sample")+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Sample id")) + scale_color_hue()
  p3 <- DimPlot(scc_integrated,label = F,group.by = "type",cols = c('Control' = '#00BFC4', 'Treatment-E' = 'orange','Treatment-L' = "#F8766D")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
  p3 <- DimPlot(scc_integrated,label = F,group.by = "type",split.by = "type",cols = c('Control' = '#00BFC4', 'Treatment-E' = 'orange','Treatment-L' = "#F8766D")) + theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Dataset type"))
  (p1|p2|p3)
  
  #（5）统计各组以及各样本在不同亚簇的比例构成
  
  scc_integrated$id <- scc_integrated@active.ident
  
  dfsam <- as.data.frame(table(scc_integrated$sample,scc_integrated$id))
  dfsam$Var1 <- factor(dfsam$Var1,levels = rev(c("LA1","LA2","LA3","LA4","LA5","LA6","LC1","LC2","LC3")))
  dfsam$Var2 <- factor(dfsam$Var2,levels = rev(levels(dfsam$Var2)))
  
  ggplot(data = dfsam,aes(x = Var2,y=Freq,fill = Var1)) +
    geom_bar(stat="identity",position = "fill",width = 0.5)+
    scale_y_continuous(labels = scales::percent)+
    theme_classic() + 
    labs(y = 'Fraction of different sample',x="") +
    coord_flip() + 
    # scale_fill_nejm() + 
    NoLegend()
  
  dfsam2 <- as.data.frame(table(scc_integrated$type,scc_integrated$id))
  dfsam2$Var2 <- factor(dfsam2$Var2,levels = rev(levels(dfsam2$Var2)))
  
  ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
    geom_bar(stat="identity",position = "fill",width = 0.5)+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = rev(c("#F8766D","orange","#00BFC4")), name = "" )+
    theme_classic() + 
    labs(y = 'Fraction of different group',x="") + 
    coord_flip()  +
    NoLegend()
  
  saveRDS(scc_integrated,file="../3.Seurat_workflow/0.Data/SCT_combine.sample.umap.rds")
  
  
  table(Idents(scc_integrated))
  