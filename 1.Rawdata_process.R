# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 40000 * 1024^2)

# （1）读取数据
LA1 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA1/outs/filtered_feature_bc_matrix")
LA2 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA2/outs/filtered_feature_bc_matrix")
LA3 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA3/outs/filtered_feature_bc_matrix")
LA4 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA4/outs/filtered_feature_bc_matrix")
LA5 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA5/outs/filtered_feature_bc_matrix")
LA6 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LA6/outs/filtered_feature_bc_matrix")

LC1 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LC1/outs/filtered_feature_bc_matrix")
LC2 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LC2/outs/filtered_feature_bc_matrix")
LC3 <- Read10X(data.dir = "/data/rawdata/P002_rd/2.Cellranger_results/LC3/outs/filtered_feature_bc_matrix")

LA1 <- CreateSeuratObject(counts = LA1, project = "LA1", min.cells = 3, min.features = 200)
LA2 <- CreateSeuratObject(counts = LA2, project = "LA2", min.cells = 3, min.features = 200)
LA3 <- CreateSeuratObject(counts = LA3, project = "LA3", min.cells = 3, min.features = 200)
LA4 <- CreateSeuratObject(counts = LA4, project = "LA4", min.cells = 3, min.features = 200)
LA5 <- CreateSeuratObject(counts = LA5, project = "LA5", min.cells = 3, min.features = 200)
LA6 <- CreateSeuratObject(counts = LA6, project = "LA6", min.cells = 3, min.features = 200)

LC1 <- CreateSeuratObject(counts = LC1, project = "LC1", min.cells = 3, min.features = 200)
LC2 <- CreateSeuratObject(counts = LC2, project = "LC2", min.cells = 3, min.features = 200)
LC3 <- CreateSeuratObject(counts = LC3, project = "LC3", min.cells = 3, min.features = 200)

LA1 # 18660 features across 12278 samples
LA1$sample <- "LA1"
LA1$type <- "Treatment-E"

LA2 # 18560 features across 14225 samples
LA2$sample <- "LA2"
LA2$type <- "Treatment-E"

LA3 # 18209 features across 9799 samples
LA3$sample <- "LA3"
LA3$type <- "Treatment-E"

LA4 # 18240 features across 8549 samples
LA4$sample <- "LA4"
LA4$type <- "Treatment-L"

LA5 # 18288 features across 8451 samples
LA5$sample <- "LA5"
LA5$type <- "Treatment-L"

LA6 # 18539 features across 9461 samples
LA6$sample <- "LA6"
LA6$type <- "Treatment-L"

LC1 # 17972 features across 11255 samples
LC1$sample <- "LC1"
LC1$type <- "Control"

LC2 # 18734 features across 15497 samples
LC2$sample <- "LC2"
LC2$type <- "Control"

LC3 # 18278 features across 13687 samples
LC3$sample <- "LC3"
LC3$type <- "Control"

# （2）统计线粒体相关基因，
LA1[["percent.mt"]] <- PercentageFeatureSet(LA1, pattern = "mt-") 
LA2[["percent.mt"]] <- PercentageFeatureSet(LA2, pattern = "mt-") 
LA3[["percent.mt"]] <- PercentageFeatureSet(LA3, pattern = "mt-") 
LA4[["percent.mt"]] <- PercentageFeatureSet(LA4, pattern = "mt-") 
LA5[["percent.mt"]] <- PercentageFeatureSet(LA5, pattern = "mt-") 
LA6[["percent.mt"]] <- PercentageFeatureSet(LA6, pattern = "mt-") 
LC1[["percent.mt"]] <- PercentageFeatureSet(LC1, pattern = "mt-") 
LC2[["percent.mt"]] <- PercentageFeatureSet(LC2, pattern = "mt-") 
LC3[["percent.mt"]] <- PercentageFeatureSet(LC3, pattern = "mt-")

# （3）提取子集：对检测基因数,reads数以及线粒体含量并进行过滤
LA1 <- subset(LA1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 25)
LA2 <- subset(LA2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 25)
LA3 <- subset(LA3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 25)
LA4 <- subset(LA4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 25)
LA5 <- subset(LA5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 25)
LA6 <- subset(LA6, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 25)
LC1 <- subset(LC1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 25)
LC2 <- subset(LC2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 25)
LC3 <- subset(LC3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 25)

# （3）绘制小提琴图：检测基因数、reads数以及线粒体含量
VlnPlot(LA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 11623 samples
VlnPlot(LA2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 13553 samples
VlnPlot(LA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9368 samples
VlnPlot(LA4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 11623 samples
VlnPlot(LA2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 13553 samples
VlnPlot(LA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 9368 samples
VlnPlot(LC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 10706 samples
VlnPlot(LC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 14438 samples
VlnPlot(LC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # 13495 samples

# （4） 运行 SCT 对各个样本进行数据标准化
LA1 <- SCTransform(LA1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LA2 <- SCTransform(LA2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LA3 <- SCTransform(LA3, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LA4 <- SCTransform(LA4, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LA5 <- SCTransform(LA5, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LA6 <- SCTransform(LA6, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LC1 <- SCTransform(LC1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LC2 <- SCTransform(LC2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
LC3 <- SCTransform(LC3, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)

#saveRDS(c(LA1,LA2,LA3,LA4,LA5,LA6,LC1,LC2,LC3),file="SCT_each.sample.rds")

# （5） 通过 CCA 对多个样本进行数据整合
object_list = list(LA1,LA2,LA3,LA4,LA5,LA6,LC1,LC2,LC3)
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)

remove(LA1,LA2,LA3,LA4,LA5,LA6,LC1,LC2,LC3)
remove(object_list,scc.list)
scc_integrated # 42576 features across 95955 samples within 3 assays
head(scc_integrated@meta.data) #查看meta信息
scc_integrated@assays$RNA@data #查看归一化后数据
saveRDS(scc_integrated,file="SCT_combine.sample.rds")

pdf(file="Vlnplot_all.pdf",height=5,width=10)
VlnPlot(scc_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()





