# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggsci)
library(ggplot2) 
options(future.globals.maxSize = 20000 * 1024^2)

# （1）绘制 heatamap plot

pdf("Heatmap.pdf",width = 30,height = 15)

# DoHeatmap(scc_integrated, slot = "data",features = c(
#   "Slc27a2","Lrp2","Slc22a8","Slc5a2","Slc5a12","Fxyd2","Slc17a3", # proximal tubule(PT)
#   "Atp11a", "Slc13a3","Slc34a1","Gpx3"# Proximal convoluted tubule cell(PTC)
#   "Aqp1","Bst1" # Descending loop of Henle (DLH)
#   "Slc12a1","Umod","Cldn8","Krt18","Krt8",  # Ascending loop of Henle(ALH)
#   "Slc12a3",  # Distal convoluted tubule(DCT)
#   "Atp6v0d2","Atp6v1g3","Slc4a1","Aqp6","Slc26a4","Hmx2" # Collecting duct intercalated cell(CD-IC)
#   "Aqp2","Hsd11b2", #  Collecting duct principal / epithelial cell (CD-PC)
#   "Rhbg","Insrr","Stmn1", # Collecting duct transitional cell (CD-TC)
#   "Cdca3","Mki67" # Novel cell
#   "Nrp1","Kdr","Ehd3","Plat","Vim","S100a4","Aqp1","Bst1", # Endothelial(Endo)
#   "Nphs1", "Nphs2", # Podocyte (Podo)
#   "Vim","S100a4" # Pericytes and vascular smooth muscle (Peri)
#   "Plac8", # Fibroblast (Fibro)
#   "C1qa","C1qb", # Macrophage (Macro)
#   "Cd79a", "Cd79b", # B lymphocyte (B lymph)
#   "Cxcr6","Ltb","Il7r","Cd3d","Cd3e","Ifng",  # T lymphocyte (T lymph)
#   "Gzma","Nkg7","Gnly", # Natural killer cell (NK)
#   "Lyz","Cd14", # Monocytes (Mono)
#   "S100a8","S100a9" # Neutrophil (Neutro)
# )) + NoLegend()

dev.off()

# （2）绘制小提琴图

scc_integrated <-readRDS("0.Data/SCT_combine.sample.anno.rds")
DimPlot(scc_integrated)

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

# my36colors <-c(
#   '#E95C59', '#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3',
#   '#E63863','#E4C755','#E59CC4','#AB3282','#23452F','#BD956A',
#   '#8C549C','#58A4C3',"#00BFC4",'#E39A35')

feature = c(
  "Kdr","Bmp2", # Endo
  "Cd79a","Cd79b", # B lymph
  "S100a8","S100a9", # Neutro
  "C1qa","Csf1r",   # Kupffer 
  "Alb","Apoa1", # Hepatocytes
  "Trbc1","Trbc2", # T lymph
  "Il2rb","Nkg7",  # NK
  "S100a4","Ccr2", # Macro
  "Runx2","Ccr9", # pDCs
  "Reln","Rbp1",# HSCs
  "Epcam","Sox9")  # Cholangiocytes

StackedVlnPlot(scc_integrated, c(feature), pt.size=0, cols=my36colors)


