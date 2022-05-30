# (0) 加载R包
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)
library(ggplot2) 
library(ggpubr)
options(future.globals.maxSize = 20000 * 1024^2)

# （1）RenameIdents


scc_integrated <- readRDS("SCT_combine.sample.umap.rds")

scc_integrated@active.ident <- scc_integrated$id
scc_integrated$id <- scc_integrated@active.ident

new.cluster.ids <- c("B lymph","Endo","Neutro","Endo",
                     "Kupffer","Endo","Endo","Endo",
                     "Hepo","Kupffer","T lymph","T lymph",
                     "LCM","Endo","T lymph","NK",
                     "Hepo","LCM","pDCs","Hepo",
                     "NK","LCM","LCM","Hepo",
                     "Endo","Neutro","Endo",
                     "T lymph","B lymph","B lymph",
                     "Endo","Cho",
                     "Kupffer","T lymph","Neutro","HSC",
                     "Neutro","LCM","Neutro","Neutro","Hepo","B lymph")

names(new.cluster.ids) <- levels(scc_integrated)
scc_integrated <- RenameIdents(scc_integrated, new.cluster.ids)

table(scc_integrated@active.ident)
# B lymph     Endo   Neutro Kupffer      Hepo  T lymph      LCM       NK     pDCs      Cho      HSC 
#   12873    35179     9734     9018     8197     9082     6822     3006     1314      457      273 

levels(scc_integrated) <- c("Hepo","Cho","HSC","Endo","Kupffer","LCM","pDCs","Neutro","B lymph","T lymph","NK")

saveRDS(scc_integrated,file="SCT_combine.anno.rds")

# #（2）Dimplot_Anno_all

pdf("Dimplot_Anno_all.pdf",width = 10,height = 10)
DimPlot(scc_integrated,label = T,label.size = 5,cols = c("Hepo" = '#E95C59',
                                                         "Cho" = '#53A85F',
                                                         "HSC" = '#F1BB72',
                                                         "Endo" = '#F3B1A0',
                                                         "Kupffer" = '#D6E7A3',
                                                         "LCM" = '#57C3F3',
                                                         "pDCs" = '#E63863',
                                                         "Neutro" = '#E4C755',
                                                         "B lymph" = '#E59CC4',
                                                         "T lymph" = '#AB3282',
                                                         "NK" = '#00BFC4')) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell annotation"))
dev.off()

# # （3）Dimplot_Anno_each_type
pdf("Dimplot_Anno_all_split.pdf",width = 20,height = 10)
DimPlot(scc_integrated,label = T,label.size = 5,split.by = "type",cols = c("Hepo" = '#E95C59',
                                                         "Cho" = '#53A85F',
                                                         "HSC" = '#F1BB72',
                                                         "Endo" = '#F3B1A0',
                                                         "Kupffer" = '#D6E7A3',
                                                         "LCM" = '#57C3F3',
                                                         "pDCs" = '#E63863',
                                                         "Neutro" = '#E4C755',
                                                         "B lymph" = '#E59CC4',
                                                         "T lymph" = '#AB3282',
                                                         "NK" = '#00BFC4')) + 
  theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15)) + labs(title =  paste0("Cell annotation"))
dev.off()

# （5）各类细胞数目与占比

table(Idents(scc_integrated),scc_integrated$type)

#           Control Treatment-E Treatment-L
# Hepo       3438        2745        2014
# Cho         166         138         153
# HSC         103          68         102
# Endo      16364       11935        6880
# Kupffer    4308        1762        2948
# LCM        1962        2793        2067
# pDCs        754         430         130
# Neutro     1486        5155        3093
# B lymph    6352        4958        1563
# T lymph    2611        3168        3303
# NK         1095        1392         519

table(Idents(scc_integrated),scc_integrated$orig.ident)

#          LA1  LA2  LA3  LA4  LA5  LA6  LC1  LC2  LC3
# Hepo    1308  783  654  536  687  791  700 1295 1443
# Cho       49   48   41   51   55   47   38   62   66
# HSC       24   29   15   24   49   29    9   65   29
# Endo    4287 3933 3715 2187 2657 2036 3482 5448 7434
# Kupffer  649  671  442  793  861 1294 1166 2182  960
# LCM      994 1153  646  764  440  863  755  742  465
# pDCs     139  172  119   50   20   60  217  297  240
# Neutro  1740 2464  951 1066  914 1113  600  679  207
# B lymph 1101 2369 1488  541  520  502 2568 1994 1790
# T lymph  875 1388  905 1144  739 1420  791 1222  598
# NK       457  543  392  181  122  216  380  452  263

data <- as.data.frame(prop.table(table(Idents(scc_integrated),scc_integrated$orig.ident),margin = 2))
# LA1         LA2         LA3         LA4         LA5         LA6         LC1         LC2
# Hepo    11.25354900  5.77731867  6.98121264  7.30543819  9.72536806  9.44928921  6.53838969  8.96938634
# Cho      0.42157791  0.35416513  0.43766012  0.69510699  0.77859570  0.56146219  0.35494115  0.42942236
# HSC      0.20648714  0.21397477  0.16011956  0.32710917  0.69365798  0.34643412  0.08406501  0.45020086
# Endo    36.88376495 29.01940530 39.65627669 29.80782336 37.61325028 24.32206427 32.52381842 37.73375814
# Kupffer  5.58375635  4.95093337  4.71818958 10.80823225 12.18856172 15.45812926 10.89108911 15.11289652
# LCM      8.55200895  8.50734155  6.89581554 10.41297533  6.22876557 10.30940151  7.05212031  5.13921596
# pDCs     1.19590467  1.26909171  1.27028181  0.68147744  0.28312571  0.71676024  2.02690080  2.05707162
# Neutro  14.97031747 18.18047665 10.15157985 14.52909909 12.93884485 13.29590252  5.60433402  4.70286743
# B lymph  9.47259744 17.47952483 15.88385995  7.37358593  7.36126840  5.99689404 23.98654960 13.81077712
# T lymph  7.52817689 10.24127499  9.66054654 15.59220390 10.46149490 16.96332577  7.38838035  8.46377615
# NK       3.93185924  4.00649303  4.18445773  2.46694834  1.72706682  2.58033688  3.54941154  3.13062751
# 
# LC3
# Hepo    10.69284920
# Cho      0.48907003
# HSC      0.21489441
# Endo    55.08706928
# Kupffer  7.11374583
# LCM      3.44572064
# pDCs     1.77843646
# Neutro   1.53390144
# B lymph 13.26417192
# T lymph  4.43127084
# NK       1.94886995

data$type <- c(rep("Treatment-E",33),rep("Treatment-L",33),rep("Control",33))
head(data)
data <- data[c(67:99,1:66),]

my_comparisons <- list(c("Treatment-E", "Control"), c("Treatment-L", "Control"), c("Treatment-E", "Treatment-L"))
levels(data$type) <- c("Control","Treatment-E","Treatment-L")

pdf("p_anno.pdf",width =16,height = 8)
p <- ggviolin(data, x="type", y="Freq",fill = "type",
               palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              add = "boxplot", add.params = list(fill="white")) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(0.6)) + 
  facet_wrap(~Var1,nrow = 2,dir = "v")
p
dev.off()

#（6）统计各组以及各样本在不同亚簇的比例构成

dfsam <- as.data.frame(table(scc_integrated$sample,scc_integrated$cell_type))
dfsam$Var1 <- factor(dfsam$Var1,levels = rev(c("LA1","LA2","LA3","LA4","LA5","LA6","LC1","LC2","LC3")))
dfsam$Var2 <- factor(dfsam$Var2,levels = levels(dfsam$Var2))

pdf("Fraction_celltype_sample.pdf",width = 8,height = 5)

ggplot(data = dfsam,aes(x = Var1,y=Freq,fill = Var2)) +
  geom_bar(stat="identity",position = "fill",width = 0.8)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic() + 
  labs(y = 'Fraction of different cell',x="") 
   #+ coord_flip()  + scale_fill_nejm() +  NoLegend()

dev.off()


dfsam2 <- as.data.frame(table(scc_integrated$type,scc_integrated$cell_type))
dfsam2$Var2 <- factor(dfsam2$Var2,levels = rev(levels(dfsam2$Var2)))

pdf("Fraction_celltype_type.pdf",width = 5,height = 8)

ggplot(data = dfsam2,aes(x = Var2,y=Freq,fill = Var1)) +
  geom_bar(stat="identity",position = "fill",width = 0.5,)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = rev(c("#00AFBB", "#E7B800", "#FC4E07")), name = "" )+
  theme_classic() + 
  labs(y = 'Fraction of different group',x="") +
  coord_flip()  #+   NoLegend()

dev.off()

pdf("Dimplot.pdf",width = 20,height = 10)
DimPlot(scc_integrated,label = T,split.by = "type",label.size = 5)
