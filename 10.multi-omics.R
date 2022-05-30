library(VennDiagram) 
library(ggplot2)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)

# 读取数据

DEG_scRNA_4w <- read.csv("scRNA/DEG between scRNA-Seq(W4 vs Con).csv",check.names = F,row.names = 1)
table(DEG_scRNA_4w$type)
DEG_scRNA_4w_up <- DEG_scRNA_4w[DEG_scRNA_4w$type == "Up(318)",]
DEG_scRNA_4w_down <- DEG_scRNA_4w[DEG_scRNA_4w$type == "Down(131)",]
DEG_scRNA_4w_all <- DEG_scRNA_4w[DEG_scRNA_4w$type != "Stable(13186)",]

DEP_4w <- read.csv("Protein/4 weeks/3.Different expressed protein between AA(4w)_vs_Con in mouse liver.csv",check.names = F,row.names = 1)
table(DEP_4w$change)
DEP_4w_up <- DEP_4w[DEP_4w$change == "Up(196)",]
DEP_4w_down <- DEP_4w[DEP_4w$change == "Down(95)",]
DEP_4w_all <- DEP_4w[DEP_4w$change != "Stable(4526)",]

overlap_4w_up_list <- intersect(row.names(DEG_scRNA_4w_up),row.names(DEP_4w_up))
# "Cyp4a14" "Gstt3"   "S100a11" "Gpnmb"   "Abca6"   "S100a9"

overlap_4w_down_list <- intersect(row.names(DEG_scRNA_4w_down),row.names(DEP_4w_down))
# [1] "Hsd3b5"    "Mup7"      "Selenbp2"  "Serpina1e" "Mup20"     "Nudt7"    
# [7] "Cyp8b1" 

overlap_4w_all_list <- intersect(row.names(DEG_scRNA_4w_all),row.names(DEP_4w_all))

scRNA_Protein_list_data <- merge(data.frame(DEG_scRNA_4w[overlap_4w_all_list,c("logFC")],rep("scRNA",14),overlap_4w_all_list),
                                 data.frame(DEP_4w[overlap_4w_all_list,c("logFC")],rep("Protein",14),overlap_4w_all_list),by = "overlap_4w_all_list")

row.names(scRNA_Protein_list_data) <- scRNA_Protein_list_data$overlap_4w_all_list
colnames(scRNA_Protein_list_data) <- c("DEP&DEGs_list","scRNA_Seq","type1","Mass_spec","type2")
head(scRNA_Protein_list_data)

ggscatter(data=scRNA_Protein_list_data, x = "scRNA_Seq", y = "Mass_spec",size = 1,
                title = "AA_4w vs Con",
                conf.int = T,) +
  geom_smooth(method = "lm",fill = "lightblue",size = 1.5)+
  stat_cor(method = "spearman")+ 
  xlab("scRNA-Seq Log2(Fold change)") + ylab("Mass spec Log2(Fold change)") + 
  geom_hline(yintercept=0, linetype=2, colour="black") +
  geom_vline(xintercept=0, linetype=2, colour="black") +
  xlim(-10,10) +ylim(-30,30) +theme(plot.title = element_text(hjust = 0.5))


Go_BP_4w_up <- enrichGO(gene = overlap_4w_up_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_4w_up@result[1:10,]

Go_BP_4w_down <- enrichGO(gene = overlap_4w_down_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

Go_BP_4w_down@result[1:10,]






###########################

DEG_scRNA_8w <- read.csv("scRNA/DEG between scRNA-Seq(W8 vs Con).csv",check.names = F,row.names = 1)
table(DEG_scRNA_8w$type)
DEG_scRNA_8w_up <- DEG_scRNA_8w[DEG_scRNA_8w$type == "Up(601)",]
DEG_scRNA_8w_down <- DEG_scRNA_8w[DEG_scRNA_8w$type == "Down(417)",]
DEG_scRNA_8w_all <- DEG_scRNA_4w[DEG_scRNA_8w$type != "Stable(12680)",]

DEP_8w <- read.csv("Protein/8 weeks/Different expressed protein between Treatment vs control in AA-Liver.csv",check.names = F,row.names = 1)
table(DEP_8w$change)
DEP_8w_up <- DEP_8w[DEP_8w$change == "Up(196)",]
DEP_8w_down <- DEP_8w[DEP_8w$change == "Down(95)",]
DEP_8w_all <- DEP_4w[DEP_8w$change != "Stable(4526)",]

overlap_8w_all_list <- intersect(row.names(DEG_scRNA_8w_all),row.names(DEP_8w_all))

overlap_8w_all_list2 <- c(overlap_8w_up_list,overlap_8w_down_list)

scRNA_Protein_list_data <- merge(data.frame(DEG_scRNA_8w[overlap_8w_all_list2,c("logFC")],rep("scRNA",19),overlap_8w_all_list2),
                                 data.frame(DEP_8w[overlap_8w_all_list2,c("logFC")],rep("Protein",19),overlap_8w_all_list2),by = "overlap_8w_all_list2")

row.names(scRNA_Protein_list_data) <- scRNA_Protein_list_data$overlap_8w_all_list2
colnames(scRNA_Protein_list_data) <- c("DEP&DEGs_list","scRNA_Seq","type1","Mass_spec","type2")
head(scRNA_Protein_list_data)

ggscatter(data=scRNA_Protein_list_data, x = "scRNA_Seq", y = "Mass_spec",size = 1,
          title = "AA_8w vs Con",
          conf.int = T,) +
  geom_smooth(method = "lm",fill = "lightblue",size = 1.5)+
  stat_cor(method = "spearman")+ 
  xlab("scRNA-Seq Log2(Fold change)") + ylab("Mass spec Log2(Fold change)") + 
  geom_hline(yintercept=0, linetype=2, colour="black") +
  geom_vline(xintercept=0, linetype=2, colour="black") +
  xlim(-10,10) +ylim(-5,5) +theme(plot.title = element_text(hjust = 0.5))





overlap_8w_up_list <- intersect(row.names(DEG_scRNA_8w_up),row.names(DEP_8w_up))
# [1] "Chil3"   "Fmo3"    "Rbm3"    "Gsr"     "Akr1b8"  "Soat2"   "Gsta1"  
# [8] "Sfxn5"   "Slco2b1" "Efhd2"   "Ikbke"   "Cyp27a1"

overlap_8w_down_list <- intersect(row.names(DEG_scRNA_8w_down),row.names(DEP_8w_down))
# [1] "Ggcx"     "Lars2"    "Hsd3b5"   "Ehd3"     "Selenbp2" "Mup17"    "Phldb2"


Go_BP_8w_up <- enrichGO(gene = overlap_8w_up_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Go_BP_8w_up@result[1:10,]

Go_BP_8w_down <- enrichGO(gene = overlap_8w_down_list, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                        ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
Go_BP_8w_down@result[1:10,]



###########################


write.csv(overlap_4w_all_list,"overlap_4w_all_list.csv")
write.csv(overlap_8w_all_list2,"overlap_8w_all_list2.csv")
