library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer) 
library(viridis)
library(wesanderson)
gene <- data.table::fread('/work/lt/LT_Project/1.ESCC_ITH/Data/Q4gene.txt')
gene=gene$id
Q4_score=list(gene)
names(Q4_score)='Q4_score'

DefaultAssay(scRNA_SCT) <- 'RNA'
scRNA_SCT <- AddModuleScore(object = scRNA_SCT, features = Q4_score, name ='Q4_score')

FeaturePlot(scRNA_SCT,features = 'Q4_score1',reduction = 'umap',label=T,pt.size = 0)

dat <- scRNA_SCT@meta.data %>% select(celltype,Q4_score1)

pdf('/work/lt/LT_Project/1.ESCC_ITH/Result/GSE160269_Q4Sig_box.pdf',width = 4,height = 3.5)
ggplot(dat, aes(x = celltype, y = Q4_score1, fill = celltype)) +
  geom_boxplot(color = "black", alpha = 0.5, width = 0.5) +  # 设置边框颜色为黑色
  scale_fill_manual(values = c("ImmuneCells" = "#1f78b4", "EpithelialCells" = "#33a02c", "StromalCells" = "#e31a1c")) +  # 设置填充颜色
  theme_minimal() +
  theme(panel.grid.major = element_blank(),   # 去掉背景中的主网格线
        panel.grid.minor = element_blank(),   # 去掉背景中的次网格线
        panel.background = element_rect(color = "black", fill = NA)) +  # 添加边框
  labs(title = "GSE160269",
       x = "",
       y = "Q4 signature score")+
  geom_signif(comparisons = list(c("ImmuneCells", "EpithelialCells"),
                                 c("ImmuneCells", "StromalCells"),
                                 c("EpithelialCells", "StromalCells")),
              map_signif_level = TRUE,y_position = c(7.5, 8.5, 9.5))+NoLegend()
dev.off()