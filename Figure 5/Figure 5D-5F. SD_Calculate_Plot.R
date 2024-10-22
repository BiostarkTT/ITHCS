library(tidyverse)
library(Seurat)
data = data.table::fread('GSE53625_riskscore_13all.csv') %>% column_to_rownames('id')

variances <- apply(data, 2, var)

variance_df <- data.frame(Signature = names(variances), Variance = variances)

# 按照方差排序
variance_df_sorted <- variance_df[order(variance_df$Variance), ]
variance_df_sorted$Variance= log(1+variance_df_sorted$Variance)

col3<-colorRampPalette(c('blue','white','firebrick2'))(14)
variance_df_sorted$Signature = factor(variance_df_sorted$Signature,levels = variance_df_sorted$Signature)
# 使用ggplot2绘制条形图
p1 <- ggplot(variance_df_sorted, aes(x = reorder(Signature, Variance), y = Variance,fill = Signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col3)+
  coord_flip() +  # 翻转坐标轴以便更好地显示长条形标签
  labs(x = "Signature", y = "Variance", title = "GSE53625") +
  theme_minimal()+
  theme(
    panel.grid = element_blank(),  # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA),  # 加上黑色边框
    axis.line = element_line(color = "black"))+
  NoLegend()
p1

data = data.table::fread('TCGA_riskscore_13all.csv') %>% column_to_rownames('id')

variances <- apply(data, 2, var)

variance_df <- data.frame(Signature = names(variances), Variance = variances)

# 按照方差排序
variance_df_sorted <- variance_df[order(variance_df$Variance), ]
variance_df_sorted$Variance= log(1+variance_df_sorted$Variance)

col3<-colorRampPalette(c('blue','white','firebrick2'))(14)
variance_df_sorted$Signature = factor(variance_df_sorted$Signature,levels = variance_df_sorted$Signature)
# 使用ggplot2绘制条形图
p2 <- ggplot(variance_df_sorted, aes(x = reorder(Signature, Variance), y = Variance,fill = Signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col3)+
  coord_flip() +  # 翻转坐标轴以便更好地显示长条形标签
  labs(x = "Signature", y = "Variance", title = "TCGA-ESCC") +
  theme_minimal()+
  theme(
    panel.grid = element_blank(),  # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA),  # 加上黑色边框
    axis.line = element_line(color = "black"))+
  NoLegend()
p2

data = data.table::fread('zhang_riskscore_13all.csv') %>% column_to_rownames('id')

variances <- apply(data, 2, var)

variance_df <- data.frame(Signature = names(variances), Variance = variances)

# 按照方差排序
variance_df_sorted <- variance_df[order(variance_df$Variance), ]
variance_df_sorted$Variance= log(1+variance_df_sorted$Variance)

col3<-colorRampPalette(c('blue','white','firebrick2'))(14)
variance_df_sorted$Signature = factor(variance_df_sorted$Signature,levels = variance_df_sorted$Signature)
# 使用ggplot2绘制条形图
p3 <- ggplot(variance_df_sorted, aes(x = reorder(Signature, Variance), y = Variance,fill = Signature)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col3)+
  coord_flip() +  # 翻转坐标轴以便更好地显示长条形标签
  labs(x = "Signature", y = "Variance", title = "Zhang et al.") +
  theme_minimal()+
  theme(
    panel.grid = element_blank(),  # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA),  # 加上黑色边框
    axis.line = element_line(color = "black"))+
      NoLegend()
p3
pdf('InterITH_Bulk.pdf',width = 8,height = 4)
cowplot::plot_grid(p1,p2,p3,cols = 3)
dev.off()
