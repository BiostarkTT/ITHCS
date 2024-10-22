library(ggplot2)

df <- readRDS("./scRNA_GSE196756_exp.rds")

# 绘制箱式图
ggplot(df, aes(x = Group, y = Q4_Signature_Score, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Q4 Signature Score") +
  theme(legend.position = "none")
