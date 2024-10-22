library(plyr)
library(ggpubr)
library(tidyverse)     
library(data.table)
rt=fread('./riskscore_boxplot.csv')   #读取输入文件
x=colnames(rt)[1]
y=colnames(rt)[3]
colnames(rt)=c("Type","Sample","sd_value")

#定义输出图片的排序方式
med=ddply(rt,"Type",summarise,med=median(sd_value))
rt$Type=factor(rt$Type, levels=med[order(med[,"med"],decreasing = F),"Type"])

#绘制
library(RColorBrewer)
library(ggsci)
Set2<-brewer.pal(7,"Set2")
Set1<-brewer.pal(7,"Set1")
colors <- c('#CA6855','#546B7E',Set2,Set1)
col=rainbow(length(levels(factor(rt$Type))))
p=ggboxplot(rt, x="Type", y="sd_value", color = "Type",
       palette = colors,
       ylab=y,
       xlab=x,
       #add = "jitter",                                            #绘制每个样品的散点
       legend = "right")
p
pdf(file="sd_value_Sig.pdf", width=10, height=5) #输出图片文件
print(p+rotate_x_text(45))#倾斜角度
dev.off()

