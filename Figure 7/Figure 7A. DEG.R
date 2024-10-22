library(tidyverse)
library(data.table)
rm(list = ls())
gc()
## 差异分析
# 加载芯片数据并定义分组
expr <- fread('../TCGA_TPM.csv') %>% aggregate(.~id,.,mean)%>% column_to_rownames(var ='id')
# expr = TPM
group <- fread('../TCGA_risk.csv') %>% as.data.frame()
table(group$risk)
expr = expr[,substr(colnames(expr),14,16) %in% "01A"]
colnames(expr) = substr(colnames(expr),1,12)
expr = expr[,colnames(expr) %in% group$id]
group = group[group$id %in% colnames(expr),]

expr <- expr[,group$id]
boxplot(expr, outline=F, notch=T, col=group,las = 2,main = "TCGA")
expr = log2(1+expr)

group <- factor(group$risk, levels = c("low", "high"))

contrast <- paste0(rev(levels(group)), collapse = "-")
contrast
# 制作分组矩阵
design <- model.matrix(~group)
colnames(design) <- levels(group)
# 制作差异表达矩阵
library(limma)
contrast.matrix <- makeContrasts(contrast, levels = design)
contrast.matrix

fit <- lmFit(expr, design)
fit <- eBayes(fit)
DEG <- topTable(fit,adjust='fdr',coef=2, number = Inf) 

diff_res=subset(DEG,DEG$P.Value<0.05)
write.table(diff_res, file='diffgene_1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#diff_res1 <- subset(res,abs(logFC)>2 & res$P.Value<0.05)
#write.table(diff_res1, file='diffgene_2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
####火山图和热图
#进行注释
all_res=DEG
logFC_cutoff <- 0
risk1 = (all_res$P.Value<0.05)&(all_res$logFC < -logFC_cutoff)
risk2 = (all_res$P.Value<0.05)&(all_res$logFC> logFC_cutoff)
all_res$change = ifelse(risk1,"DOWN",ifelse(risk2,"UP","NOT"))
table(all_res$change)
write.csv(all_res,'TCGA_diff_all.csv')