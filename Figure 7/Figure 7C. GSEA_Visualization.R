library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
rm(list = ls())
###输入文件###
genelist_input <- fread(file="./diff_all.csv") %>% as.data.frame()%>% dplyr::select(1:2) %>% arrange(desc(logFC))
colnames(genelist_input)[1] = "Gene"
genename <- as.character(genelist_input$Gene) #提取第一列基因名

gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<-"Gene"
# write.csv(as.data.frame(gene_map),"gene_transfer.csv",row.names =F)#导出结果至默认路径下

aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)

geneList = genelist_input[,2]
names(geneList) = as.character(genelist_input[,1])
geneList

HALLMARK_gmt <- read.gmt("h.all.v7.5.1.symbols.gmt")
dim(HALLMARK_gmt)
res  <- GSEA(geneList,TERM2GENE = HALLMARK_gmt,nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
hm = as.data.frame(res)

### 第二步：GSEA可视化
# devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
# KEGG_gseresult
# retain curve
gseaNb(object = res,
       geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
       newGsea = T)

# add gene name
# gene = c('PTX3','FOXC2',"SFRP1")

gseaNb(object = res,
       geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
       subPlot = 2,
       termWidth = 5,
       addGene = gene,
       legend.position = c(0.8,0.8),
       newGsea = T,
       addPval = T
       )
