library(forcats)
library(ggplot2)
library(tidyverse)
library(data.table)
source(file = 'Published_Signature_Formula.R')
# input: title
# title <- paste0(plot.title, "\nin multiple region LIHC cohort")
dat <- data.table::fread('GSE33426_exp.csv') %>% column_to_rownames('id')
group <- data.table::fread('GSE33426_multiRegion_T.txt')
dat <- dat[,group$id]
colnames(dat) <- group$sample
gene <- fread('Q1_gene_coef.csv')

{
  Sig1 <- gene$id
Sig1_coef <- gene$coef


clonal.sig <- data.frame(Sig1,Sig1_coef)
colnames(clonal.sig) <- c('Gene',"beta")
# input: signature 
prognostic_signature <- clonal.sig

# input: expression matrix
gene_matrix <- t(dat) %>% as.data.frame()
gene_matrix$SampleID <- rownames(gene_matrix)
gene_matrix <- gene_matrix %>% select(SampleID,everything())

cut.off <- "median"

plot.title <- 'Clonal-base signature'

file.name <- 'clonal_risk_sig_losic'



# select signature genes
inter.sect <- intersect(colnames(gene_matrix), prognostic_signature$Gene)

if(length(prognostic_signature$Gene) > length(inter.sect)){
  print('The following genes are missed:')
  
  print(setdiff(prognostic_signature$Gene, inter.sect))
}

prognostic_signature <- subset(prognostic_signature, Gene %in% inter.sect)
tmp <- gene_matrix[, prognostic_signature$Gene] %>% t() %>% as.data.frame()
colnames(tmp) <- gene_matrix$SampleID
gene_matrix <- tmp

# calculate RiskScore
RiskScore <- colSums(gene_matrix * prognostic_signature$beta)
RiskScore <- data.frame(SampleID=names(RiskScore), RiskScore=RiskScore)

#RiskScore <- RiskScore %>% filter(RiskScore != 0)

RiskScore$PatientID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\-"))[1]})
RiskScore$RegionID <- sapply(X=RiskScore$SampleID, FUN=function(x) {unlist(strsplit(as.character(x), split="\\-"))[2]})


# define riskscore_thresh
if(is.numeric(cut.off)){
  riskscore_thresh <- cut.off
  
}else{
  riskscore_thresh <- median(RiskScore$RiskScore)
}
# riskscore_thresh <- -5
 }

{
# classify patients as high-risk or low-risk
RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")

## get risk classes
risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
risk_class <- data.frame(High=as.matrix(risk_class)[,"High"], Low=as.matrix(risk_class)[,"Low"]) %>% tibble::rownames_to_column("PatientID")

# assign as "low", "discordant" or "high"
risk_class$class <- NA
risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
risk_class$class <- gsub(x=risk_class$class, pattern="NA", replacement="")
risk_class$class <- gsub(x=risk_class$class, pattern="HighLow", replacement="Discordant")

# join to RiskScore df
RiskScore <- dplyr::left_join(x=RiskScore, y=risk_class[,c("PatientID", "class")], by="PatientID")
RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))

# no. per class
tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()


# scatter-plot
sc.plot <- ggplot(RiskScore, aes(x=fct_reorder(PatientID, RiskScore + as.numeric(class), .fun=min), y=RiskScore)) 
sc.plot <- sc.plot + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(floor(min(RiskScore$RiskScore)), ceiling(max(RiskScore$RiskScore))))
sc.plot <- sc.plot + theme(axis.text.x = element_text(angle=90, vjust=0.5))
sc.plot <- sc.plot + geom_hline(yintercept = riskscore_thresh, col="black", lty="dotted")
sc.plot <- sc.plot + geom_line(col="black")
sc.plot <- sc.plot + ggtitle(label="", subtitle = paste0(length(unique(RiskScore$PatientID)), " ESCC patients = ", table(tmp$class)["Low"], " low + ", table(tmp$class)["High"], " high + ", table(tmp$class)["Discordant"], " discordant")) + theme(plot.title = element_text(hjust=0.5, face="bold")) + xlab("PatientID") + ylab(plot.title)
sc.plot <- sc.plot + theme(legend.position = "bottom") + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
sc.plot <- sc.plot + theme(aspect.ratio = 0.5)
sc.plot <- sc.plot + geom_point(pch=16, aes(col=class), size=4, alpha=0.5) + scale_color_manual(values = c("#3B4992FF", "azure4", "#EE0000FF")) + theme(legend.position = "none")
sc.plot }
ggsave('sc.plot_Sig1.pdf',plot = sc.plot,width = 6,height = 4)
# calculate percentages for risk classes
tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
tmp <- data.frame(table(tmp$class))
tmp$class <- tmp$Var1 %>% as.character()
tmp$class[c(1,3)] <- paste0("Concordant\n", tmp$class[c(1,3)], " Risk")
tmp <- tmp[c(1,3,2),]
tmp$class <- factor(tmp$class, levels=tmp$class)
tmp$Perc <- round(tmp$Freq/sum(tmp$Freq)*100) 


#  bar-plot
bar.plot <- ggplot(tmp, aes(x=class, y=Perc, fill=class)) + geom_bar(stat="identity")
bar.plot <- bar.plot + scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title=NULL))
bar.plot <- bar.plot + geom_text(size=5, aes(y = (Perc+2.5), label = paste0(Perc, "%"))) 
bar.plot <- bar.plot + theme_classic() + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")), axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),  axis.ticks.length = unit(-1.4, "mm")) + scale_y_continuous(expand = c(0,0), limits=c(0,70)) + theme(legend.position = "none", aspect.ratio = 1)
bar.plot <- bar.plot + xlab("") + ylab("Survival risk classification (%)")
bar.plot <- bar.plot + ggtitle(label=plot.title)
bar.plot
ggsave('bar.plot_Sig1.pdf',plot = bar.plot,width = 4,height = 4)

write.csv(RiskScore,"Riskscore_Sig2.csv")
