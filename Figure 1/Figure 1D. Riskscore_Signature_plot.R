library(forcats)
library(ggplot2)
library(tidyverse)
library(data.table)
source(file = 'Published_Signature_Formula.R')

# Load expression and clinical data
dat <- data.table::fread('GSE33426_exp.csv') %>% column_to_rownames('id')
group <- data.table::fread('GSE33426_multiRegion_T.txt')
dat <- dat[, group$id]
colnames(dat) <- group$sample

# Define signature genes and coefficients
clonal.sig <- data.frame(Sig1, Sig1_coef)
colnames(clonal.sig) <- c('Gene', "beta")
prognostic_signature <- clonal.sig

# Prepare the expression matrix
gene_matrix <- t(dat) %>% as.data.frame()
gene_matrix$SampleID <- rownames(gene_matrix)
gene_matrix <- gene_matrix %>% select(SampleID, everything())

# Define cut-off for risk classification
cut.off <- "median"
plot.title <- 'Clonal-based signature'
file.name <- 'clonal_risk_sig_losic'

# Select signature genes present in the dataset
inter.sect <- intersect(colnames(gene_matrix), prognostic_signature$Gene)
if(length(prognostic_signature$Gene) > length(inter.sect)){
  print('The following genes are missed:')
  print(setdiff(prognostic_signature$Gene, inter.sect))
}

# Subset signature to available genes and calculate the risk score
prognostic_signature <- subset(prognostic_signature, Gene %in% inter.sect)
tmp <- gene_matrix[, prognostic_signature$Gene] %>% t() %>% as.data.frame()
colnames(tmp) <- gene_matrix$SampleID
gene_matrix <- tmp
RiskScore <- colSums(gene_matrix * prognostic_signature$beta)
RiskScore <- data.frame(SampleID = names(RiskScore), RiskScore = RiskScore)

# Extract PatientID and RegionID from SampleID
RiskScore$PatientID <- sapply(X = RiskScore$SampleID, FUN = function(x) { unlist(strsplit(as.character(x), split = "\\-"))[1] })
RiskScore$RegionID <- sapply(X = RiskScore$SampleID, FUN = function(x) { unlist(strsplit(as.character(x), split = "\\-"))[2] })

# Define threshold for risk classification
if(!is.numeric(cut.off)){
  riskscore_thresh <- median(RiskScore$RiskScore)
}
write.csv(RiskScore, "riskscore_Sig1.csv")
RiskScore <- fread('riskscore_Sig1.csv')

# Classify patients as high or low risk
RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")

# Create risk class based on the presence of both high and low scores
risk_class <- table(RiskScore$PatientID, RiskScore$RiskScore_bin)
risk_class <- data.frame(High = as.matrix(risk_class)[, "High"], Low = as.matrix(risk_class)[, "Low"]) %>%
  tibble::rownames_to_column("PatientID")

# Assign patient risk class: "Low", "High", or "Discordant"
risk_class$class <- NA
risk_class$class <- ifelse(risk_class$High > 0, paste(risk_class$class, "High", sep=""), risk_class$class)
risk_class$class <- ifelse(risk_class$Low > 0, paste(risk_class$class, "Low", sep=""), risk_class$class)
risk_class$class <- gsub(x = risk_class$class, pattern = "NA", replacement = "")
risk_class$class <- gsub(x = risk_class$class, pattern = "HighLow", replacement = "Discordant")

# Merge risk class with RiskScore
RiskScore <- dplyr::left_join(x = RiskScore, y = risk_class[, c("PatientID", "class")], by = "PatientID")
RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))

# Scatter plot of risk scores
sc.plot <- ggplot(RiskScore, aes(x = fct_reorder(PatientID, RiskScore + as.numeric(class), .fun = min), y = RiskScore)) +
  theme_classic() + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(size = 8, margin = unit(c(0, 2.5, 0, 0), "mm")),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, margin = unit(c(2.5, 0, 0, 0), "mm")),
        axis.ticks.length = unit(-1.4, "mm")) +
  scale_y_continuous(expand = c(0, 0), limits = c(floor(min(RiskScore$RiskScore)), ceiling(max(RiskScore$RiskScore)))) +
  geom_hline(yintercept = riskscore_thresh, col = "black", lty = "dotted") +
  geom_line(col = "black") +
  geom_point(pch = 16, aes(col = class), size = 1.8, alpha = 0.5) +
  scale_color_manual(values = c("#3B4992FF", "azure4", "#EE0000FF")) +
  theme(legend.position = "none", aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("PatientID") + ylab(plot.title)

# Save the scatter plot
ggsave('sc.plot_Sig1.pdf', plot = sc.plot, width = 3.3, height = 4)

# Calculate percentages for risk classes
tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
tmp <- data.frame(table(tmp$class))
tmp$class <- tmp$Var1 %>% as.character()
tmp$class[c(1, 3)] <- paste0("Concordant\n", tmp$class[c(1, 3)], " Risk")
tmp <- tmp[c(1, 3, 2), ]
tmp$class <- factor(tmp$class, levels = tmp$class)
tmp$Perc <- round(tmp$Freq / sum(tmp$Freq) * 100)

# Bar plot for risk classification
bar.plot <- ggplot(tmp, aes(x = class, y = Perc, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#3B4992FF", "#EE0000FF", "gray75"), guide = guide_legend(title = NULL)) +
  geom_text(size = 5, aes(y = (Perc + 2.5), label = paste0(Perc, "%"))) +
  theme_classic() +
  theme(axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text.y = element_text(margin = unit(c(0, 2.5, 0, 0), "mm")),
        axis.text.x = element_text(margin = unit(c(2.5, 0, 0, 0), "mm")),
        axis.ticks.length = unit(-1.4, "mm"),
        legend.position = "none", aspect.ratio = 1) +
  xlab("") + ylab("Survival risk classification (%)") + ggtitle(plot.title)

# Save the bar plot
ggsave('bar.plot_Sig1.pdf', plot = bar.plot, width = 3, height = 3.5)