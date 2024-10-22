library(forcats)
library(ggplot2)
library(tidyverse)
library(data.table)
source(file = 'Published_Signature_Formula.R')

{
  # Load expression and clinical data
  dat <- data.table::fread('GSE33426_exp.csv') %>% column_to_rownames('id')
  group <- data.table::fread('GSE33426_multiRegion_T.txt')
  dat <- dat[, group$id]
  colnames(dat) <- group$sample
  
  # Define clonal signature and coefficients
  clonal.sig <- data.frame(Sig13, Sig13_coef)
  colnames(clonal.sig) <- c('Gene', "beta")
  prognostic_signature <- clonal.sig
  
  a = colnames(dat)
  
  # Prepare gene expression matrix
  gene_matrix <- t(dat) %>% as.data.frame()
  rownames(gene_matrix) <- a
  gene_matrix$SampleID <- rownames(gene_matrix)
  gene_matrix <- gene_matrix %>% select(SampleID, everything())
  
  # Define cut-off for risk classification
  cut.off <- "median"
  
  # Filter the signature genes present in the dataset
  inter.sect <- intersect(colnames(gene_matrix), prognostic_signature$Gene)
  if (length(prognostic_signature$Gene) > length(inter.sect)) {
    print('The following genes are missed:')
    print(setdiff(prognostic_signature$Gene, inter.sect))
  }
  
  # Subset the prognostic signature and calculate RiskScore
  prognostic_signature <- subset(prognostic_signature, Gene %in% inter.sect)
  tmp <- gene_matrix[, prognostic_signature$Gene] %>% t() %>% as.data.frame()
  colnames(tmp) <- gene_matrix$SampleID
  gene_matrix <- tmp
  
  RiskScore <- colSums(gene_matrix * prognostic_signature$beta)
  RiskScore <- data.frame(SampleID = names(RiskScore), RiskScore = RiskScore)
  
  # Extract PatientID and RegionID from SampleID
  RiskScore$PatientID <- sapply(X = RiskScore$SampleID, FUN = function(x) { unlist(strsplit(as.character(x), split = "\\-"))[1] })
  RiskScore$RegionID <- sapply(X = RiskScore$SampleID, FUN = function(x) { unlist(strsplit(as.character(x), split = "\\-"))[2] })
  
  # Define risk score threshold
  if (!is.numeric(cut.off)) {
    riskscore_thresh <- median(RiskScore$RiskScore)
  }
  
  # Classify patients as high-risk or low-risk
  RiskScore$RiskScore_bin <- ifelse(RiskScore$RiskScore > riskscore_thresh, "High", "Low")
  
  # Get risk classes
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
  RiskScore <- dplyr::left_join(RiskScore, risk_class[, c("PatientID", "class")], by = "PatientID")
  RiskScore$class <- factor(RiskScore$class, levels = c("Low", "Discordant", "High"))
  
  # Calculate percentages for risk classes
  tmp <- RiskScore %>% dplyr::select(PatientID, class) %>% dplyr::distinct()
  tmp <- data.frame(table(tmp$class))
  tmp$class <- tmp$Var1 %>% as.character()
  tmp$class[c(1,3)] <- paste0("Concordant\n", tmp$class[c(1,3)], " Risk")
  tmp <- tmp[c(1, 3, 2), ]
  tmp$class <- factor(tmp$class, levels = tmp$class)
  tmp$Perc <- round(tmp$Freq / sum(tmp$Freq) * 100)
}

# Generate risk class data for multiple signatures
signature_names <- paste0("Sig ", 1:13)
merge <- bind_rows(lapply(1:13, function(i) {
  tmp$SigName <- signature_names[i]
  tmp
}))

# Bar plot for survival risk classification across signatures
library(ggpubr)
bar.plot <- ggbarplot(data = merge, x = 'SigName', y = 'Perc', fill = 'class', color = 'class', label = TRUE, 
                      lab.col = "white", lab.pos = "in", x.text.angle = 45, 
                      palette = c("#3B4992FF", "#EE0000FF", "gray75"), lab.nb.digits = 0, 
                      xlab = FALSE, ylab = "Survival risk classification (%)")

# Save bar plot as PDF
pdf("MultiRisk_barplot.pdf", width = 5.5, height = 4)
print(bar.plot)
dev.off()