###################### Risk Plot Generation ######################
# Clear environment variables
rm(list = ls())   
options(stringsAsFactors = FALSE)

# 1. Load necessary libraries and data
library(pheatmap)
library(tidyverse)
library(data.table)

# 1.1 Load the dataset
rt <- fread("GSE53625_riskscore.pdf") %>% column_to_rownames('id')
rt <- rt[order(rt$ITHRS),]  # Sort by risk score (ITHRS)
str(rt)

# 2. Generate Risk Curve
point <- rt$ITHRS
range(point)

pdf("GSE53625_riskscore.pdf", width = 5.5, height = 3)
colors <- colorRampPalette(c("#013565", "firebrick2"))(length(point))

# Plot risk scores
plot(point, type = "p", pch = 20, cex = 0.55, xlab = "Patients (increasing risk score)", ylab = "log2(1+ITHRS Score)", col = colors)
abline(h = median(point), v = round(nrow(rt) / 2), lty = 2)  # Add median and midpoint lines
dev.off()

# 3. Generate Survival Status Plot
pdf("Test_risk_median.pdf", width = 10, height = 5)

# Assign colors based on survival status (OS)
color <- ifelse(rt$OS == 1, "firebrick2", "#013565")
table(color)

# Plot survival times
plot(rt$OS.time, pch = 19, xlab = "Patients (increasing risk score)", ylab = "Survival time (years)", col = color)
abline(v = round(nrow(rt) / 2), lty = 2)  # Add line separating risk groups
legend("topleft", c("Dead", "Alive"), bty = "n", pch = 19, col = c("firebrick2", "#013565"), cex = 1.2)
dev.off()

# 4. Generate Risk Heatmap
# 4.1 Extract expression data (excluding survival columns)
rt1 <- rt[c(3:(ncol(rt)-2))] %>% t() %>% as.data.frame()

# 4.2 Create annotation for risk categories
annotation <- data.frame(type = rt[, ncol(rt)])
rownames(annotation) <- rownames(rt)

# 4.3 Draw heatmap
ann_colors <- list(risk = c("high" = 'firebrick2', "low" = "#013565"))
bk <- c(seq(-3, -0.1, by = 0.01), seq(0, 3, by = 0.01))

pdf('Test_riskplot.pdf', width = 6, height = 4)
pheatmap(rt1, 
         annotation_colors = ann_colors,
         annotation = annotation,
         scale = "row",
         cluster_cols = FALSE,
         fontsize_row = 11,
         show_colnames = FALSE,
         fontsize_col = 5,
         color = c(colorRampPalette(colors = c("#013565", "white"))(length(bk) / 2), 
                   colorRampPalette(colors = c("white", "firebrick2"))(length(bk) / 2)),
         breaks = bk)
dev.off()