# Load necessary libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggsci)

# Load the C-index data and format it
dat <- fread('C-index_all.txt') %>% 
  as.data.frame() %>% 
  column_to_rownames('id')

# Prepare data for circos plot
bardata <- as.data.frame(t(dat))  # Transpose the data
bardata$cohort <- rownames(bardata)  # Add cohort column

# Reshape the data to long format for plotting
bardata_m <- melt(bardata, value.name = 'y', variable.name = 'x')
colnames(bardata_m)[1] <- 'sectors'
bardata_m$y[is.na(bardata_m$y)] <- 0  # Replace NA values with 0

# Set up colors for the plot
Set2 <- brewer.pal(7, "Set2")
Set1 <- brewer.pal(8, "Set1")
col <- c('#CA6855', '#546B7E', Set2, Set1)

# Generate x-axis positions for each sector
x1 <- rep(1:(ncol(bardata)-1), nrow(bardata))
x1 <- x1[order(x1)]
bardata_m$x1 <- x1
bardata_m$sectors <- factor(bardata_m$sectors, levels = unique(bardata_m$sectors))

# Create circos plot
pdf(file = 'Circos_C-Index.pdf', height = 6, width = 12)
circos.clear()

circos.par(gap.degree = c(rep(2, nrow(bardata)-1), 15), cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.01, 0.01), start.degree = 90)
circos.initialize(bardata_m$sectors, xlim = c(0, ncol(bardata)))

# Set background colors for each sector
bgcol <- c('#1A1A1A', rep('#999999', 5))

# Track for sector labels
circos.track(bardata_m$sectors, y = bardata_m$y, track.height = 0.03, bg.col = bgcol, bg.border = 'white',
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(3), 
                           CELL_META$sector.index,
                           facing = 'inside',
                           niceFacing = TRUE,
                           cex = 1,
                           font = 2)
             })

# Main barplot track
circos.track(bardata_m$sectors, ylim = c(0, 0.85), track.height = 0.78, bg.border = 'white')
for (i in unique(bardata_m$sectors)) {
  circos.barplot(bardata_m[bardata_m$sectors == i, 'y'], pos = bardata_m[bardata_m$sectors == i, 'x1'], sector.index = i, col = col[1:14])
}

# Y-axis on the left for the first sector
circos.yaxis(side = 'left', at = seq(0, 1, 0.2), labels = TRUE, sector.index = unique(bardata_m$sectors)[1], track.index = 2)

# Draw inner sector for AUC label
draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, rou1 = 0.15, col = "black", border = "#EEEEEE")
text(0, 0, "AUC", col = 'white', cex = 1.5, font = 2)

# Create legend for signatures
sig <- col[1:(ncol(bardata)-1)]
names(sig) <- unique(bardata_m$x)
lgd_cancer <- Legend(title = "Signature", at = names(sig), legend_gp = gpar(fill = sig))
draw(lgd_cancer, x = unit(1.3, "snpc"), just = "left")

dev.off()

# === Heatmap for C-index ===

# Load the C-index data again
hdata <- fread('C-index_all.txt') %>% 
  as.data.frame() %>% 
  column_to_rownames('id')

# Sort data by the first column in descending order
hdata <- hdata[order(hdata[,1], decreasing = TRUE), ]

# Generate heatmap for C-index values
pdf(file = 'C_index_mean.pdf', height = 20, width = 4)

Heatmap(as.matrix(hdata), name = "AUC",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(0.5, 0.65, 0.80), c("#1f78b4", "#FFFFE3", "#f76262")),
        rect_gp = gpar(col = "black", lwd = 5),
        width = nrow(hdata) * unit(6, "mm"),
        height = nrow(hdata) * unit(6, "mm"),
        column_split = c('Combined', rep('Individual', c(ncol(hdata)-1)))
)

dev.off()
