library(tidyverse)
library(data.table)

# Load expression data
normalised.vsd <- fread('GSE33426_exp.csv') %>% column_to_rownames('id')

# Top 100 variably expressed genes based on standard deviation
exp.sd <- apply(normalised.vsd, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing = TRUE)[1:100]]
normalised.vsd <- normalised.vsd[top.sd.genes, ]

# Load clinical data
group <- fread('GSE33426_multiRegion_T.txt') %>% select(1, 3)
colnames(group)[1:2] <- c('SampleID', 'PatientID')

# Match expression data with clinical sample IDs
normalised.vsd <- normalised.vsd[, group$SampleID]

# Transpose matrix and calculate Z-scores
mat1 <- normalised.vsd %>% t() %>% as.data.frame()
mat1 <- as.data.frame(sapply(X = mat1, FUN = scale))
rownames(mat1) <- rownames(normalised.vsd)

# Prepare annotation for the heatmap (dummy Etiology data added)
mat2 <- group
mat2$Etiology <- "a"
mat2 <- spread(mat2, SampleID, Etiology)
mat2 <- t(mat2[, -1]) %>% as.data.frame()
colnames(mat2) <- group$PatientID

# Align the order of mat1 and mat2
mat2 <- mat2[rownames(mat1), ]

# Define hierarchical clustering functions
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(dist(t(mat1)))
mat_cluster_rows <- sort_hclust(hclust(dist(mat1)))

# Plot heatmap 1 (expression heatmap)
ht1 <- Heatmap(
  matrix = as.matrix(mat1), 
  name = "ht1", 
  col = colorRamp2(seq(-3, 3), viridis(7, option = "magma")),
  column_dend_height = unit(30, "mm"),
  row_dend_width = unit(10, "mm"),
  row_dend_reorder = rev(mat_cluster_rows$order), 
  column_dend_reorder = rev(mat_cluster_cols$order),
  show_column_names = TRUE,
  show_row_names = TRUE,
  width = 4,
  heatmap_legend_param = list(title = NULL, color_bar = "continuous")
)
ht1

# Plot heatmap 2 (annotation heatmap)
colors <- structure(c("green4", "darkorchid3", "darkorange"), names = c("a", "b", "c"))
ht2 <- Heatmap(
  matrix = as.matrix(mat2),
  col = colors,
  name = "ht2",
  na_col = "white",
  rect_gp = gpar(col = "white"),
  show_row_names = TRUE, 
  show_column_names = TRUE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(title = NULL, color_bar = "discrete"),
  column_names_gp = gpar(fontsize = 4),
  width = 1
)
ht2

# Save the heatmaps to a PDF file
pdf('RNA-ITH_heatmap.pdf')
print(ht1 + ht2)

# Add borders to heatmaps
decorate_heatmap_body("ht1", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
decorate_heatmap_body("ht2", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1))})
dev.off()
