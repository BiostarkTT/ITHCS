# Load necessary R packages
library(Seurat)          # For single-cell RNA sequencing analysis
library(SeuratData)      # Seurat data package, includes spatial transcriptomics datasets
library(ggplot2)         # For data visualization
library(patchwork)       # To combine multiple plots
library(dplyr)           # For data manipulation

# Load a standard spatial transcriptomics dataset "anterior1"
stdat <- readRDS("ESCC_st.RDS")

# Plot nCount_Spatial as a violin plot and spatial feature plot
plot1 <- VlnPlot(stdat, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()  # Violin plot, without legend
plot2 <- SpatialFeaturePlot(stdat, features = "nCount_Spatial") + theme(legend.position = "right")  # Spatial feature plot, with legend on the right
wrap_plots(plot1, plot2)  # Combine the two plots

# Apply SCT normalization (standardized process) to the spatial assay data
stdat <- SCTransform(stdat, assay = "Spatial", verbose = FALSE)

# Plot spatial expression of riskscore
SpatialFeaturePlot(stdat, features = c("riskscore"))

# Customize the legend by adjusting text size, title size, and key size
plot <- SpatialFeaturePlot(stdat, features = c("riskscore")) + 
  theme(legend.text = element_text(size = 0),      # Set legend text size to 0, effectively hiding it
        legend.title = element_text(size = 20),    # Set legend title font size
        legend.key.size = unit(1, "cm"))           # Set legend key size to 1 cm

# Save the plot as a JPEG file with specified dimensions and quality
jpeg(filename = "spatial_vignette_riskscore.jpg", height = 700, width = 1200, quality = 50)
print(plot)  # Print the plot
dev.off()    # Close the device and save the file

# Visualize spatial expression of the riskscore gene with varying point size and transparency
p1 <- SpatialFeaturePlot(stdat, features = "riskscore", pt.size.factor = 1)  # Set point size to 1
p2 <- SpatialFeaturePlot(stdat, features = "riskscore", alpha = c(0.1, 1))   # Set transparency range from 0.1 to 1
p1 + p2  # Combine the two plots for comparison

# Perform PCA (Principal Component Analysis) using the SCT normalized data
stdat <- RunPCA(stdat, assay = "SCT", verbose = FALSE)

# Construct a nearest neighbor graph based on the first 30 principal components
stdat <- FindNeighbors(stdat, reduction = "pca", dims = 1:30)

# Perform clustering based on the nearest neighbor graph
stdat <- FindClusters(stdat, verbose = FALSE)

# Run UMAP (Uniform Manifold Approximation and Projection) for dimensionality reduction based on the first 30 principal components
stdat <- RunUMAP(stdat, reduction = "pca", dims = 1:30)

# Plot UMAP dimensional reduction with cluster labels
p1 <- DimPlot(stdat, reduction = "umap", label = TRUE)  # UMAP plot with cluster labels
p2 <- SpatialDimPlot(stdat, label = TRUE, label.size = 3)  # Spatial dimensional plot with cluster labels, label size set to 3
p1 + p2  # Combine the UMAP and spatial dimensional plots
