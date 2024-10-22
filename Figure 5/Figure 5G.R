rm(list = ls())  # Clear environment
library(Seurat)
library(tidyverse)

### 1. Read and prepare the data
lung_T <- readRDS('scRNA_ESCC.RDS')

# Create Seurat object
scRNA <- CreateSeuratObject(counts = lung_T, project = "scRNA", min.cells = 3, min.features = 200)

### 2. Quality Control (QC) ###
# Calculate percentage of mitochondrial genes
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

# Calculate percentage of red blood cell-related genes
HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]] <- PercentageFeatureSet(scRNA, features = HB.genes)

# Violin plot before QC
col.num <- length(levels(scRNA@active.ident))
violin <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), 
                  cols = rainbow(col.num), pt.size = 0.01, ncol = 4) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
violin

# Scatter plots for QC checks
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow = 1, legend = "none")
pearplot

### 3. Data Filtering ###
# Set QC thresholds
minGene <- 200
maxGene <- 3000
pctMT <- 20
pctHB <- 3

# Filter cells based on QC criteria
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

# Violin plot after QC
violin <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), 
                  cols = rainbow(col.num), pt.size = 0.1, ncol = 4) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
violin

### 4. Data Normalization ###
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

### 5. Dimensionality Reduction and Clustering ###
# Identify the top 3000 variable genes
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(scRNA), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 2.5)
plot <- CombinePlots(plots = list(plot1, plot2), legend = "bottom")
plot

# Scale the data for all genes
scale.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

# Perform PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot1 <- DimPlot(scRNA, reduction = "pca", group.by = "sample")
plot2 <- ElbowPlot(scRNA, ndims = 20, reduction = "pca")
plotc <- plot1 + plot2
plotc

# Select principal components (PCs) based on the ElbowPlot
pc.num <- 1:10
scRNA <- FindNeighbors(scRNA, dims = pc.num)

# Perform clustering
scRNA <- FindClusters(scRNA)
table(scRNA@meta.data$seurat_clusters)

# Save cluster information
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID = rownames(metadata), cluster_ID = metadata$seurat_clusters)
write.csv(cell_cluster, 'cluster/cell_cluster.csv', row.names = FALSE)

### 6. UMAP and t-SNE Visualization ###
scRNA <- RunUMAP(scRNA, dims = pc.num)
scRNA <- RunTSNE(scRNA, dims = pc.num)

# Save UMAP embeddings
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap, 'cluster/embed_umap.csv')

# Plot UMAP with cluster labels
plot2 <- DimPlot(scRNA, reduction = "umap", label = TRUE)
ggsave("cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)

### 7. Cell Type Identification ###
# Define marker genes for immune and other cell types (from Habermann et al.)
habermann_imm <- c('CD274', "CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", 
                   "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", 
                   "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", 
                   "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", 
                   "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", 
                   "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")

# Plot marker expression for other cell types
DotPlot(scRNA, features = habermann_oth, group.by = "seurat_clusters") + coord_flip()

# Plot marker expression for immune cells
DotPlot(scRNA, features = habermann_imm, group.by = "seurat_clusters") + coord_flip()

# Assign cell types based on marker expression (manual annotation needed)
# scRNA[['celltype']] <- unname(cluster_celltype[scRNA@meta.data$seurat_clusters])

# Visualize cell types with UMAP
DimPlot(scRNA, reduction = 'umap', group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()