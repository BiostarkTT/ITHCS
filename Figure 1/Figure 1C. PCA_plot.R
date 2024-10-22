library(tidyverse)
library(data.table)

# Load normalized expression data
normalised.vsd <- fread('GSE33426_exp.csv') %>% column_to_rownames('id')

# Select top 100 variably expressed genes based on standard deviation
exp.sd <- apply(normalised.vsd, 1, sd)
top.sd.genes <- names(exp.sd)[order(exp.sd, decreasing = TRUE)[1:100]]
normalised.vsd <- normalised.vsd[top.sd.genes, ]

# Load clinical data and select relevant columns
group <- fread('GSE33426_multiRegion_T.txt') %>% select(1, 3)

# Match expression data with sample IDs from the clinical data
normalised.vsd <- normalised.vsd[, group$id]

# Perform Principal Components Analysis (PCA)
pca.res <- prcomp(t(normalised.vsd), scale = TRUE)
pca.pc <- as.data.frame(pca.res$x)

# Prepare sample information and merge with PCA results
sam.info <- group %>% as.data.frame()
sam.info$Etiology <- "a"
colnames(sam.info) <- c('id', 'Patient', 'Etiology')
pca.pc$id <- rownames(pca.pc)
pca.pc <- inner_join(pca.pc, sam.info, by = 'id')

# Plot PCA results
library(ggpubr)
scatter.plot <- ggscatter(pca.pc, 'PC1', 'PC2', shape = 'Etiology', color = 'Patient')
scatter.plot
ggsave(scatter.plot, file = 'PCA100gene.pdf', width = 4, height = 4)
