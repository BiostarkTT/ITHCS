# Suppress messages from libraries
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))

# Load expression data and metadata
expr <- readr::read_rds(file = "./ESCC_196756.rds.gz")
meta <- readr::read_rds(file = "./ESCC_196756_metadata.rds.gz")

# Filter expression matrix based on minimum cells and genes
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

# Perform entropy analysis
ent.res <- SE_fun(expr)
head(ent.res)  # Display the first few rows of the entropy results

# Plot entropy results
SEplot(ent.res)

# Calculate ROGUE score (robustness of gene expression) based on entropy results
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value  # Output the ROGUE score

# Generate boxplot for ROGUE values
rogue.boxplot(rogue.value)
