# Load necessary libraries
library(tidyverse)
library(data.table)
library(reshape2)
library(ggpubr)

# 1. Load and process data
dat <- fread('./zhang_ICB.csv')
dat <- dat %>% select(1, 16, 4, 12) %>% column_to_rownames('Patient')

# Rename columns
colnames(dat)[1] <- "Type"

# 2. Convert data to long format for ggplot2
data <- melt(dat, id.vars = c("Type"))
colnames(data) <- c("Type", "Gene", "Expression")

# Sort data by Gene
data <- data %>% arrange(Gene)

# Ensure "Type" factor levels are ordered for plot consistency
data$Type <- factor(data$Type, levels = c('low', 'high'))

# 3. Plot boxplot comparing gene expression by type (low vs. high)
p <- ggboxplot(data, x = "Gene", y = "Expression", fill = "Type",
               ylab = "Gene expression", xlab = "",
               palette = c('navy', 'firebrick3'),  # Custom colors
               width = 0.6)

# Add statistical significance annotations
p1 <- p + stat_compare_means(aes(group = Type),
                             symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                symbols = c("***", "**", "*", " ")),
                             label = "p.signif")

# Display the final plot
p1
