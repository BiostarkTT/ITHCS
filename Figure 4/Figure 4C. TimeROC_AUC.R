# Load necessary libraries
library(tidyverse)
library(data.table)
library(ggplot2)

# Load the GEO dataset and create the bar plot
dat1 <- fread('../GEO/GEO_timeROC.csv')
dat1$type <- factor(dat1$type, levels = paste0("Sig", seq(1:14)))

p1 <- ggplot(dat1, aes(x = type, y = AUC, fill = year)) + 
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() +
  labs(x = "Signature Type", y = "AUC", fill = "Year")  # Add axis labels and legend title

p1
# ggsave("GEO_TimeROC.pdf", width = 12, height = 5, plot = p1)

# Load the TCGA dataset and create the bar plot
dat2 <- fread('../TCGA/TCGA_timeROC.csv')
dat2$type <- factor(dat2$type, levels = c("ITHR", paste0("Sig", seq(1:13))))

p2 <- ggplot(dat2, aes(x = type, y = AUC, fill = year)) + 
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() +
  labs(x = "Signature Type", y = "AUC", fill = "Year")  # Add axis labels and legend title

p2
# ggsave("TCGA_TimeROC.pdf", width = 12, height = 5, plot = p2)

# Load the Zhang dataset and create the bar plot
dat3 <- fread('../Zhang/Zhang_timeROC.csv')
dat3$type <- factor(dat3$type, levels = paste0("Sig", seq(1:14)))

p3 <- ggplot(dat3, aes(x = type, y = AUC, fill = year)) + 
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() +
  labs(x = "Signature Type", y = "AUC", fill = "Year")  # Add axis labels and legend title

p3
# ggsave("Zhang_TimeROC.pdf", width = 12, height = 5, plot = p3)

# Create a line plot for AUC over time for TCGA data
ggplot(dat2, aes(x = year, y = AUC, group = type, color = type)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Year", y = "AUC", color = "Signature Type") +
  scale_x_discrete(limits = c("1y", "2y", "3y", "4y", "5y"))  # Ensure correct year order

# Advanced line plot with customized colors
colors <- c('#CA6855','#546B7E', palette("Set2"), palette("Set1"))

ABC <- ggplot(dat2, aes(x = year, y = AUC, group = type, color = type)) +
  geom_line(size = 1) +  # Set line size
  geom_point(size = 3, shape = 21, fill = "white") +  # Customize points
  scale_color_manual(values = setNames(colors, c("ITHR", paste0("Sig", seq(1:13))))) +
  theme_bw() +  # Use white background theme
  theme(
    text = element_text(size = 12),  # Set text size
    plot.title = element_text(face = "bold", hjust = 0.5),  # Title styling
    plot.subtitle = element_text(hjust = 0.5),  # Subtitle styling
    legend.position = "bottom",  # Position legend
    legend.box = "horizontal",  # Horizontal legend box
    axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt X-axis labels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  labs(
    title = "Signature AUC Over Time in TCGA-ESCC cohort",  # Plot title
    subtitle = "AUC values for different signatures across 5 years",  # Plot subtitle
    x = "Year",  # X-axis label
    y = "AUC",  # Y-axis label
    color = "Signature"  # Legend title
  ) +
  scale_x_discrete(limits = c("1y", "2y", "3y", "4y", "5y"))  # Ensure correct year order

ABC
# Save the advanced plot
pdf("TCGA_timeROC.pdf", width = 4, height = 6)
print(ABC)
dev.off()