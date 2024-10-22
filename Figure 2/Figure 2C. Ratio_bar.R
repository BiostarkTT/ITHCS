library(tidyverse)
library(data.table)

# Load ssGSEA scores and set row names
rt <- fread('In_house_ssGSEA_scores.csv') %>% column_to_rownames('id')

# Select relevant columns: 'group_cut' and the variable of interest
variate <- "Type"
rt <- rt[, c('group_cut', variate)]

# Display cell counts for each group
table(rt$group_cut)
# Display proportion of each type
prop.table(table(rt$Type))
# Display cell counts for each group by type
table(rt$Type, rt$group_cut)

# Calculate proportions of each cell type within groups
Cellratio <- prop.table(table(rt$Type, rt$group_cut), margin = 2) %>% as.data.frame()
colnames(Cellratio)[1] <- variate
Cellratio$Percentage <- round(Cellratio$Freq * 100, 1)

# Define colors for the plot
allcolour <- c("#EF6563", "#6BBAAF", "#6CA3D3", "#8F708F", "#874942", "#7E59A8")

# Create bar plot
plot <- ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = `Type`), stat = "identity", width = 0.7, size = 0.5, colour = '#333333') + 
  geom_text(aes(x = Var2, y = Freq, label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white") +  
  theme_classic() +
  labs(x = 'Sample', y = 'Ratio') +
  scale_fill_manual(values = allcolour) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        legend.position = 'top')

# Display the plot
plot

# Save the plot as a PDF
ggsave("Type_subcluster_inhouse.pdf", plot, width = 3, height = 3.5)