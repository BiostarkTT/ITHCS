# Loading necessary libraries
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggsci)

# Load the dataset for AUC
load("C—index.Rdata")

# Select the cohort of interest
cohort <- rownames(AUC[['ESCC']]$CAD.Sig)
auc <- AUC[['ESCC']][1:6]

# Extract AUC data for different cohorts
ls_AUC <- lapply(cohort, function(c) {
  cc <- lapply(auc, function(auc) { auc[c, 'AUC'] }) %>% do.call(rbind, .)
  names(cc) <- names(auc)
  cc <- as.data.frame(cc)
  colnames(cc) <- 'AUC'
  return(cc)
}) 

# Circos plot function for comparing cohorts
circos_cormpare <- function(ls_AUC) {
  bardata <- do.call(cbind, ls_AUC) %>% as.data.frame() %>% t() %>% as.data.frame()
  bardata$cohort <- rownames(bardata)
  bardata_m <- melt(bardata, value.name = 'y', variable.name = 'x')
  colnames(bardata_m)[1] <- 'sectors'
  bardata_m$y[is.na(bardata_m$y)] <- 0
  
  # Color palette setup
  Set2 <- brewer.pal(8, "Set2")
  Set1 <- brewer.pal(9, "Set1")
  col <- c(Set2, Set1)
  col[1:2] <- c("#CA6855", '#546B7E')
  
  df <- bardata_m
  x1 <- rep(1:(ncol(bardata)-1), nrow(bardata))
  x1 <- x1[order(x1)]
  df$x1 <- x1
  df$sectors <- factor(df$sectors, levels = unique(df$sectors))
  
  circos.clear()
  circos.par(gap.degree = c(rep(2, nrow(bardata)-1), 15), cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0.01), start.degree = 90)
  circos.initialize(df$sectors, xlim = c(0, ncol(bardata)))
  
  bgcol <- c('#1A1A1A', rep('#999999', 5))
  circos.track(df$sectors, y = df$y, track.height = 0.03, bg.col = bgcol, bg.border = 'white', panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(3), CELL_META$sector.index, facing = 'inside', niceFacing = TRUE, cex = 1, font = 2)
  })
  
  circos.track(df$sectors, ylim = c(0, 0.85), track.height = 0.78, bg.border = 'white')
  for (i in unique(df$sectors)) {
    circos.barplot(df[df$sectors == i, 'y'], pos = df[df$sectors == i, 'x1'], sector.index = i, col = col[1:7])
  }
  
  circos.yaxis(side = 'left', at = seq(0, 1, 0.2), labels = TRUE, sector.index = unique(df$sectors)[1], track.index = 2)
  draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, rou1 = 0.15, col = "black", border = "#EEEEEE")
  text(0, 0, "AUC", col = 'white', cex = 1.5, font = 2)
  
  sig <- col[1:(ncol(bardata)-1)]
  names(sig) <- unique(df$x)
  
  lgd_cancer <- Legend(title = "Signature", at = names(sig), legend_gp = gpar(fill = sig))
  draw(lgd_cancer, x = unit(1.3, "snpc"), just = "left")
}

# Generate and save the circos plot
pdf(file = 'Figure_CAD_circos_compare.pdf', height = 6, width = 12)
circos_cormpare(ls_AUC)
dev.off()

###====Fig 5B - Heatmap====####

# Prepare data for the heatmap
hdata <- ls_AUC %>% do.call(cbind, .) %>% `colnames<-`(names(ls_AUC))
hdata$MeanAUC <- rowMeans(hdata)
hdata <- hdata[order(hdata[,1], decreasing = TRUE), ]

# Generate and save the heatmap
pdf(file = 'Figure_CAD_Heatmap.pdf', height = 4, width = 4)
Heatmap(as.matrix(hdata), name = "AUC",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(0.5, 0.80), c("#FFFFE3", "navy")),
        rect_gp = gpar(col = "black", lwd = 5),
        width = nrow(hdata) * unit(6, "mm"),
        height = nrow(hdata) * unit(6, "mm"),
        column_split = c('Combined', rep('Individual', ncol(hdata) - 1))
)
dev.off()