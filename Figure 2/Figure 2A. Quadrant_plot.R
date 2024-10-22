# Load ITH score
load(file = './ESCC_intra.inter.ith.score.RData')

QuadrantPlot <- function(ith.score, file.name){
  # Ensure correct column names for intra and inter scores
  if (!all(c('losic.intra.sd', 'losic.inter.sd') %in% colnames(ith.score))) {
    stop("ITH score data should contain 'losic.intra.sd' and 'losic.inter.sd' columns")
  }
  
  # Rename columns for clarity
  colnames(ith.score) <- c('intra.score', 'inter.score')
  
  # Assign quadrants based on intra and inter scores
  within.ith <- ifelse(ith.score$intra.score > mean(ith.score$intra.score), 'top', 'bottom')
  between.ith <- ifelse(ith.score$inter.score > mean(ith.score$inter.score), 'right', 'left')
  ith.score$quadrant <- paste(within.ith, between.ith, sep = '_')
  
  # Table for quadrant scaling
  quadrant_count <- table(factor(ith.score$quadrant, levels = c("top_left", "bottom_left", "top_right", "bottom_right")))
  scaling_factors <- quadrant_count / min(quadrant_count)
  
  # Open PDF for saving plot
  pdf(paste0(file.name, '.pdf'))
  
  # Plot scaled circles for quadrants
  plot(c(100, 100, 400, 400), c(400, 100, 400, 100), pch = 16, cex = scaling_factors * 2, 
       ylim = c(0, 500), xlim = c(0, 500), col = c("firebrick1", "darkorchid2", "gold1", "turquoise2"),
       axes = FALSE, ylab = '', xlab = '')
  abline(h = 250, v = 250, col = 1, lty = 2, lwd = 4)
  text(c(100, 100, 400, 400), c(450, 200, 450, 200), paste(quadrant_count, 'genes'))
  box()
  
  # Density plots for 'inter.score'
  hist(ith.score$inter.score, breaks = 100, prob = TRUE, col = '#6baed699', axes = FALSE,
       ylab = '', xlab = '', main = 'Between tumour, LIHC', border = 'white')
  lines(density(ith.score$inter.score), lwd = 3, col = '#08519c')
  abline(v = mean(ith.score$inter.score), col = 1, lwd = 4, lty = 2)
  
  # Density plots for 'intra.score'
  hist(ith.score$intra.score, breaks = 100, prob = TRUE, col = '#fd8d3c99', axes = FALSE,
       ylab = '', xlab = '', main = 'Within tumour, LIHC', border = 'white')
  lines(density(ith.score$intra.score), lwd = 3, col = '#a63603')
  abline(v = mean(ith.score$intra.score), col = 1, lwd = 4, lty = 2)
  
  # Close PDF device
  dev.off()
}

# Run the function with appropriate data and file name
QuadrantPlot(losic.ith.score[, c('losic.intra.sd', 'losic.inter.sd')], 'losic_qu
