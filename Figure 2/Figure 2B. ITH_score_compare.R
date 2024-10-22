
##  scatter plot and correlation function
corr_plot_DB <- function(data, x, y, point_colour="deepskyblue", colour_var=NA, 
 point_border_col = "deepskyblue", title="", x_axis="x", y_axis="y", corr_method="spearman", 
 pch_setting=21, point_size=5, point_alpha=0.5, x_limits = NA, y_limits = NA, aspect_ratio=1, best_fit_line=FALSE, 
 line_colour="black", line_type="dashed") {
  library(ggplot2)
  #import
  corr <- data
  
  # spearman correlation
  if(corr_method=="spearman") {
    #calculate spearman's rho
    gg_cor <- cor.test(as.numeric(corr[,x]), as.numeric(corr[,y]), method = "spearman")$estimate %>% as.numeric() %>% signif(digits=3)
    #calculate spearman's p-value  
    gg_cor_p <- cor.test(corr[,x], corr[,y], method = "spearman")$p.value %>% as.numeric() %>% signif(digits=3)
  }
  # spearman correlation
  if(corr_method=="pearson") {
    #calculate PMCC
    gg_cor <- cor.test(as.numeric(corr[,x]), as.numeric(corr[,y]), method = "pearson")$estimate %>% as.numeric() %>% signif(digits=3)
    #calculate pearson's p-value  
    gg_cor_p <- cor.test(corr[,x], corr[,y], method = "pearson")$p.value %>% as.numeric() %>% signif(digits=3)
  }
  
  # fix P-value if v small
  if(gg_cor_p == 0) {gg_cor_p <- "< 2.2e-16"}
  
  # draw plot
  gg <- ggplot(corr, aes_string(x, y)) 
  # add points
  if(is.na(colour_var)) { 
    gg <- gg + geom_point(pch=pch_setting, col=point_border_col, fill=point_colour, size=point_size, alpha=point_alpha) 
  }
  if(!is.na(colour_var)) { 
    gg <- gg + geom_point(aes_string(col=colour_var), size=point_size, alpha=point_alpha) 
  }
  
  # add themes
  gg <- gg + theme_classic() 
  gg <- gg + theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) #border
  gg <- gg + theme(aspect.ratio = aspect_ratio)
  
  # add spearman/pearson correlation outputs
  if(corr_method=="spearman") {
    #gg <- gg + ggtitle(label=title, subtitle=paste("Spearman's Rho = ", gg_cor, "\nSpearman's P-value = ", gg_cor_p, sep=""))
    subtitle <- paste0("Rs = ", gg_cor, "\nP = ", gg_cor_p)
  }
  if(corr_method=="pearson") { 
    #gg <- gg + ggtitle(label=title, subtitle=paste("PMCC = ", gg_cor, "\nPearson's P-value = ", gg_cor_p, sep=""))  
    subtitle <- paste0("PMCC = ", gg_cor, "\nP = ", gg_cor_p)
  }
  
  # add title and subtitle
  gg <- gg + ggtitle(label=title, subtitle=subtitle)  
  
  # add best-fit line
  if(best_fit_line) {
    gg <- gg + geom_smooth(method = "lm", se = F, colour=line_colour, linetype=line_type)
  }
  
  gg <- gg + xlab(x_axis) + ylab(y_axis) + theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  if(!is.na(x_limits)) {gg <- gg + xlim(c(x_limits[1], x_limits[2]))}
  if(!is.na(y_limits)) {gg <- gg + ylim(c(y_limits[1], y_limits[2]))}
  
  return(gg)
}


#load(file='./data/intra.inter.ith.score.RData')

## plot
# ITH correlation claculated by using losic data set
{
  
  # correlate gene-wise MAD scores with SD scores
  p1 <- corr_plot_DB(data = losic.ith.score, x = "losic.intra.sd", y = "losic.intra.mad", title = "", point_size = 2, 
   x_axis = "Standard deviation", y_axis = "Median absolute deviation", best_fit_line = T, aspect_ratio = NULL)
  
  # correlate gene-wise CV scores with SD scores
  
  p2 <- corr_plot_DB(data = losic.ith.score, x = "losic.intra.sd", y = "losic.intra.cv", title = "", point_size = 2, 
   x_axis = "Standard deviation", y_axis = "Coefficient of variation", best_fit_line = T, aspect_ratio = NULL)
   
  # correlate gene-wise CV scores with MAD scores 
  p3 <- corr_plot_DB(data = losic.ith.score, x = "losic.intra.cv", y = "losic.intra.mad", title = "", point_size = 2, 
   x_axis = "Coefficient of variation", y_axis = "Median absolute deviation", best_fit_line = T, aspect_ratio = NULL)
   
}


{
  
  # correlate gene-wise MAD scores with SD scores
  p4 <- corr_plot_DB(data = losic.ith.score, x = "losic.inter.sd", y = "losic.inter.mad", title = "", point_size = 2, 
   x_axis = "Standard deviation", y_axis = "Median absolute deviation", best_fit_line = T, aspect_ratio = NULL)
  
  # correlate gene-wise CV scores with SD scores
  
  p5 <- corr_plot_DB(data = losic.ith.score, x = "losic.inter.sd", y = "losic.inter.cv", title = "", point_size = 2, 
   x_axis = "Standard deviation", y_axis = "Coefficient of variation", best_fit_line = T, aspect_ratio = NULL)
   
  # correlate gene-wise CV scores with MAD scores 
  p6 <- corr_plot_DB(data = losic.ith.score, x = "losic.inter.cv", y = "losic.inter.mad", title = "", point_size = 2, 
   x_axis = "Coefficient of variation", y_axis = "Median absolute deviation", best_fit_line = T, aspect_ratio = NULL)
   
}


# # ITH correlation between losic and shi data set
# {
#   
#   # correlate gene-wise MAD scores with SD scores
#   p7 <- corr_plot_DB(data = ith.score, x = "losic.intra.sd", y = "shi.intra.sd", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\nintra ith scores(SD)", y_axis = "SHI\nintra ith scores(SD)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#   
#   # correlate gene-wise CV scores with SD scores 
#   p8 <- corr_plot_DB(data = ith.score, x = "losic.intra.mad", y = "shi.intra.mad", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\nintra ith scores(MAD)", y_axis = "SHI\nintra ith scores(MAD)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#    
#   # correlate gene-wise CV scores with MAD scores 
#   p9 <- corr_plot_DB(data = ith.score, x = "losic.intra.cv", y = "shi.intra.cv", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\nintra ith scores(CV)", y_axis = "SHI\nintra ith scores(CV)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#    
# }
# 
# 
# {
#   
#   # correlate gene-wise MAD scores with SD scores
#   p10 <- corr_plot_DB(data = ith.score, x = "losic.inter.sd", y = "shi.inter.sd", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\ninter ith scores(SD)", y_axis = "SHI\ninter ith scores(SD)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#   
#   # correlate gene-wise CV scores with SD scores 
#   p11 <- corr_plot_DB(data = ith.score, x = "losic.inter.mad", y = "shi.inter.mad", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\ninter ith scores(MAD)", y_axis = "SHI\ninter ith scores(MAD)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#    
#   # correlate gene-wise CV scores with MAD scores 
#   p12 <- corr_plot_DB(data = ith.score, x = "losic.inter.cv", y = "shi.inter.cv", corr_method = "pearson", point_size = 2, 
#    x_axis = "LOSIC\ninter ith scores(CV)", y_axis = "SHI\ninter ith scores(CV)", 
#    best_fit_line = T, point_colour = "gold", point_border_col = "gold")
#    
# }

library(ggpubr)
setwd('/result/Section2/')
ggsave(ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2), file='ith_correlation_losic.pdf')

ggsave(ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2), file='ith_correlation_between_losic_shi.pdf')


# quadrant gene compare
{
within.ith <- ifelse(ith.score$losic.intra.sd > mean(ith.score$losic.intra.sd), 'top', 'bottom')
between.ith <- ifelse(ith.score$losic.inter.sd > mean(ith.score$losic.inter.sd), 'right', 'left')
ith.score$losic.quadrant.sd <- paste(within.ith, between.ith, sep='_')

within.ith <- ifelse(ith.score$losic.intra.mad > mean(ith.score$losic.intra.mad), 'top', 'bottom')
between.ith <- ifelse(ith.score$losic.inter.mad > mean(ith.score$losic.inter.mad), 'right', 'left')
ith.score$losic.quadrant.mad <- paste(within.ith, between.ith, sep='_')

within.ith <- ifelse(ith.score$losic.intra.cv > mean(ith.score$losic.intra.cv), 'top', 'bottom')
between.ith <- ifelse(ith.score$losic.inter.cv > mean(ith.score$losic.inter.cv), 'right', 'left')
ith.score$losic.quadrant.cv <- paste(within.ith, between.ith, sep='_')

write.csv(ith.score,"ith_Q4gene.csv")

within.ith <- ifelse(ith.score$shi.intra.sd > mean(ith.score$shi.intra.sd), 'top', 'bottom')
between.ith <- ifelse(ith.score$shi.inter.sd > mean(ith.score$shi.inter.sd), 'right', 'left')
ith.score$shi.quadrant.sd <- paste(within.ith, between.ith, sep='_')

within.ith <- ifelse(ith.score$shi.intra.mad > mean(ith.score$shi.intra.mad), 'top', 'bottom')
between.ith <- ifelse(ith.score$shi.inter.mad > mean(ith.score$shi.inter.mad), 'right', 'left')
ith.score$shi.quadrant.mad <- paste(within.ith, between.ith, sep='_')

within.ith <- ifelse(ith.score$shi.intra.cv > mean(ith.score$shi.intra.cv), 'top', 'bottom')
between.ith <- ifelse(ith.score$shi.inter.cv > mean(ith.score$shi.inter.cv), 'right', 'left')
ith.score$shi.quadrant.cv <- paste(within.ith, between.ith, sep='_')

}

# table(ith.score[, c('losic.quadrant.sd', 'shi.quadrant.sd')])
                 # shi.quadrant.sd
# losic.quadrant.sd bottom_left bottom_right top_left top_right
     # bottom_left         4988          805      718       529
     # bottom_right         672          212      221       255
     # top_left             500          196      236       388
     # top_right            736          440      570      1931




