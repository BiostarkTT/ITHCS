library("survival") 
library("survminer")

SurvivalPlot <- function(survival.data, sample.class, filename, out.file.path){
 # Plot survival curve for different group of patients.
 #
 # Args:
 #   survival.data: the first three column of data frame must include patient, survival time and survival status.
 #   sample.class: the first two column of data frame must include patient and label.
 #   filename: the name of file, for instance 'cyp.system.mut.survival.pdf'.
 #   out.file.path: the path where file deposited.
 
 # Returns:
 #  null
 
 survival.data <- survival.data[, 1:3, drop = FALSE]
 colnames(survival.data) <- c('patient', 'times', 'status')
 
 sample.class <- sample.class[, 1:2, drop = FALSE]
 colnames(sample.class) <- c('patient', 'subtype')
 
 clinical.data <- merge(x = sample.class, y = survival.data, by = 'patient')
 
 cols <- c("#E7B800", "#2E9FDF", "#D1751E","#572875","#D8908F","#E0B665")[1:length(unique(clinical.data$subtype))]
 fit <- survfit(Surv(times, status) ~ subtype, data = clinical.data)
 
 
 sur.plot <- ggsurvplot(data = clinical.data, fit, pval = TRUE, conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata", linetype = "strata", 
                        surv.median.line = "hv", xscale = 365.25, break.time.by = 365.25, ggtheme = theme_bw(), 
                        palette = cols, xlab = 'Years', ylab = 'Survival probability')
 
 sur <- arrange_ggsurvplots(list(sur.plot), print = TRUE, ncol = 1, nrow = 1)
 
 ggsave(filename, plot = sur, path = out.file.path)
}

