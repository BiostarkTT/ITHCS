library(tidyverse)
library(data.table)
library(survival)
library(survminer)
TCGA = fread('TCGA-ESCC_Risk_input.tsv') %>% column_to_rownames('id') %>% as.data.frame()

TCGA$Alcohol = factor(TCGA$Alcohol,levels = c("No","Yes"))
TCGA$Gender = factor(TCGA$Gender,levels = c('female','male'))
TCGA$`TNM stage` = factor(TCGA$`TNM stage` ,levels= c('I','II','III',"IV"))

multiCox=coxph(Surv(OS.time, OS) ~ ., data = TCGA)
multiCoxSum=summary(multiCox)
multiCox$concordance
#????ģ????????Ϣ
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab <- as.data.frame(outTab)


pdf(file="forest_TCGA.pdf",
    width = 10,             #ͼƬ?Ŀ???
    height = 3,            #ͼƬ?ĸ߶?
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4),
         fontsize = 0.7,
         refLabel = "reference",
         noDigits = 2)
dev.off()
