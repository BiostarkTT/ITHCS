library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(limma)
library(tidyverse)
library(dplyr)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(data.table)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
rt1 <- fread('./GSE53625_train.csv') %>% as.data.frame()
rt2 <- fread('./TCGA_test.csv')%>% as.data.frame()
rt3 <- fread('./zhang_test.csv')%>% as.data.frame()
# 生成包含三个数据集的列表
mm <- list(GSE53625 = rt1,
           TCGA = rt2, zhang = rt3)

result <- data.frame()
# GSE53625作为训练集
est_data <- mm$GSE53625
# GEO作为验证集
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
# est_data <- est_data %>% as.data.frame()
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]

val_dd_list <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS',pre_var)]})

# 设置种子数和节点数，其中节点数可以调整
rf_nodesize <- 5
seed <- 100

##################################
#### 1-1.RSF ####
##################################
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)



##################################
#### 2.Enet ####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##################################
#### 3.StepCox ####
##################################
for (direction in c("both", "backward", "forward")) {
  # direction = "forward"
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(val_dd_list, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ ALPL+ARL8A+BIK+C3orf80+CAP1+CCND1+COX10+CT83+EIF4H+FBXL19+FOSB+GSTO2+HSD17B2+IL24+KLHDC8B+LIPE+LRRN4CL+MED31+MEN1+MFHAS1+MPC1+MPDZ+MTRF1+NAA16+POLR3D+PSME2+RPL36AL+RPL3L+SNAI2+SRSF9+SUGCT+SYT15+TMEM223+TRIM4+TTC17+TTC39C+UBL4A+ZDHHC4+ZNF251, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both", "backward", "forward")) {
  # direction = "forward"
  fit <- step(coxph(Surv(OS.time, OS)~., est_dd), direction = direction)
  rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因
  est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 500, parallel = T)
  cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime=x[, 1], newstatus=x[,2], type="lp")))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ ., x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
    cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  set.seed(seed)
  cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
  fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time,
                 event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1,2)])))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  result <- rbind(result, cc)
  x1 <- as.matrix(est_dd2[, rid])
  x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "cox", alpha = 0)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$OS.time,
               censoring.status = est_dd2$OS,
               featurenames = colnames(est_dd2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = 5,
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  rs <- lapply(val_dd_list2, function(w){
    test <- list(x = t(w[, -c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2], RS = rr)
    return(rr2)
  })
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}

##################################
#### 4-1.CoxBoost ####
##################################
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd[, 'OS.time'], est_dd[, 'OS'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

##################################
#### 5.plsRcox####
##################################
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$OS.time, status = est_dd$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$OS.time, event = est_dd$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

##################################
#### 6.superpc####
##################################
data <- list(x = t(est_dd[, -c(1,2)]), y = est_dd$OS.time, censoring.status = est_dd$OS, featurenames = colnames(est_dd)[-c(1, 2)])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)

##################################
#### 7.GBM####
##################################
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time, OS)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)

##################################
#### 8.survival-SVM####
##################################
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

##################################
#### 9.Ridge####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "cox")

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cvfit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time,OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

##################################
####10.Lasso####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 5, #例文描述：10-fold cross-validation
                family = 'cox', type.measure = 'C')
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

##################################
#### 10.1.Lasso + CoxBoost####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

##################################
#### 10.2.Lasso + GBM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

##################################
#### 10.3.Lasso + plsRcox####
##################################

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

##################################
#### 10.4.Lasso + RSF####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid<-rid[-1]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

##################################
#### 10.5.Lasso + stepcox####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[, c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

##################################
#### 10.6.Lasso + superPC####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', rid)]})
data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$OS.time, censoring.status = est_dd2$OS,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$OS.time, censoring.status = w$OS, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

##################################
#### 10.7.Lasso + survival-SVM####
##################################
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "cox", alpha = 1, type.measure = 'C')
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
fit = survivalsvm(Surv(OS.time,OS)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)

#将得到的结果赋给result2变量进行操作
result2 <- result

###将结果的长数据转换为宽数据
dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()
# write.table(dd2,"output_C_test_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#将C指数定义为数值型
dd2[,-1] <- apply(dd2[,-1], 2, as.numeric)
#求每个模型的C指数在三个数据集的均值
dd2$All <- apply(dd2[,2:4], 1, mean)
#求每个模型的C指数在GEO验证集的均值
dd2$GEO <- apply(dd2[,3:4], 1, mean)
###查看每个模型的C指数
head(dd2)

#输出C指数结果
write.table(dd2,"output_C_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)

dd2=read.table('output_C_test_index.txt', header=T, sep="\t", check.names=F)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
# 根据C指数排序
dd2 <- dd2[order(dd2$All, decreasing = T),]
# 仅绘制GEO验证集的C指数热图
dt <- dd2[, 1:5]
rownames(dt) <- dd2$Model

##热图绘制
Cindex_mat=dt
Cindex_mat <- Cindex_mat[,-1]
avg_Cindex <- apply(Cindex_mat, 1, mean)           # 计算每种算法在所有队列中平均C-index
avg_Cindex <- sort(avg_Cindex, decreasing = T)     # 对各算法C-index由高到低排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      # 对C-index矩阵排序

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
CohortCol <- allcolour[1:4]
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)


cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c("#807DBA", "#FFFFE0", "#FFB6C1"), # 黄绿配色
              #col = c("#4195C1", "#FFFFFF", "#CB5746"), # 红蓝配色
              rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
              cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
             # width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
             # height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path( "Cindex_test.pdf"), width = 6.65, height = 8.5)
draw(hm)
invisible(dev.off())
