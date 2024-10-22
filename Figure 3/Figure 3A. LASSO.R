# Load necessary libraries
library("glmnet")
library("survival")

# Load the expression data
rt <- read.table("uniSigExp_GSE53625.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Prepare input data for Lasso-Cox model
x <- as.matrix(rt[, 3:ncol(rt)])  # Independent variables (gene expression data)
y <- data.matrix(Surv(rt$futime, rt$fustat))  # Survival data (futime and fustat)

# Fit a Lasso-Cox model
fit <- glmnet(x, y, family = "cox")
plot(fit, xvar = "lambda", label = TRUE)  # Plot model with respect to lambda

# Perform cross-validation to select optimal lambda value
set.seed(8)
cvfit <- cv.glmnet(x, y, type.measure = 'C', nfolds = 5, family = "cox")
plot(cvfit)
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")  # Highlight optimal lambdas

# Extract non-zero coefficients at lambda.min
coef <- coef(cvfit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]  # Extract active coefficients
lassoGene <- row.names(coef)[index]  # Gene names corresponding to active coefficients

# Include futime and fustat columns for output
lassoGene <- c("fustat", "futime", lassoGene)
lassoExp <- rt[, lassoGene]
lassoExp <- cbind(id = row.names(lassoExp), lassoExp)

# Save the results to a file
write.table(lassoExp, file = "lassoExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)