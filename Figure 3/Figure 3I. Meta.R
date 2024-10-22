# Load necessary library
library(meta)

# Load all files with ".cox.txt" in their name from the working directory
files <- dir(pattern = ".cox.txt")  # Retrieve only files with '.cox.txt' in their name

# Combine all files into a single data frame
data <- data.frame()
for (i in files) {
  rt <- read.table(i, header = TRUE, sep = "\t", check.names = FALSE)
  data <- rbind(data, rt)
}

# Alternatively, load pre-saved HR data
# data <- data.table::fread('HR_all.txt')

# Extract necessary columns for meta-analysis
study <- data$id
HR <- data$HR
lower.HR <- data$HR.95L
upper.HR <- data$HR.95H

# Perform meta-analysis (fixed-effects model)
meta <- metagen(log(HR),
                lower = log(lower.HR),
                upper = log(upper.HR),
                studlab = study,
                sm = "HR",
                comb.random = FALSE,    # TRUE for random-effects model
                comb.fixed = TRUE)      # TRUE for fixed-effects model
print(meta)

# Plot the forest plot
pdf(file = "forest.pdf", width = 10, height = 5)
forest(meta,
       col.square = "black",           # Color of squares
       col.diamond = "red",            # Color of diamond
       col.diamond.lines = "black")    # Color of diamond lines
dev.off()
