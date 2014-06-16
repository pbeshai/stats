source("util.R")

# parametric data taken from each block of the experiment (within-subjects)
all_data <- read.csv("sample_data.csv")
columns <- c("conditionA","conditionB","conditionC")

util.printBigHeader("Running Parametric Within-subjects Analysis on A vs B vs C");
results <- util.withinSubjectsAnalysis(all_data, columns)