source("util.R")

# Data for the same likert-style question asked after each block in the experiment
all_data <- read.csv("sample_data.csv")
columns <- c("likertA1","likertA2","likertA3")

util.printBigHeader("Running Non-parametric Within-subjects Analysis")
likert_results <- util.withinSubjectsNonParametricAnalysis(all_data, columns, dvName="rating")
