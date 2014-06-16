source("util.R")

# Data for the same likert-style question asked after each block in the experiment
all_data <- read.csv("sample_data.csv")
columns <- c("likertA1","likertA2","likertA3")

util.printBigHeader("Running Non-parametric Within-subjects Analysis")
userDifficulty_results <- util.withinSubjectsNonParametricAnalysis(all_data, columns, dvName="rating")


# Note that running this code will often produce warnings from R about ties and
# zero values preventing exact computation of a p-value. The tie warning 
# happens when multiple participants in a block give the same response. 
# The zero warning happens when a single participant gives the same response 
# in two blocks (e.g., P1 slow = 4, P1 med = 4). It seems many people ignore 
# these warnings and still use the results.
