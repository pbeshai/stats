library("ez")

util.pToString <- function (p) {
  if (p < 0.001) {
    pStr <- "p < .001"
  } else if (p < 0.01) {
    pStr <- "p < .01" 
  } else if (p < 0.05) {
    pStr <- "p < .05"
  } else {
    pStr <- paste("p =", round(p, 3))
  }
  return(pStr)
}

util.anovaBriefPrint <- function (briefResults) {
  writeLines("descriptive")
  print(briefResults$descriptive)
  
  writeLines("\nANOVA")
  writeLines(briefResults$anova)
  
  if ("posthoc" %in% names(briefResults)) {
    writeLines("\nPost-Hoc")
    print(briefResults$posthoc)
  }
}


util.anovaToString <- function (anova_results) {
  DFn <- anova_results$ANOVA$DFn
  DFd <- anova_results$ANOVA$DFd
  F <- round(anova_results$ANOVA$F, 3)
  ges <- round(anova_results$ANOVA$ges, 3)
  output <- ""
  
  # determine the p value
  mauchly <- anova_results$`Mauchly's Test for Sphericity`
  sphericity <- anova_results$`Sphericity Corrections`
  
  # check for sphericity violation
  if (mauchly$p<.05) {
    output <- paste0("Mauchly's test showed assumption of sphericity was violated: W(", 
                     anova_results$ANOVA$DFn, ") = ", round(mauchly$W, 3), ", ", 
                     util.pToString(mauchly$p), "\n")
    
    # check if we use Greenhouse-Geisser (epsilon < .75) or Huynh-Feldt correction
    if (sphericity$GGe < 0.75) {
      output <- paste0(output, "Using Greenhouse-Geisser correction.\n")
      p <- sphericity$`p[GG]`
    } else {
      output <- paste0(output, "Using Huynh-Feldt correction.\n")
      p <- sphericity$`p[HF]`
    }
  } else { # sphericity not violated, do not change p
    p <- anova_results$ANOVA$p
  }
  
  output <- paste0(output, "F(", DFn, ", ", DFd, ") = ", F, ", ", util.pToString(p), ", ges = ", ges)
  return(output)
}

util.mode <- function (x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

util.printHeader <- function (header) {
  writeLines(paste0("\n\n", header, "\n-------------------------------------------------------------------------------"))
}

util.printBigHeader <- function (header) {
  writeLines("\n===============================================================================")
  writeLines(header)
  writeLines("===============================================================================")
}

# convert from 1 row per participant to n rows per participant where n is the number of conditions
# This is required to run repeated measures/within-subjects ANOVA
util.withinSubjectsStack <- function (data, columns, participantColumn = "Participant", dvName = "value", ivName = "condition", participantName = "participant", omitNA=FALSE) {
  # Take a subset of the data only slow vs med vs fast user correct
  subset <- data[c(participantColumn, columns)]
  
  # cannot have NA values for Friedman test
  if (nrow(subset) > nrow(na.omit(subset))) {
    warning(paste0("Data contained NA values and was pruned. ", nrow(subset) - nrow(na.omit(subset)), " rows affected."))
    subset <- na.omit(subset)
  }
  
  # Create the participant column for each of the conditions (used when stacked)
  participant <- rep(subset[[participantColumn]], length(subset) - 1) # -1 to subtract participant column
  
  # stack the data for repeated-measures anova (1 row per condition)
  subset_stack <- stack(subset[2:length(subset)]) # ignore Participant column here
  subset_stack[3] <- participant # 1 = value (DV), 2 = condition (IV), 3 = participant
  
  # Name the data
  colnames(subset_stack) <- c(dvName, ivName, participantName)
  
  
  # ensure we have the participant factor set up properly in case we pruned some 
  # (otherwise get the error: not an unreplicated complete block design)
  subset_stack[[participantName]] <- factor(subset_stack[[participantName]])
  
  return(list(data=subset, stack=subset_stack))
}

util.withinSubjectsNonParametricAnalysis <- function (data, columns, participantColumn = "Participant", dvName = "value", ivName = "condition", participantName = "participant") {
  # do not use scientific notation
  saved_scipen = getOption("scipen")
  saved_digits = getOption("digits")
  options(scipen=100,digits=4)
  
  results_summary = list()
  
  # prepare the data
  prepared_data <- util.withinSubjectsStack(data, columns, participantColumn, dvName, ivName, participantName, omitNA=TRUE) 
  stacked_data <- prepared_data$stack
  subset_data <- prepared_data$data
    
  # summary stats
  util.printHeader("Summary Statistics")
  summary_results <- list(modes=apply(subset_data[2:length(subset_data)], 2, util.mode), summary=summary(subset_data[2:length(subset_data)]))
  
  writeLines("Modes")
  print(summary_results$modes)
  writeLines("")
  print(summary_results$summary)
  
  results_summary$modes = summary_results$modes
  
  # non-parametric anova equivalent: Friedman test
  util.printHeader("Friedman Rank Sum Test Results")
  friedman_results <- friedman.test(rating ~ condition | participant, data=stacked_data)
  
  # pretty print
  results_summary$friedman = paste0("Chi^2(", as.numeric(friedman_results$parameter), ") = ", round(as.numeric(friedman_results$statistic), 3), ", ", util.pToString(friedman_results$p.value))
  writeLines(results_summary$friedman)
  print(friedman_results)
  
  if (friedman_results$p.value > 0.05) {
    writeLines("==> Friedman test not significant.")
    posthoc_results <- NULL
  } else {
    # post-hoc tests
    util.printHeader("Post-hoc Test Results (Pairwise Wilcoxon Test with Holm correction)")
    posthoc_results <- pairwise.wilcox.test(stacked_data[[dvName]], stacked_data[[ivName]], p.adjust.method="holm", paired=T)
    print(posthoc_results)
    results_summary$posthoc <- posthoc_results$p.value
  }
  
  return(list(brief=results_summary, data=subset_data, stacked_data=stacked_data, summary=summary_results, friedman=friedman_results, posthoc=posthoc_results))
}


util.withinSubjectsAnalysis <- function (data, columns, participantColumn = "Participant", dvName = "value", ivName = "condition", participantName = "participant") {
  saved_scipen = getOption("scipen")
  saved_digits = getOption("digits")
  options(scipen=100,digits=4)
  
  results_summary = list()
  
  # prepare the data
  prepared_data <- util.withinSubjectsStack(data, columns, participantColumn, dvName, ivName, participantName) 
  stacked_data <- prepared_data$stack
  subset_data <- prepared_data$data
  
  
  util.printHeader("Summary Statistics")
  
  # summary_results <- ezStats(data=stacked_data, dv=.(dvName), wid=.(participantName), within=.(ivName))
  # sad work-around since ez package does strange eval of parameters
  summary_results <- eval(parse(text=paste0("ezStats(data=stacked_data, dv=", dvName, ", wid=", participantName, ", within=", ivName, ")")))
  print(summary_results)
  results_summary$descriptive <- summary_results
  
  # within-subjects ANOVA
  util.printHeader("ANOVA Results")
  #anova_results <- ezANOVA(data=stacked_data, dv=.(dvName), wid=.(participantName), within=.(ivName))
  # sad work-around since ez package does strange eval of parameters
  anova_results <- eval(parse(text=paste0("ezANOVA(data=stacked_data, dv=", dvName, ", wid=", participantName, ", within=", ivName, ")")))
  
  
  # pretty print the results
  results_summary$anova <- util.anovaToString(anova_results)
  writeLines(results_summary$anova)
  writeLines("\n")
  print(anova_results)
  
  # boxplot the data
  boxplot(value~condition,data=stacked_data)
  
  if (anova_results$ANOVA$p > 0.05) {
    writeLines("==> ANOVA not significant.")
    posthoc_results <- NULL
  } else {
    util.printHeader("Post-hoc Test Results (Pairwise t-Test with Holm correction)")
    posthoc_results <- pairwise.t.test(stacked_data[[dvName]], stacked_data[[ivName]], p.adjust.method="holm", paired=T)
    print(posthoc_results)
    results_summary$posthoc <- posthoc_results$p.value
  }
  
  # revert options for scientific notation
  options(scipen=saved_scipen,digits=saved_digits)
  rm(saved_scipen,saved_digits)
  
  return(list(brief=results_summary, plot=plot, data=subset_data, stacked_data=stacked_data, summary=summary_results, anova=anova_results, posthoc=posthoc_results))
}