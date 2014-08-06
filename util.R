library("ez")
library("coin")
library("ggplot2")
library("scales")
library("gridExtra")


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
  numDecimals <- 2
  DFn <- anova_results$ANOVA$DFn
  DFd <- anova_results$ANOVA$DFd
  F <- round(anova_results$ANOVA$F, numDecimals)
  ges <- round(anova_results$ANOVA$ges, numDecimals)
  output <- ""

  # determine the p value
  mauchly <- anova_results$`Mauchly's Test for Sphericity`
  sphericity <- anova_results$`Sphericity Corrections`

  # check for sphericity violation
  if (mauchly$p<.05) {
    output <- paste0("Mauchly's test showed assumption of sphericity was violated: W(",
                     anova_results$ANOVA$DFn, ") = ", round(mauchly$W, numDecimals), ", ",
                     util.pToString(mauchly$p), "\n")

    # check if we use Greenhouse-Geisser (epsilon < .75) or Huynh-Feldt correction
    if (sphericity$GGe < 0.75) {
      output <- paste0(output, "Using Greenhouse-Geisser correction.\n")
      p <- sphericity$`p[GG]`
      DFn <- DFn * sphericity$GGe
      DFd <- DFd * sphericity$GGe
    } else {
      output <- paste0(output, "Using Huynh-Feldt correction.\n")
      p <- sphericity$`p[HF]`
      DFn <- DFn * sphericity$HFe
      DFd <- DFd * sphericity$HFe
    }
  } else { # sphericity not violated, do not change p
    p <- anova_results$ANOVA$p
  }

  output <- paste0(output, "F(", round(DFn, numDecimals), ", ", round(DFd, numDecimals), ") = ", F, ", ", util.pToString(p), ", ges = ", ges)
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

util.descriptiveStats <- function (data, ignoreFirstColumn = FALSE) {
  if (ignoreFirstColumn) { # typically the Participants column
    data <- data[2:length(data)]
  }
  
  util.printHeader("Descriptive Statistics")
  descriptive <- list(modes=apply(data, 2, util.mode), 
                      summary=summary(data))
  
  writeLines("Modes")
  print(descriptive$modes)
  writeLines("")
  print(descriptive$summary)
  
  descriptive
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

  dput(subset_data)
  # summary stats
  descriptive <- util.descriptiveStats(subset_data, TRUE)

  results_summary$modes <- descriptive$modes
  results_summary$medians <- apply(subset_data[2:length(subset_data)], 2, median)

  # check if only two conditions, then only do wilcoxon test
  if (length(levels(stacked_data[[ivName]])) == 2) {
    writeLines("\nOnly 2 levels, using Wilcoxon test.\n")
    x<-subset_data[[levels(stacked_data$condition)[1]]]
    y<-subset_data[[levels(stacked_data$condition)[2]]]
    
    wilcoxon_results <- util.wilcoxonSignedRankTest(x,y)
    print(wilcoxon_results)
    
    # Pretty print
    writeLines(paste0("W = ", wilcoxon_results$statistic[[1]], ", Z = ", round(wilcoxon_results$statistic[[2]], 4), ", ", 
                      util.pToString(wilcoxon_results$p.value), ", r = ", round(wilcoxon_results$parameter[[2]], 4)))
    
    if (wilcoxon_results$p.value > 0.05) {
      writeLines("==> Wilcoxon test not significant.")
    } else {
      writeLines("Wilcoxon test significant.")
    }
    
    return(list(brief=results_summary, data=subset_data, stacked_data=stacked_data, descriptive=descriptive, wilcoxon=wilcoxon_results))
  }
  
  writeLines("\nMore than 2 levels, using Friedman test.\n");
  
  # non-parametric anova equivalent: Friedman test
  util.printHeader("Friedman Rank Sum Test Results")
  # sample formula: rating ~ condition | participant
  friedman_results <- friedman.test(as.formula(paste(dvName, "~", ivName, "|", participantName)), data=stacked_data)

  # pretty print
  results_summary$friedman = paste0("chi^2(", as.numeric(friedman_results$parameter), ") = ", round(as.numeric(friedman_results$statistic), 3), ", ", util.pToString(friedman_results$p.value))
  writeLines(results_summary$friedman)
  print(friedman_results)

  if (friedman_results$p.value > 0.05) {
    writeLines("==> Friedman test not significant.")
    posthoc_results <- NULL
  } else {
    # post-hoc tests
    util.printHeader("Post-hoc Test Results (Pairwise Wilcoxon Test with Bonferroni correction)")
    posthoc_results <- util.pairwise.wilcoxonSignedRankTest(stacked_data[[dvName]], stacked_data[[ivName]], p.adjust.method="bonferroni")
    print(posthoc_results)
    writeLines("Effect sizes (r)")
    print(posthoc_results$r)
    results_summary$posthoc <- list(p.value=posthoc_results$p.value, r=posthoc_results$r)
  }

  return(list(brief=results_summary, data=subset_data, stacked_data=stacked_data, descriptive=descriptive, friedman=friedman_results, posthoc=posthoc_results))
}


# Inspired by http://yatani.jp/HCIstats/WilcoxonSigned
# Can be used on Likert type responses, within-subjects.
# Provides effect size (r) and the W and Z statistics.
util.wilcoxonSignedRankTest <- function (x, y) {
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  # get the W statistic using wilcox.test. exact = F since we expect 
  # ties (two participants in same set with same value) and 
  # zeroes (one participant in x and y having the same value in both sets)
  W <- as.numeric(wilcox.test(x, y, paired = T, exact = F)$statistic)
  names(W) <- "W"
  
  # compute exact p value and get Z value using wilcoxsign_test from coin package
  results <- wilcoxsign_test(x ~ y, distribution="exact", zero.method="Pratt")
  p <- pvalue(results)
  
  Z <- as.numeric(statistic(results))
  names(Z) <- "Z"
  
  N <- 2*length(x)
  names(N) <- "N"
  
  r <- Z/sqrt(N) # effect size: 0.1 small, 0.3 medium, 0.5 large
  names(r) <- "r"
  
  # use htest to make it pretty
  ans <- list(statistic=c(W, Z), parameter=c(N, r), p.value=p, method="Wilcoxon Signed Rank Test with effect size (r)", data.name=DNAME)
  class(ans) <- "htest"
  ans
}

util.pairwise.wilcoxonSignedRankTest <- function (x, g, p.adjust.method = p.adjust.methods, 
          ...) 
{
  p.adjust.method <- match.arg(p.adjust.method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  METHOD <- "Wilcoxon signed rank test"
  test <- "123"
  compare.levels <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]
    util.wilcoxonSignedRankTest(xi, xj)$p.value
  }
  
  # lazy way to get all the effect size calculations
  compare.levelsR <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]
    util.wilcoxonSignedRankTest(xi, xj)$parameter[["r"]]
  }
  
  
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  RVAL <- pairwise.table(compare.levelsR, levels(g), "none")
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, r = RVAL,
              p.adjust.method = p.adjust.method, test = test)
  class(ans) <- "pairwise.htest"
  ans
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
    util.printHeader("Post-hoc Test Results (Pairwise t-Test with Bonferroni correction)")
    posthoc_results <- pairwise.t.test(stacked_data[[dvName]], stacked_data[[ivName]], p.adjust.method="bonferroni", paired=T)
    print(posthoc_results)
    results_summary$posthoc <- posthoc_results$p.value
  }

  # revert options for scientific notation
  options(scipen=saved_scipen,digits=saved_digits)
  rm(saved_scipen,saved_digits)

  return(list(brief=results_summary, plot=plot, data=subset_data, stacked_data=stacked_data, summary=summary_results, anova=anova_results, posthoc=posthoc_results))
}

util.posthocAnalysis <- function (data, dvName="value", ivName="condition", numTrials=NULL, participantName="participant", numDecimals=1, paired=F) {
  results <- list()
  aggr_data <- aggregate(as.formula(paste0(dvName, "~", participantName, "*", ivName)), data=data, FUN=sum)
  results$title <- paste(dvName, "vs", ivName);
  results$t.test <- pairwise.t.test(aggr_data[[dvName]], aggr_data[[ivName]], p.adjust.method="bonferroni", paired=paired)
  
  results$numTrials <- numTrials
  
  dvivFormula <- as.formula(paste0(dvName, "~", ivName))
  results$descriptive <- aggregate(dvivFormula, data=aggr_data, FUN=function (x) { round(mean(x), numDecimals) })
  results$descriptive[,3] <- aggregate(dvivFormula, data=aggr_data, FUN=function (x) { round(sd(x), numDecimals) })[,2]
  colnames(results$descriptive) <- c(ivName, "mean", "sd") 
  if (!is.null(numTrials)) {
    results$descriptive[,4] <- aggregate(dvivFormula, data=aggr_data, FUN=function (x) { paste0(round(100*mean(x)/numTrials, numDecimals), "%") })[,2]
    results$descriptive[,5] <- aggregate(dvivFormula, data=aggr_data, FUN=function (x) { paste0(round(100*sd(x)/numTrials, numDecimals), "%") })[,2]
    colnames(results$descriptive) <- c(ivName, "mean", "sd", "mean%", "sd%")
  } 
  results
}

util.withinSubjectsMultifactorAnalysis <- function (stacked_data, dvName = "value", ivs = "condition", participantName = "participant", numTrials=NULL, dvLabel="value", ivLabels=NULL) {
  # do not use scientific notation
  saved_scipen = getOption("scipen")
  saved_digits = getOption("digits")
  options(scipen=100,digits=4)
  
  util.printBigHeader("Running Within-Subjects Multifactor Analysis");
  if (is.null(ivLabels)) {
    ivLabels <- ivs
  }
  ivNames <- paste0(".(", paste(ivs, collapse=","), ")")
  
  options(contrasts=c("contr.sum", "contr.poly"))
  #anova_results <- ezANOVA(data=stacked_data, dv=.(dvName), wid=.(participantName), within=.(ivs))
  anova_results <- eval(parse(text=paste0("ezANOVA(data=stacked_data, dv=", dvName, ", wid=", participantName, ", within=", ivNames, ")")))
  
  
  # plot significant interaction effects
  interactions <- anova_results$ANOVA[grepl(":",anova_results$ANOVA[["Effect"]]),]
  sig_interactions <- interactions[interactions$p<.05,]
  
  if (nrow(sig_interactions) > 0) {
    apply(sig_interactions, 1, function (interaction) {
      interactionIVs <- strsplit(interaction[["Effect"]], ":")[[1]]
      if (length(interactionIVs) > 2) {
        # TODO this should probably handle the case whenthere are >2 IVs in the interaction
        warning(paste("More than 2 IVs in the interaction effect. Probably not plotting what you want.", interaction[["Effect"]]))
      }
      
      f <- as.formula(paste0("value ~ ",paste0(c(participantName, interactionIVs), collapse="*")))
      aggr_data <- aggregate(f, data=stacked_data, FUN=sum)
      
      # plot both ways in case one is more helpful than another
      interaction.plot(aggr_data[[interactionIVs[1]]], aggr_data[[interactionIVs[2]]], aggr_data[[dvName]],
                       xlab=ivLabels[[interactionIVs[1]]], ylab=dvLabel, trace.label=ivLabels[[interactionIVs[2]]])
      interaction.plot(aggr_data[[interactionIVs[2]]], aggr_data[[interactionIVs[1]]], aggr_data[[dvName]],
                       xlab=ivLabels[[interactionIVs[2]]], ylab=dvLabel, trace.label=ivLabels[[interactionIVs[1]]])
    })  
  }
  
  
  # for each significant main effect, run t-tests
  posthoc <- lapply(1:length(ivs), function (i) {
    iv <- ivs[[i]]
    if (is.null(numTrials)) {
      nt <- NULL
    } else {
      nt <- numTrials[[i]]  
    }
    
    # get the row and run the analysis if significant
    ivResults <- anova_results$ANOVA[which(anova_results$ANOVA$Effect == iv),]
    if (ivResults$p<.05) {
      f <- as.formula(paste0("value ~ ",paste0(c(participantName, iv), collapse="*")))
      aggr_data <- aggregate(f, data=stacked_data, FUN=sum)
     
      boxplot(as.formula(paste0("value ~ ", iv)), data=aggr_data, xlab=ivLabels[[iv]], ylab=dvLabel)
      
      return(util.posthocAnalysis(stacked_data, ivName=iv, numTrials=nt, paired=T))
    } 
    
    # Not significant:
    
    # lazy way to run this to get the descriptive stats.
    results <- util.posthocAnalysis(stacked_data, ivName=iv, numTrials=nt, paired=T)
    results$t.test <- FALSE
    results
  });
  #   print();
  #   print(util.posthocAnalysis(stacked_data, ivName="dSpeed", numTrials=144, paired=T));
  #   print(util.posthocAnalysis(stacked_data, ivName="dLocation", numTrials=216, paired=T));
  #   

  print(anova_results);
  print(posthoc);
  
  # revert options for scientific notation
  options(scipen=saved_scipen,digits=saved_digits)
  rm(saved_scipen,saved_digits)
  
  return(list(anova=anova_results, posthoc=posthoc))
}

util.barPlot <- function(data, x="condition", y="value", fill="type", xlabel=NULL, ylabel=NULL, legend_label=NULL, title=NULL, grouped=F) {
  # spacing between bars determined by difference from 1 of bar width
  barWidth <- .75
  
  barPosition <- "stack" #default
  if (grouped) {
    # space between groups determined by difference from 1 of bar width + bar spacing
    barSpacing <-.08
    
    barPosition <- position_dodge(barWidth + barSpacing)
  }
  
  if (nlevels(data[[fill]]) <= 2) {
    scaleColours <- c("#88CD7F", "#2c7fb8")
  } else if (nlevels(data[[fill]]) <= 4) {
    scaleColours <- c("#a1dab4", "#41b6c4",  "#2c7fb8", "#253494")
  } else {
    scaleColours <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac")
  }
  
  
  the_plot <- ggplot(data, aes_string(x=x, y=y, fill=fill)) + 
    # draw the bars with the legend and no outline to prevent diagonal line on legend colours
    geom_bar(position=barPosition, 
             stat="identity",
             width=barWidth) +
    
    # draw the bars with outline and no legend (draws over the bars with no outline)
    geom_bar(position=barPosition, 
             stat="identity", 
             width=barWidth,
             colour=alpha("black", .5),
             show_guide=F) +
    
    # add in error bars based on confidence interval (ci)
    geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=barPosition) +
    
    # add in axis labels and title
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title) +
    
    # legend 
    scale_fill_manual(name=legend_label, values=scaleColours) + 
    
    # y-axis. expand=c(0,0) prevents padding at bottom, limit=c(0,1) ensures top stays at 100%
    scale_y_continuous(labels = percent_format(), breaks=seq(0,1,.2), expand=c(0,0), limit=c(0,1)) +
    
    theme_bw() +
    theme(axis.title.y=element_text(vjust=1),  # adjust axis and title spacing
          axis.title.x=element_text(vjust=-.5), 
          plot.title=element_text(vjust=1, face="bold"),
          # hide vertical grid lines
          panel.grid.major.x=element_blank(),
          # add a nice outline around the legend colours
          legend.key = element_rect(colour=alpha("black", .75)))
  
  the_plot
}

util.frequencyBarPlot <- function(data, x="condition", y="value", fill="type", xlabel=NULL, ylabel=NULL, legend_label=NULL, title=NULL, grouped=F, ylimit=NULL) {
  # spacing between bars determined by difference from 1 of bar width
  barWidth <- .75
  
  barPosition <- "stack" #default
  if (grouped) {
    # space between groups determined by difference from 1 of bar width + bar spacing
    barSpacing <-.08
    
    barPosition <- position_dodge(barWidth + barSpacing)
  }
  
  if (nlevels(data[[fill]]) <= 2) {
    scaleColours <- c("#88CD7F", "#2c7fb8")
  } else if (nlevels(data[[fill]]) <= 4) {
    scaleColours <- c("#a1dab4", "#41b6c4",  "#2c7fb8", "#253494")
  } else {
    scaleColours <- c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac") # http://colorbrewer2.org
    scaleColours <- c("#FED9D9", "#FC9191",  "#7bccc4", "#43a2ca", "#0868ac") # red at bottom 2
  }
  
  the_plot <- ggplot(data, aes_string(x=x, y=y, fill=fill)) + 
    # draw the bars with the legend and no outline to prevent diagonal line on legend colours
    geom_bar(position="dodge", 
             stat="identity",
             width=barWidth) +
    
    # draw the bars with outline and no legend (draws over the bars with no outline)
    geom_bar(position="dodge", 
             stat="identity", 
             width=barWidth,
             colour=alpha("black", .5),
             show_guide=F) +
        
    # add in axis labels and title
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title) +
    
    # legend 
    scale_fill_manual(name=legend_label, values=scaleColours) + 
    
    # y-axis. expand=c(0,0) prevents padding at bottom, limit=c(0,1) ensures top stays at 100%
#     scale_y_discrete(expand=c(0,0)) +
  scale_y_continuous(labels = percent_format(), breaks=seq(0,1,.2), expand=c(0,0), limit=ylimit) +
  
    
    theme_bw() +
    theme(axis.title.y=element_text(vjust=1),  # adjust axis and title spacing
          axis.title.x=element_text(vjust=-.5), 
          plot.title=element_text(vjust=1, face="bold"),
          # hide vertical grid lines
          panel.grid.major.x=element_blank(),
          # add a nice outline around the legend colours
          legend.key = element_rect(colour=alpha("black", .75)))
  
  the_plot
}





####### FROM: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

