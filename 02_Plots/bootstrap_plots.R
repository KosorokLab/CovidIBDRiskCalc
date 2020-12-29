################################################################################
##  Plot wrappers
##
################################################################################

#' Create Coefficient Summary Plots, either box plots or plots showing the proportions
#' for the estimated coefficient sign distributions. 
CreateCoefSummaryPlots <- function(contrast.matrix,
                                   hosp.coef,
                                   icu.coef,
                                   death.coef,
                                   plot.type = "sign",
                                   top.n.effects.flag = FALSE,
                                   top.n.cutoff = 10,
                                   sort.coef.plots = FALSE,
                                   var.axis.label = NULL,
                                   plot.lower.lim = NULL,
                                   plot.upper.lim = NULL){
  plot.type <- tolower(plot.type)
  stopifnot(plot.type %in% c("sign", "box"))
  
  hosp_contrasts <- t(contrast.matrix) %*% hosp.coef
  icu_contrasts <- t(contrast.matrix) %*% icu.coef 
  death_contrasts <- t(contrast.matrix) %*% death.coef
  
  # Long data frames for plotting
  hosp_cont_long <- CreateCoefLongDfForPlotting(hosp_contrasts)
  icu_cont_long <- CreateCoefLongDfForPlotting(icu_contrasts)
  death_cont_long <- CreateCoefLongDfForPlotting(death_contrasts)
  
  if (is.null(plot.lower.lim)) {
    lower_lim <- min(c(min(hosp_cont_long$Coefficient),
                       min(icu_cont_long$Coefficient),
                       min(death_cont_long$Coefficient)))
  } else lower_lim <- plot.lower.lim
  
  if (is.null(plot.upper.lim)) {
    upper_lim <- max(c(max(hosp_cont_long$Coefficient),
                       max(icu_cont_long$Coefficient),
                       max(death_cont_long$Coefficient)))
  } else upper_lim <- plot.upper.lim
  
  
  if (top.n.effects.flag == TRUE) {
    hosp_cont_long_vars <- hosp_cont_long %>% group_by(Variable) %>% 
      summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
      arrange(desc(AbsMeanCoef)) %>% 
      slice_head(., n = top.n.cutoff) %>% 
      pull(Variable)
    
    icu_cont_long_vars <- icu_cont_long %>% group_by(Variable) %>% 
      summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
      arrange(desc(AbsMeanCoef)) %>% 
      slice_head(., n = top.n.cutoff) %>% 
      pull(Variable)
    
    death_cont_long_vars <- death_cont_long %>% group_by(Variable) %>% 
      summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
      arrange(desc(AbsMeanCoef)) %>% 
      slice_head(., n = top.n.cutoff) %>% 
      pull(Variable)
    
    hosp_cont_long <- hosp_cont_long %>% filter(Variable %in% hosp_cont_long_vars)
    icu_cont_long <- icu_cont_long %>% filter(Variable %in% icu_cont_long_vars)
    death_cont_long <- death_cont_long %>% filter(Variable %in% death_cont_long_vars)
  }
  
  # Create a single data frame with all the variables to facilitate using grid_arrange
  long_w_outcomes <- bind_rows(hosp_cont_long %>% mutate(Outcome = "Hospitalization+"),
                               icu_cont_long %>% mutate(Outcome = "ICU+"),
                               death_cont_long %>% mutate(Outcome = "Death"),) %>% 
    mutate_at(., .vars = "Outcome", ~factor(., levels = c("Hospitalization+",
                                                          "ICU+",
                                                          "Death"),
                                            ordered = TRUE))
  
  
  if (plot.type == "sign"){
    coef_sign_plots <- CreateCoefSignPlot(coef.long =  long_w_outcomes, 
                                          sort.plot = sort.coef.plots,
                                          all.outcomes = TRUE)
    
    if (is.null(var.axis.label) == FALSE){
      coef_sign_plots <- coef_sign_plots + labs(x = var.axis.label)
    }
    
    arranged_plots <- coef_sign_plots
  }
  
  
  if (plot.type == "box"){
    
    coef_box_plots <- CreateCoefBoxplot(long_w_outcomes, 
                                        sort.plot = sort.coef.plots,
                                        all.outcomes = TRUE) +
      lims(y = c(lower_lim, upper_lim)) 
    
    
    if (is.null(var.axis.label) == FALSE){
      coef_box_plots <- coef_box_plots + labs(x = var.axis.label)
    }
    
    arranged_plots <- coef_box_plots
  }
  
  
  return(arranged_plots)
}

################################################################################
##  Create Plots from a list of bootstrap model fits
##
################################################################################

#' Create plots summarizing the coefficient results across boootstrap replicates
#' Currently includes two plots:
#' 1) A boxplot summarizing the distribution of the fitted coefficients for each covariate
#' 2) A bar chart summarizing the proportion of time each variable is equal to zero
CreateBootCoefPlots <- function(boot.list, 
                                lambda.to.use = "lambda.min", 
                                nice.var.key = NULL,
                                nonzero.threshold = NULL){
  stopifnot(is.null(nonzero.threshold) | (nonzero.threshold >= 0 & nonzero.threshold <= 1))
  boot_coefs <- map(boot.list, ~as.matrix(coef(., s = lambda.to.use)))
  
  boot_coef_long <- tibble(Variable = rep(rownames(boot_coefs[[1]]),
                                          times = length(boot.list)),
                           Coefficient = unlist(boot_coefs)) %>% 
    filter(Variable != "(Intercept)")
  
  if (is.null(nice.var.key) == FALSE) {
    boot_coef_long <- left_join(boot_coef_long, nice.var.key, by = c("Variable" = "VarName")) %>% 
      mutate(AllNamesNice = if_else(is.na(NiceName), Variable, NiceName)) %>% 
      select(-Variable, -NiceName) %>% 
      rename(Variable = AllNamesNice)
  }
  
  boot_coef_zero_sum <- boot_coef_long %>% 
    group_by(Variable) %>% 
    summarise(ProbZero = sum(Coefficient == 0)/length(Coefficient))
  
  if (is.null(nonzero.threshold) == FALSE){
    variables_above_threshold <- boot_coef_zero_sum %>% 
      filter(1 - ProbZero >= nonzero.threshold) %>% 
      pull(Variable)
    
    boot_coef_zero_sum <- boot_coef_zero_sum %>% 
      filter(Variable %in% variables_above_threshold)
    
    boot_coef_long <- boot_coef_long %>% 
      filter(Variable %in% variables_above_threshold)
  }
  
  coef_box <- ggplot(data = boot_coef_long, 
                     aes(x = fct_reorder(Variable, Coefficient, .fun = median), 
                         y = Coefficient)) +
    geom_boxplot() + 
    labs(y = "Log Odds",
         x = "Contrast") + 
    coord_flip()
  
  
  coef_zero_bar <- ggplot(data = boot_coef_zero_sum, 
                          aes(x = fct_reorder(Variable, ProbZero, .fun = median,
                                              .desc = TRUE), 
                              y = ProbZero)) +
    geom_col() + 
    labs(y = "Probability Coefficient = 0",
         x = "Contrast") +
    coord_flip()
  
  b_boot_reps <- nrow(boot_coef_long)/length(unique(boot_coef_long$Variable))
  
  boot_coef_sign_sum <- boot_coef_long %>% 
    group_by(Variable) %>% 
    summarise(nZero = sum(Coefficient == 0),
              nPositive = sum(Coefficient > 0),
              nNegative = sum(Coefficient < 0),
              SortHelperPos = sum(Coefficient > 0),
              SortHelperZero = sum(Coefficient == 0),
              SortHelperNeg = sum(Coefficient < 0),
              SortHelper = b_boot_reps*SortHelperPos - SortHelperNeg) %>% 
    pivot_longer(., cols = starts_with("n"),
                 names_to = "Sign",
                 names_prefix = "n",
                 values_to = "Count") %>% 
    mutate_at(., .vars = "Sign", 
              ~factor(., 
                      levels = c("Negative", "Zero", "Positive"),
                      ordered = TRUE)) %>% 
    arrange(SortHelper)
  # Colorblind friendly palette 
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # Remixing so Zero is grey, negative is blue, and positive is orange
  cbp2 <- c("#56B4E9", "#999999", "#E69F00",  "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  coef_sign_plot <- ggplot(boot_coef_sign_sum, 
                           aes(x = fct_reorder(Variable, SortHelper, .fun = median), 
                               y = Count, fill = Sign)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values=cbp2)+
    labs(y = "Coefficient Estimated Sign Proportion",
         x = "Contrast") +
    coord_flip()
  
  coef_plots <- list(coef_box, coef_zero_bar, coef_sign_plot)
  return(coef_plots)
}

ArrangeBootstrapPlotPair <- function(gplot.lmin, gplot.l1se){
  plot_lmin_with_label <- gplot.lmin + 
    labs(subtitle = "Lambda min (less shrinkage)") + 
    theme(plot.subtitle = element_text(size = 9, hjust = 0),
          axis.title.x = element_text(size = 9),
          axis.text.y = element_text(size = 7.5),
          axis.title.y = element_text(size = 9))
  
  plot_l1se_with_label <- gplot.l1se + 
    labs(subtitle = "Lambda 1 s.e. (more shrinkage)") + 
    theme(plot.subtitle = element_text(size = 9, hjust = 0),
          axis.title.x = element_text(size = 9),
          axis.text.y = element_text(size = 7.5),
          axis.title.y = element_text(size = 9))
  
  paired_plot <- gridExtra::arrangeGrob(
    plot_lmin_with_label, plot_l1se_with_label, nrow = 2)
  
  return(paired_plot)
}

################################################################################
##  Individual plot functions
##
################################################################################

CreateCoefSignPlot <- function(coef.long, sort.plot = TRUE,
                               all.outcomes = FALSE){
  
  stopifnot(sort.plot == TRUE | sort.plot == FALSE)
  if (all.outcomes == TRUE) stopifnot("Outcome" %in% colnames(coef.long))
  
  b_boot_reps <- nrow(coef.long)/length(unique(coef.long$Variable))
  
  
  if(all.outcomes == TRUE) grouping_vars <- c("Variable", "Outcome") else grouping_vars <- c("Variable") 
  
  boot_coef_sign_sum <- coef.long %>% 
    group_by_at(vars(all_of(grouping_vars))) %>% 
    summarise(nZero = sum(Coefficient == 0),
              nPositive = sum(Coefficient > 0),
              nNegative = sum(Coefficient < 0),
              SortHelperPos = sum(Coefficient > 0),
              SortHelperZero = sum(Coefficient == 0),
              SortHelperNeg = sum(Coefficient < 0),
              SortHelper = b_boot_reps*SortHelperPos - SortHelperNeg) %>% 
    pivot_longer(., cols = starts_with("n"),
                 names_to = "Sign",
                 names_prefix = "n",
                 values_to = "Count") %>% 
    mutate_at(., .vars = "Sign", 
              ~factor(., 
                      levels = c("Negative", "Zero", "Positive"),
                      ordered = TRUE)) 
  
  # Remixing so Zero is grey, negative is blue, and positive is orange
  cbp2 <- c("#56B4E9", "#999999", "#E69F00",  "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  if (sort.plot == TRUE) {
    coef_plot_root <- ggplot(boot_coef_sign_sum, 
                             aes(x = fct_reorder(Variable, SortHelper, .fun = median), 
                                 y = Count, fill = Sign))
  } else{
    coef_plot_root <- ggplot(boot_coef_sign_sum, 
                             aes(x = Variable, 
                                 y = Count, fill = Sign))
  }
  
  coef_sign_plot <- coef_plot_root + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values=cbp2)+
    labs(y = "Estimated Sign Proportion",
         x = "Contrast") +
    coord_flip()
  
  if (all.outcomes == TRUE){
    coef_sign_plot <- coef_sign_plot + facet_grid(rows = vars(Outcome)) +
      theme(legend.position="bottom")
  }
  
  return(coef_sign_plot)
}

CreateProbZeroPlot <- function(coef.long, sort.plot = TRUE){
  stopifnot(sort.plot == TRUE | sort.plot == FALSE)
  
  boot_coef_zero_sum <- coef.long %>% 
    group_by(Variable) %>% 
    summarise(ProbZero = sum(Coefficient == 0)/length(Coefficient))
  
  if (sort.plot == TRUE){
    zero_plot_root <- ggplot(data = boot_coef_zero_sum, 
                             aes(x = fct_reorder(Variable, ProbZero, .fun = median,
                                                 .desc = TRUE), 
                                 y = ProbZero))
  } else{
    zero_plot_root <- ggplot(data = boot_coef_zero_sum, 
                             aes(x =Variable, y = ProbZero))
  }
  
  coef_zero_bar <- zero_plot_root +
    geom_col() + 
    labs(y = "Probability Coefficient = 0",
         x = "Variable") +
    coord_flip()
  
  return(coef_zero_bar)
}


CreateCoefBoxplot <- function(coef.long, sort.plot = TRUE,
                              all.outcomes = FALSE){
  
  stopifnot(sort.plot == TRUE | sort.plot == FALSE)
  
  if (sort.plot == TRUE){
    box_plot_root <- ggplot(data = coef.long, 
                            aes(x = fct_reorder(Variable, Coefficient, .fun = median), 
                                y = Coefficient)) 
  } else{
    box_plot_root <- ggplot(data = coef.long, 
                            aes(x = Variable, 
                                y = Coefficient)) 
  }
  
  coef_box <- box_plot_root +
    geom_boxplot() + 
    labs(y = "Log Odds",
         x = "Contrast") + 
    coord_flip()
  
  
  if (all.outcomes == TRUE){
    coef_box <- coef_box + facet_grid(rows = vars(Outcome))
  }
  
  return(coef_box)
}


######################################################
## Helper Functions for Preparing Plot Inputs
##
######################################################

CreateLongDfAllOutcomes <- function(contrast.matrix,
                                    hosp.coef,
                                    icu.coef,
                                    death.coef){
  hosp_contrasts <- t(contrast.matrix) %*% hosp.coef
  icu_contrasts <- t(contrast.matrix) %*% icu.coef 
  death_contrasts <- t(contrast.matrix) %*% death.coef
  
  # Long data frames for plotting
  hosp_cont_long <- CreateCoefLongDfForPlotting(hosp_contrasts)
  icu_cont_long <- CreateCoefLongDfForPlotting(icu_contrasts)
  death_cont_long <- CreateCoefLongDfForPlotting(death_contrasts)
  
  # Order the levels for plot ordering purposes
  long_w_outcomes <- bind_rows(hosp_cont_long %>% mutate(Outcome = "Hospitalization+"),
                               icu_cont_long %>% mutate(Outcome = "ICU+"),
                               death_cont_long %>% mutate(Outcome = "Death"),) %>% 
    mutate_at(., .vars = "Outcome", ~factor(., levels = c("Hospitalization+",
                                                          "ICU+",
                                                          "Death"),
                                            ordered = TRUE))
  
  return(long_w_outcomes)
  
}

CreateCoefLongDfForPlotting <- function(coef.mat, coefs.as.rows = TRUE){
  if (coefs.as.rows == FALSE) coef.mat <- t(coef.mat)
  
  boot_coef_long <- tibble(Variable = rep(rownames(coef.mat), times = ncol(coef.mat)),
                           Coefficient = c(unlist(coef.mat))) %>% 
    mutate(Variable = factor(Variable, levels = rownames(coef.mat))) %>% 
    filter(Variable != "(Intercept)")
  
  return(boot_coef_long)
}

#' Helper function for calculating the weights to use when averaging over binary 
#' interaction terms. 
#' 
InteractionContrastWeights <- function(var1.name, var2.name, in.dat = iph){
  #Want the conditional probability given the first variable is true
  pt <- in.dat %>% select(all_of(c(var1.name, var2.name))) %>% 
    table %>% 
    prop.table(., margin = 1) %>% 
    .[2,]
  
  prob_weights <- c(-pt[1], pt[1], -pt[2], pt[2])
  names(prob_weights) <- c(paste0(var1.name, "FALSE:", var2.name, "FALSE"),
                           paste0(var1.name, "FALSE:", var2.name, "TRUE"),
                           paste0(var1.name, "TRUE:", var2.name, "FALSE"),
                           paste0(var1.name, "TRUE:", var2.name, "TRUE"))
  return(prob_weights)
}
