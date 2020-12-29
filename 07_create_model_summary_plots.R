################################################################################
##  Read in relevant data
##
################################################################################

source("./01_Regression/lasso.R")
source("./01_Regression/glmnet_bootstrap.R")
source("./02_Plots/construct_contrasts.R")
source("./02_Plots/bootstrap_plots.R")

if(exists("hosp_formula") == FALSE){
  hosp_formula <- readRDS("./Work_Products/Model_Formulas/hosp_formula.rds")
  icu_formula <- readRDS("./Work_Products/Model_Formulas/icu_formula.rds")
  death_formula <- readRDS("./Work_Products/Model_Formulas/death_formula.rds")
}

if(exists("coef_mat_list") == FALSE){
  coef_mat_list <- readRDS("./Work_Products/Model_Coefficients/coef_mat_list.rds")
  
}
hosp_coef <- coef_mat_list$hosp_min_coef_mat
icu_coef <- coef_mat_list$icu_min_coef_mat
death_coef <- coef_mat_list$death_min_coef_mat

################################################################################
## Magic constants  
##
################################################################################

m_imputed <- 30
b_bootstrap_reps <- ncol(hosp_coef) / m_imputed

################################################################################
## Settings  
##
################################################################################
target_contrast_dir <- "./Work_Products/Contrast_Matrices/"

state_contrast_dir <- "./Work_Products/Contrast_Matrices/Geographic/States"

country_contrast_dir <- "./Work_Products/Contrast_Matrices/Geographic/Countries"

# Subset the contrasts to only include the top effects (in terms of absolute value of the mean estimate)
top_n_effects_flag <- FALSE
top_n_cutoff <- 10

sort_coef_plots <- FALSE

save_plots_flag <- FALSE

average_over_imputations_first <- FALSE

################################################################################
##  Read in contrasts and calculate estimates
##
################################################################################

contrast_matrix <- ConstructSingleMatrixFromContrastDir(target_contrast_dir)

country_contrast_matrix <- ConstructSingleMatrixFromContrastDir(country_contrast_dir)

state_contrast_matrix <- ConstructSingleMatrixFromContrastDir(state_contrast_dir)



hosp_contrasts <- t(contrast_matrix) %*% hosp_coef
icu_contrasts <- t(contrast_matrix) %*% icu_coef 
death_contrasts <- t(contrast_matrix) %*% death_coef

if(average_over_imputations_first == TRUE){
  hosp_contrasts_avg <- matrix(nrow = nrow(hosp_contrasts),
                               ncol = ncol(hosp_contrasts)/m_imputed)
  icu_contrasts_avg <- matrix(nrow = nrow(hosp_contrasts),
                              ncol = ncol(hosp_contrasts)/m_imputed)
  death_contrasts_avg <- matrix(nrow = nrow(hosp_contrasts),
                                ncol = ncol(hosp_contrasts)/m_imputed)
  
  for(i in 1:b_bootstrap_reps){
    cur_col_indices <- 1:m_imputed + m_imputed*(i-1)
    print(cur_col_indices[1])
    hosp_contrasts_avg[,i] <- rowMeans(hosp_contrasts[,cur_col_indices])
    icu_contrasts_avg[,i] <- rowMeans(icu_contrasts[,cur_col_indices])
    death_contrasts_avg[,i] <- rowMeans(death_contrasts[,cur_col_indices])
  }
  rownames(hosp_contrasts_avg) <- rownames(hosp_contrasts)
  rownames(icu_contrasts_avg) <- rownames(icu_contrasts)
  rownames(death_contrasts_avg) <- rownames(death_contrasts)
  
  hosp_contrasts <- hosp_contrasts_avg
  icu_contrasts <- icu_contrasts_avg
  death_contrasts <- death_contrasts_avg
}

# Long data frames for plotting
hosp_cont_long <- CreateCoefLongDfForPlotting(hosp_contrasts)
icu_cont_long <- CreateCoefLongDfForPlotting(icu_contrasts)
death_cont_long <- CreateCoefLongDfForPlotting(death_contrasts)

primary_comp_long_w_outcomes <- CreateLongDfAllOutcomes(contrast.matrix = contrast_matrix,
                                                        hosp.coef = hosp_coef,
                                                        icu.coef = icu_coef,
                                                        death.coef = death_coef)

country_comp_long_w_outcomes <- CreateLongDfAllOutcomes(contrast.matrix = country_contrast_matrix,
                                                        hosp.coef = hosp_coef,
                                                        icu.coef = icu_coef,
                                                        death.coef = death_coef)

state_comp_long_w_outcomes <- CreateLongDfAllOutcomes(contrast.matrix = state_contrast_matrix,
                                                      hosp.coef = hosp_coef,
                                                      icu.coef = icu_coef,
                                                      death.coef = death_coef)

if (top_n_effects_flag == TRUE) {
  hosp_cont_long_vars <- hosp_cont_long %>% group_by(Variable) %>% 
    summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
    arrange(desc(AbsMeanCoef)) %>% 
    slice_head(., n = top_n_cutoff) %>% 
    pull(Variable)
  
  icu_cont_long_vars <- icu_cont_long %>% group_by(Variable) %>% 
    summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
    arrange(desc(AbsMeanCoef)) %>% 
    slice_head(., n = top_n_cutoff) %>% 
    pull(Variable)
  
  death_cont_long_vars <- death_cont_long %>% group_by(Variable) %>% 
    summarise(AbsMeanCoef = abs(mean(Coefficient))) %>% 
    arrange(desc(AbsMeanCoef)) %>% 
    slice_head(., n = top_n_cutoff) %>% 
    pull(Variable)
  
  hosp_cont_long <- hosp_cont_long %>% filter(Variable %in% hosp_cont_long_vars)
  icu_cont_long <- icu_cont_long %>% filter(Variable %in% icu_cont_long_vars)
  death_cont_long <- death_cont_long %>% filter(Variable %in% death_cont_long_vars)
}

################################################################################
##  Create summary dfs for referencing
##  
################################################################################

#' Create a data frame with summary statistics across bootstrap and imputations
#' Expects 
CreateResultsSummaryDF <- function(long.coef.df, 
                                   single.outcome = FALSE){
  
  grouping_vars <- c("Variable", "Outcome")
  if (single.outcome == TRUE)  grouping_vars <- "Variable"
  
  results_summary <- long.coef.df %>% 
    group_by_at(vars(all_of(grouping_vars))) %>% 
    summarise(AbsMeanCoef = abs(mean(Coefficient)),
              MeanCoef = mean(Coefficient),
              OR = exp(MeanCoef),
              ProbZero = mean(Coefficient == 0),
              ProbPos = mean(Coefficient > 0),
              ProbNeg = mean(Coefficient < 0)) %>% 
    arrange(desc(AbsMeanCoef))
  
  return(results_summary)
}

main_sum <- CreateResultsSummaryDF(long.coef.df = primary_comp_long_w_outcomes,
                                   single.outcome = FALSE)

hosp_sum <- hosp_cont_long %>% group_by(Variable) %>% 
  summarise(AbsMeanCoef = abs(mean(Coefficient)),
            MeanCoef = mean(Coefficient),
            OR = exp(MeanCoef)) %>% 
  arrange(desc(AbsMeanCoef))

icu_sum <- icu_cont_long %>% group_by(Variable) %>% 
  summarise(AbsMeanCoef = abs(mean(Coefficient)),
            MeanCoef = mean(Coefficient),
            OR = exp(MeanCoef)) %>% 
  arrange(desc(AbsMeanCoef)) 

death_sum <- death_cont_long %>% group_by(Variable) %>% 
  summarise(AbsMeanCoef = abs(mean(Coefficient)),
            MeanCoef = mean(Coefficient),
            OR = exp(MeanCoef)) %>% 
  arrange(desc(AbsMeanCoef))

################################################################################
##  Create the Plots
##
################################################################################

sign_plots_top <- CreateCoefSummaryPlots(contrast.matrix = contrast_matrix,
                                         hosp.coef = hosp_coef,
                                         icu.coef = icu_coef,
                                         death.coef = death_coef,
                                         plot.type = "sign",
                                         top.n.effects.flag = TRUE,
                                         top.n.cutoff = 10,
                                         sort.coef.plots = TRUE)

box_plots_top <- CreateCoefSummaryPlots(contrast.matrix = contrast_matrix,
                                        hosp.coef = hosp_coef,
                                        icu.coef = icu_coef,
                                        death.coef = death_coef,
                                        plot.type = "box",
                                        top.n.effects.flag = TRUE,
                                        top.n.cutoff = 10,
                                        sort.coef.plots = TRUE)

sign_plots_all <- CreateCoefSummaryPlots(contrast.matrix = contrast_matrix,
                                         hosp.coef = hosp_coef,
                                         icu.coef = icu_coef,
                                         death.coef = death_coef,
                                         plot.type = "sign",
                                         top.n.effects.flag = FALSE,
                                         sort.coef.plots = FALSE)

box_plots_all <- CreateCoefSummaryPlots(contrast.matrix = contrast_matrix,
                                        hosp.coef = hosp_coef,
                                        icu.coef = icu_coef,
                                        death.coef = death_coef,
                                        plot.type = "box",
                                        top.n.effects.flag = FALSE,
                                        sort.coef.plots = FALSE)

sign_plots_country <- CreateCoefSummaryPlots(contrast.matrix = country_contrast_matrix,
                                             hosp.coef = hosp_coef,
                                             icu.coef = icu_coef,
                                             death.coef = death_coef,
                                             plot.type = "sign",
                                             top.n.effects.flag = FALSE,
                                             sort.coef.plots = FALSE,
                                             var.axis.label = "Country")

box_plots_country <- CreateCoefSummaryPlots(contrast.matrix = country_contrast_matrix,
                                            hosp.coef = hosp_coef,
                                            icu.coef = icu_coef,
                                            death.coef = death_coef,
                                            plot.type = "box",
                                            top.n.effects.flag = FALSE,
                                            sort.coef.plots = FALSE,
                                            var.axis.label = "Country")

sign_plots_state <- CreateCoefSummaryPlots(contrast.matrix = state_contrast_matrix,
                                           hosp.coef = hosp_coef,
                                           icu.coef = icu_coef,
                                           death.coef = death_coef,
                                           plot.type = "sign",
                                           top.n.effects.flag = FALSE,
                                           sort.coef.plots = FALSE,
                                           var.axis.label = "U.S. Region")

box_plots_state <- CreateCoefSummaryPlots(contrast.matrix = state_contrast_matrix,
                                          hosp.coef = hosp_coef,
                                          icu.coef = icu_coef,
                                          death.coef = death_coef,
                                          plot.type = "box",
                                          top.n.effects.flag = FALSE,
                                          sort.coef.plots = FALSE,
                                          var.axis.label = "U.S. Region")
