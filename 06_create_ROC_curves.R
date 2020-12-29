################################################################################
##  Read in relevant data
##
################################################################################
if(exists("iph") == FALSE){
  source("01_determine_variables_for_modeling.R") #This will load the data as well
}

if(exists("BootPredictMatrix") == FALSE){
  source("./01_Regression/lasso.R")
  source("./01_Regression/glmnet_bootstrap.R")
}

if(exists("iph_train_imputed") == FALSE){
  iph_train_imputed <- readRDS(paste0(data_path, "/Imputed Datasets/data_set(10-20-2020)/train_mice_std.rds"))
  iph_test_imputed <- readRDS(paste0(data_path, "/Imputed Datasets/data_set(10-20-2020)/test_mice_std.rds"))
}

if(exists("hosp_formula") == FALSE){
  hosp_formula <- readRDS("./Work_Products/Model_Formulas/hosp_formula.rds")
  icu_formula <- readRDS("./Work_Products/Model_Formulas//icu_formula.rds")
  death_formula <- readRDS("./Work_Products/Model_Formulas/death_formula.rds")
}

if(exists("coef_mat_list") == FALSE){
  coef_mat_list <- readRDS("./Work_Products/Model_Coefficients/coef_mat_list.rds")
}

iph_train_imputed_list <- mice::complete(iph_train_imputed, "all")

iph_test_imputed_list <- mice::complete(iph_test_imputed, "all")


################################################################################
##  Settings
##
################################################################################
save_roc_curves_and_summaries <- TRUE

################################################################################
##  ROC Package Wrappers
##
################################################################################


proc_wrapper <- function(pred.mat, outcome.name, test.data,
                         model.name) {
  roc_pred_dat <- data.frame(apply(pred.mat, 1, mean), 
                             test.data %>% pull(outcome.name) %>% as.numeric)
  colnames(roc_pred_dat) <- c("Predicted", "Label")
  
  pROC_obj <- pROC::roc(
    response = roc_pred_dat$Label,
    predictor = roc_pred_dat$Predicted,
    smoothed = FALSE,
    # arguments for ci
    ci=TRUE, ci.alpha=0.9, stratified=FALSE,
    # arguments for plot
    plot=TRUE, grid=TRUE,
    legacy.axes=TRUE, 
    print.auc=TRUE, show.thres=TRUE, 
    main = paste(model.name, outcome.name, sep = " - "))
  
  print(pROC_obj)
  return(pROC_obj)
}

ExtractCIFromProcList <- function(proc.list){
  imputed_aucs <- map_dbl(proc.list, ~.$auc)
  
  # Element one is lower 95% limit, Two is mean, Three is Upper 95% limit
  imputed_lower_lims <- map_dbl(proc.list, ~.$ci[1])
  imputed_upper_lims <- map_dbl(proc.list, ~.$ci[3])
  
  ci_summary <- c("Mean AUC" = mean(imputed_aucs),
                  "Lower 95% CI" = mean(imputed_lower_lims),
                  "Upper 95% CI" = mean(imputed_upper_lims))
  return(ci_summary)
}

################################################################################
##  Create Predicted probabilities
##
################################################################################

if(exists("hosp_test_imputed_predicted") == FALSE){
  
  hosp_test_imputed_predicted <- map(iph_test_imputed_list, 
                                     ~BootPredictMatrix(coef_mat_list$hosp_min_coef_mat,
                                                        hosp_formula,
                                                        .))
  
  icu_test_imputed_predicted <- map(iph_test_imputed_list, 
                                    ~BootPredictMatrix(coef_mat_list$icu_min_coef_mat,
                                                       icu_formula,
                                                       .))
  
  death_test_imputed_predicted <- map(iph_test_imputed_list, 
                                      ~BootPredictMatrix(coef_mat_list$icu_min_coef_mat,
                                                         death_formula,
                                                         .))
}



################################################################################
##  Create ROC Plots - imputed data
##
################################################################################



hosp_test_imputed_proc_curvs <- map2(.x = hosp_test_imputed_predicted,
                                     .y = iph_test_imputed_list,
                                     ~proc_wrapper(
                                       pred.mat = .x, 
                                       outcome.name = "HospOrDeath", 
                                       test.data = .y,
                                       model.name = "Imputed Test Data "))

icu_test_imputed_proc_curvs <- map2(.x = icu_test_imputed_predicted,
                                    .y = iph_test_imputed_list,
                                    ~proc_wrapper(
                                      pred.mat = .x, 
                                      outcome.name = "IcuVentDeath", 
                                      test.data = .y,
                                      model.name = "Imputed Test Data "))

death_test_imputed_proc_curvs <- map2(.x = death_test_imputed_predicted,
                                      .y = iph_test_imputed_list,
                                      ~proc_wrapper(
                                        pred.mat = .x, 
                                        outcome.name = "death", 
                                        test.data = .y,
                                        model.name = "Imputed Test Data "))

hosp_auc_summary <- ExtractCIFromProcList(hosp_test_imputed_proc_curvs)

icu_auc_summary <- ExtractCIFromProcList(icu_test_imputed_proc_curvs)

death_auc_summary <- ExtractCIFromProcList(death_test_imputed_proc_curvs)

auc_sums <- bind_rows(hosp_auc_summary, icu_auc_summary, death_auc_summary) %>% 
  mutate(Outcome = c("Hosp+", "ICU+", "Death")) %>% 
  select(Outcome, everything())

################################################################################
##  Write out files
##
################################################################################

if (save_roc_curves_and_summaries == TRUE){

  imp_roc_curves_list <- list( hosp_test_imputed_rocs = hosp_test_imputed_proc_curvs,
                               icu_test_imputed_rocs = icu_test_imputed_proc_curvs,
                               death_test_imputed_rocs = death_test_imputed_proc_curvs)
  
  auc_sum_list <- list(hosp_auc_summary = hosp_auc_summary,
                       icu_auc_summary = icu_auc_summary,
                       death_auc_summary = death_auc_summary,
                       auc_sums =auc_sums)
  
  saveRDS(imp_roc_curves_list, "./Work_Products/ROC_Curves/imp_roc_curves_list.rds")
  saveRDS(auc_sum_list, "./Work_Products/ROC_Curves/auc_sum_list.rds")
  
}