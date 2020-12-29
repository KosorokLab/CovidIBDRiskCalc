if (!exists("iph"))  source("01_determine_variables_for_modeling.R") #This will load the data as well
source("./00_Data_Processing/impute_data.R")
source("./01_Regression/lasso.R")

################################################################################
##
## Settings
##
################################################################################

run_train_imputation <- TRUE
run_test_imputation <- TRUE

save_imputations_flag <- TRUE

impute_then_transform_root_vars <- c("preddose_daily",
                                     "mcpdose",
                                     "azadose",
                                     "heightCombinedCm",
                                     "weightBMIAbv60Rm",
                                     "HospIcuDthOrd")

vars_for_derivation_only <- c("immod___1",
                              "immod___4",
                              "immod___5",
                              "numeric_studyid",
                              "comorbid___12",
                              "comorbid___13")

impute_then_transform_child_vars <- c("preddose_dailyByWt",
                                     "mcpdoseByWt",
                                     "azadoseByWt",
                                     "bmiAbv60Rm",
                                     "HospOrDeath",
                                     "IcuVentDeath")

# These country indicators are used to prevent colinearity between the 
# Other countries indicator in the state variable and the US level of the country variable
# These are indicators for countries with >= 30 observations excluding the US
country_indicators <- iph %>% select(starts_with("countryInd")) %>% colnames

vars_to_keep <- c(candidate_vars[!(candidate_vars %in%impute_then_transform_child_vars)],
                  impute_then_transform_root_vars,
                  vars_for_derivation_only,
                  country_indicators)

vars_exclude_from_imp_pred <- c("countryGrouped",
                                "Other_CS_flg",
                                vars_for_derivation_only,
                                "race___2", 
                                "race___4", 
                                impute_then_transform_child_vars)

################################################################################
##
## Training set imputation
##
################################################################################

if (run_train_imputation == TRUE) {
  #Seed is set within the MICE function
  iph_train_imputed <- ImputeDatasetDefaultSettings(
    in.data = iph_train,
    vars.to.include = vars_to_keep,
    exclude.pred.vars = vars_exclude_from_imp_pred,
    impute.then.transform = TRUE)
  
  if (save_imputations_flag == TRUE){
  
    WriteImputedDataList(imputed.data = iph_train_imputed, 
                         out.name = "train_mice_std.rds", 
                         iph.path = data_path, 
                         analytic.filename = analytic_filename)
  }
}


################################################################################
##
## Test set imputation
##
################################################################################

if (run_test_imputation == TRUE){
  iph_test_imputed <- ImputeDatasetDefaultSettings(
    in.data = iph_test,
    vars.to.include = vars_to_keep,
    exclude.pred.vars = c(vars_exclude_from_imp_pred, "countryIndIran"),
    impute.then.transform = TRUE)
  
  if (save_imputations_flag == TRUE){
    
    WriteImputedDataList(iph_test_imputed, 
                         "test_mice_std.rds", 
                         data_path, 
                         analytic_filename)
  }
}