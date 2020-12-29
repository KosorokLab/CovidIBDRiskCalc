#' Wrapper around mice 
#' 
#' @param in.data original data set
#' @param vars.to.include
#' @param exclude.pred.vars
#' @param impute.then.transform logical flag indicating whether
ImputeDatasetDefaultSettings <- function(in.data, vars.to.include, exclude.pred.vars, impute.then.transform = TRUE,
                                         ...){
  data_to_impute <- in.data %>% select(all_of(vars.to.include))
  
  mice_pred_vars <- mice::quickpred(data_to_impute, 
                                    mincor = 0.025,
                                    method = "spearman",
                              exclude = exclude.pred.vars)
  
  # Split, impute, train to prevent data leakage 
  # Seed is https://nclottery.com/Pick4 Daytime drawing Monday September 28 2020
  imputed_data <- mice::mice(data = data_to_impute,
                             predictorMatrix = mice_pred_vars,
                             m = 30,
                             maxit = 10,
                             seed = 7210,
                             ...)
  
  
  if (impute.then.transform == TRUE) {
    mice_long <- mice::complete(imputed_data, "long", include = TRUE)
    mice_long_w_transformed <- .ImputeTransformedVariables(mice_long)
    imputed_data <- mice::as.mids(long = mice_long_w_transformed)
  }
  
  return(imputed_data)
}

#' Helper function to calculate the values for transformed variables e.g. BMI
#' @param imputed_long a data frame with all of the imputed datasets stacked on top of each other
#' @return dataframe with the derived variables appended
.ImputeTransformedVariables <- function(imputed_long){
  
  imputed_long <- .DeriveDrugVars(imputed_long)
  
  imputed_long$bmiAbv60Rm <- with(imputed_long,  weightBMIAbv60Rm/((heightCombinedCm/100)^2)) 
  imputed_long$HospOrDeath <- with(imputed_long,  HospIcuDthOrd != "None")
  imputed_long$IcuVentDeath <- with(imputed_long,  HospIcuDthOrd %in% c("ICU", "Ventilated", "Died"))
  
  imputed_long$hosp <- with(imputed_long,  HospIcuDthOrd == "Hospitalized")
  imputed_long$IcuOrVent <- with(imputed_long,  HospIcuDthOrd == "ICU")
  imputed_long$death <- with(imputed_long,  HospIcuDthOrd == "Died")
  
  imputed_long %<>% mutate_at(.vars = paste0("race___", 0:4), factor)
  
  return(imputed_long)
}

#' Write the imputed data out to an RDS file. Essentially a wrapper around 
#' saveRDS that also creates the relevant directory if it doesn't currently exist. 
WriteImputedDataList <- function(imputed.data, 
                                 out.name, 
                                 iph.path, 
                                 analytic.filename) {
  
  dataset_name <- str_extract(analytic.filename, ".*(?=\\.csv)")
  out_dir <- paste0(iph.path, "Imputed Datasets/", dataset_name)
  
  if (dir.exists(out_dir) == FALSE) dir.create(out_dir, recursive = TRUE)
  
  out_file_name <- paste(out_dir, out.name, sep = "/")
  saveRDS(imputed.data, out_file_name)
}