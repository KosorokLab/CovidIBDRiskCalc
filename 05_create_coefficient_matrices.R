source("./01_Regression/glmnet_bootstrap.R")
library(tidyverse)
library(glmnet)

# Date is based off of the entire data set, models only fit to training data which is a 
# random subsample of the data through September 15th
dir_string <- "./data/Fitted Models/data_set(10-20-2020)/"
hosp_string <- "2020-11-30 14_27_55_hosp_1alpha_1000iter_2000seed.rds"
icu_string <-  "2020-11-30 14_27_55_icu_1alpha_1000iter_2000seed.rds"
death_string <- "2020-11-30 14_27_55_death_1alpha_1000iter_2000seed.rds"

hosp_glmnet_fit <- readRDS(paste0(dir_string, hosp_string))
icu_glmnet_fit <- readRDS(paste0(dir_string, icu_string))
death_glmnet_fit <- readRDS(paste0(dir_string, death_string))

######################################################
## Settings
######################################################
cv_fit_flag <- FALSE
# Indicates whether the fits are a list of lists or a single list of fits
# A list of list results from having each bootstrap sample of indices re-used across all imputed data sets
nested_bootstrap_flag <- TRUE

lambda_val_list <- readRDS("./Work_Products/Lambda_Vals/lambda_vals.rds") 

hosp_lambda_1se <- lambda_val_list$hosp_lambda_1se
hosp_lambda_min <- lambda_val_list$hosp_lambda_min 

icu_lambda_1se <- lambda_val_list$icu_lambda_1se
icu_lambda_min <- lambda_val_list$icu_lambda_min

death_lambda_1se <- lambda_val_list$death_lambda_1se
death_lambda_min <- lambda_val_list$death_lambda_min


######################################################
## Settings
######################################################

if (nested_bootstrap_flag == TRUE){
  hosp_glmnet_fit <- unlist(hosp_glmnet_fit, recursive = FALSE)
  icu_glmnet_fit <- unlist(icu_glmnet_fit, recursive = FALSE)
  death_glmnet_fit <- unlist(death_glmnet_fit, recursive = FALSE)
}

if (cv_fit_flag == TRUE) {
  hosp_min_coef_mat <- CreateBootCoefMatrix(hosp_glmnet_fit, "lambda.min")
  hosp_1se_coef_mat <- CreateBootCoefMatrix(hosp_glmnet_fit, "lambda.1se")
  
  icu_min_coef_mat <- CreateBootCoefMatrix(icu_glmnet_fit, "lambda.min")
  icu_1se_coef_mat <- CreateBootCoefMatrix(icu_glmnet_fit, "lambda.1se")
  
  death_min_coef_mat <- CreateBootCoefMatrix(death_glmnet_fit, "lambda.min")
  death_1se_coef_mat <- CreateBootCoefMatrix(death_glmnet_fit, "lambda.1se")
  
}

if (cv_fit_flag == FALSE) {
  hosp_min_coef_mat <- CreateBootCoefMatrix(hosp_glmnet_fit, hosp_lambda_min)
  hosp_1se_coef_mat <- CreateBootCoefMatrix(hosp_glmnet_fit, hosp_lambda_1se)
  
  icu_min_coef_mat <- CreateBootCoefMatrix(icu_glmnet_fit, icu_lambda_min)
  icu_1se_coef_mat <- CreateBootCoefMatrix(icu_glmnet_fit, icu_lambda_1se)
  
  death_min_coef_mat <- CreateBootCoefMatrix(death_glmnet_fit, death_lambda_min)
  death_1se_coef_mat <- CreateBootCoefMatrix(death_glmnet_fit, death_lambda_1se)
  
}

coef_mat_list <- list(source_files  = c("hosp" = hosp_string, "icu" = icu_string, "death" = death_string),
                      hosp_min_coef_mat = hosp_min_coef_mat,
                      hosp_1se_coef_mat = hosp_1se_coef_mat,
                      icu_min_coef_mat = icu_min_coef_mat,
                      icu_1se_coef_mat = icu_1se_coef_mat,
                      death_min_coef_mat = death_min_coef_mat,
                      death_1se_coef_mat = death_1se_coef_mat)

saveRDS(coef_mat_list, file = "./Work_Products/Model_Coefficients/coef_mat_list.rds")
