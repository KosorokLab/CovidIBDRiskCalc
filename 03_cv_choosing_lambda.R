library(doParallel)

cl <- makePSOCKcluster(12)
registerDoParallel(cl)

library(tidyverse)
library(glmnet)
library(glmnetUtils)

source("./01_Regression/lasso.R")
source("./01_Regression/glmnet_bootstrap.R")

################################################################################
##
## Data
##
################################################################################


if (!exists("iph_train_imputed")) iph_train_imputed <- readRDS("data/Imputed Datasets/data_set(10-20-2020)/train_mice_std.rds")

imp_data_list <- mice::complete(iph_train_imputed, "all")

hosp_formula <- readRDS("./Work_Products/Model_Formulas/hosp_formula.rds")
icu_formula <- readRDS("./Work_Products/Model_Formulas//icu_formula.rds")
death_formula <- readRDS("./Work_Products/Model_Formulas/death_formula.rds")

################################################################################
##
## Settings
##
################################################################################

save_cv_model_fits <- FALSE
save_cv_model_dir <- NULL

################################################################################
##
## Fit Models
##
################################################################################


hosp_cv_deviance_list <- map(imp_data_list, 
                              ~CvGlmnet(in.formula = hosp_formula, 
                                         in.data = .,
                                         na.action = na.omit,
                                         family = "binomial", 
                                         nfolds = 10,
                                         type.measure = "deviance",
                                        parallel=TRUE))


icu_cv_deviance_list <- map(imp_data_list, 
                             ~CvGlmnet(in.formula = icu_formula, 
                                       in.data = .,
                                       na.action = na.omit,
                                       family = "binomial", 
                                       nfolds = 10,
                                       type.measure = "deviance",
                                       parallel=TRUE))

death_cv_deviance_list <- map(imp_data_list, 
                              ~CvGlmnet(in.formula = death_formula, 
                                        in.data = .,
                                        na.action = na.omit,
                                        family = "binomial", 
                                        nfolds = 10,
                                        type.measure = "deviance",
                                        parallel=TRUE))

stopCluster(cl)


if (save_cv_model_fits == TRUE){
  saveRDS(hosp_cv_deviance_list, "hosp_cv_list.rds")
  saveRDS(icu_cv_deviance_list, "icu_cv_list.rds")
  saveRDS(death_cv_deviance_list, "death_cv_list.rds")
  
}

################################################################################
##
## Functions for accessing the deviance ratio and lambda values
##
################################################################################


GetLambdaAndDeviancesFromCVFit <- function(in.cv.fit){
  fit_results <- cbind("value" = c("Min", "1se"),
                       "lambda" = c(in.cv.fit$lambda.min,
                                    in.cv.fit$lambda.1se),
                       "dev.ratio" = c(in.cv.fit$glmnet.fit$dev.ratio[which.min(in.cv.fit$cvm)],
                                       in.cv.fit$glmnet.fit$dev.ratio[which(in.cv.fit$lambda == in.cv.fit$lambda.1se)]))
  return(fit_results)
}

GetDeviancesFromCvList <- function(cv.list){
  GetLambdaAndDeviancesFromCVFit <- function(in.cv.fit){
    fit_results <- cbind("value" = c("Min", "1se"),
                         "lambda" = c(in.cv.fit$lambda.min,
                                      in.cv.fit$lambda.1se),
                         "dev.ratio" = c(in.cv.fit$glmnet.fit$dev.ratio[which.min(in.cv.fit$cvm)],
                                         in.cv.fit$glmnet.fit$dev.ratio[which(in.cv.fit$lambda == in.cv.fit$lambda.1se)]))
    return(fit_results)
  }
  fit_results_cv <- map(cv.list, GetLambdaAndDeviancesFromCVFit)
  lambda_types <- map(fit_results_cv, ~(.[,1])) %>% unlist
  lambda_vals <- map(fit_results_cv, ~(.[,2])) %>% unlist
  dev_ratios <- map(fit_results_cv, ~(.[,3])) %>% unlist
  
  fit_df <- tibble(LambdaType = lambda_types,
                   Lambda = lambda_vals,
                   DevRatio = dev_ratios) %>% 
    mutate_at(., .vars = c("Lambda", "DevRatio"), as.numeric) %>% 
    arrange(desc(DevRatio))
  return(fit_df)
}

################################################################################
## Get Lambda values
## 
##
################################################################################

hosp_dev <- GetDeviancesFromCvList(hosp_cv_deviance_list)
icu_dev <- GetDeviancesFromCvList(icu_cv_deviance_list)
death_dev <- GetDeviancesFromCvList(death_cv_deviance_list)

hosp_summary <- hosp_dev %>% 
  group_by(LambdaType) %>% 
  summarise(MeanLambda = mean(Lambda), 
            MeanDevRatio = mean(DevRatio))

icu_summary <- icu_dev %>% 
  group_by(LambdaType) %>% 
  summarise(MeanLambda = mean(Lambda), 
            MeanDevRatio = mean(DevRatio))

 death_summary <- death_dev %>% 
  group_by(LambdaType) %>% 
  summarise(MeanLambda = mean(Lambda), 
            MeanDevRatio = mean(DevRatio))

 lambda_summary <- bind_rows(hosp_summary,
                             icu_summary,
                             death_summary)  %>% 
   select(-MeanDevRatio)

 lambda_val_list <- lambda_summary$MeanLambda
 names(lambda_val_list) <- c("hosp_lambda_1se",
                             "hosp_lambda_min",
                             "icu_lambda_1se",
                             "icu_lambda_min",
                             "death_lambda_1se",
                             "death_lambda_min")
 
 lambda_val_list <- as.list(lambda_val_list)
 
 saveRDS(lambda_val_list, "./Work_Products/Lambda_Vals/lambda_vals.rds") 
 