library(tidyverse)
library(future)
library(glmnet)
library(glmnetUtils)
library(tictoc)

source("./01_Regression/lasso.R")
source("./01_Regression/glmnet_bootstrap.R")


iph_train_imputed <- readRDS("./data/Imputed Datasets/data_set(10-20-2020)/train_mice_std.rds")

hosp_formula <- readRDS("./Work_Products/Model_Formulas/hosp_formula.rds")
icu_formula <- readRDS("./Work_Products/Model_Formulas/icu_formula.rds")
death_formula <- readRDS("./Work_Products/Model_Formulas/death_formula.rds")

lambda_seq_list <- readRDS("./Work_Products/Lambda_Vals/lambda_seq_list.rds")

formula_list <- list(hosp = hosp_formula,
                     icu = icu_formula,
                     death = death_formula)

outcome.names <- c("HospOrDeath", "IcuVentDeath", "death")
names(formula_list) <- outcome.names

plan(multisession)

######################################################
## Settings
######################################################
save_flag <- TRUE

iters_to_run <- 10

# https://nclottery.com/Pick4
# Wednesday Oct 13 daytime draw

sim_seed <- 2000

start_time <- lubridate::now()

warning("Alpha values not currently used except for naming the file. Function will use Lasso penalty at the moment")
hosp_alpha <- 1
icu_alpha <- 1
death_alpha <- 1

time_string <- as.character(start_time) %>% str_replace_all(., pattern = "\\:", replacement = "_")
dir_string <- "./data/Fitted Models/data_set(10-20-2020)/"

hosp_string <- paste0("_hosp_", hosp_alpha, "alpha_",  
                      iters_to_run, "iter_", 
                      sim_seed, "seed")

icu_string <- paste0("_icu_", icu_alpha, "alpha_",  
                     iters_to_run, "iter_", 
                     sim_seed, "seed")

death_string <- paste0("_death_", death_alpha, "alpha_",  
                       iters_to_run, "iter_", 
                       sim_seed, "seed")

######################################################
## Run - Hospital
######################################################
set.seed(sim_seed)

tic()
all_mods_parallel <- BootMiceSimultaneous(
  mice.data = iph_train_imputed, 
  n.iter = iters_to_run, 
  formula.list = formula_list,
  outcome.names = c("HospOrDeath", "IcuVentDeath", "death"),
  lambda.seq.list = lambda_seq_list,
  parallel = TRUE)

toc()


reordered_mods <- TransformBootstrapListToModelList(all_mods_parallel)

hosp_mods <- reordered_mods[[1]]
icu_mods <- reordered_mods[[2]]
death_mods <- reordered_mods[[3]]
######################################################
## Save
######################################################


if(save_flag == TRUE){
  saveRDS(hosp_mods, paste0(dir_string, time_string, hosp_string, ".rds"))
  saveRDS(icu_mods, paste0(dir_string, time_string, icu_string,  ".rds"))
  saveRDS(death_mods, paste0(dir_string, time_string, death_string, ".rds"))
} 

