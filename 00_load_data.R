################################################################################
##
## This script reads in the data. By sourcing this script the full data, 
## training set, and test set should be available in the global environment. 
##
################################################################################

library(tidyverse)
library(magrittr)

source("./00_Data_Processing/import_iph_data.R")

################################################################################
##
## Settings
##
################################################################################

analytic_filename <- "data_set(10-20-2020).csv"

# If TRUE, then a list of numeric study ids will be read in and used for the test-train
# split instead of a genearting a test-train split using random number generation.
use_existing_train_split <- TRUE

# If TRUE, create a test-train split using random number generation and save the 
# study ids to a csv
save_random_train_split <- FALSE

# Split, impute, train to prevent data leakage 
# Seed is https://nclottery.com/Pick4 Evening drawing Sunday September 27 2020
test_train_split_seed <- 1420

train_proportion <- .85

################################################################################
##
## Read in the data
##
################################################################################

data_path <- "./data/"

colspec_data <- read_csv("./00_Data_Processing/column_spec.csv",
                         col_types = "cccc") 

state_region_map_path <- "./00_Data_Processing/state_region_map.csv"

iph <- ReadIPHData(path = paste0(data_path, analytic_filename), 
                   col_spec_df = colspec_data,
                   clean = TRUE,
                   state.map.path = state_region_map_path,
                   hardcode.dir = "./data/",
                   country.inc.dir = data_path)


################################################################################
##
## Check that nothing has gone terribly wrong
## This will only catch extreme issues
##
################################################################################

iph_default_read <- ReadIPHData(path = paste0(data_path, analytic_filename), 
                   col_spec_df = colspec_data,
                   clean = FALSE)

cols_to_check <- c("numeric_studyid", "sex", "physician")

stopifnot(all_equal(iph %>% select(all_of(cols_to_check)),
          iph_default_read %>% select(all_of(cols_to_check))))


################################################################################
##
## Create the train and test data sets
##
################################################################################

if (use_existing_train_split == FALSE){
  set.seed(test_train_split_seed)
  
  iph_train <- iph %>% 
    group_by(HospIcuVentDthOrd) %>% 
    sample_frac(., size = train_proportion, replace = FALSE) %>% 
    ungroup
  
  if (save_random_train_split == TRUE){
    iph_train %>% select(numeric_studyid) %>% 
      write_csv(., paste0(data_path, "training_split_ids.csv"))
  }
}



if (use_existing_train_split == TRUE){
  iph_train_ids <- read_csv(paste0(data_path, "training_split_ids.csv"))
  
  iph_train <- iph %>% filter(numeric_studyid %in% iph_train_ids$numeric_studyid)
}

iph_test <- iph %>% 
  filter(!(numeric_studyid %in% iph_train$numeric_studyid))