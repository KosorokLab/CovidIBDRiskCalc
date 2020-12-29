################################################################################
################################################################################
## Functions related to model-fitting on bootstrap samples
##
################################################################################

#' Single iteration of fitting glmnet to a bootstrapped data sample
#' @param i not used, to facilitate working with lapply
#' @param covariate.data covariate data to use
#' @param resp.mat matrix where each column is an outcome to fit a (separate) model to
#' @param outcome.names Character vector whose entries are the variable names for the outcomes in the data set
#' @param lambda.seq.list List whose entries are numeric vectors with lambda values to use for model fitting
#' @return list whose entries are glmnet fitted models
BootGlmnetSimultaneousIter <- function(bootstrap.indices, 
                                       covariate.data, 
                                       resp.mat, 
                                       outcome.names,
                                       lambda.seq.list = NULL,
                                       ...){
  if (is_null(lambda.seq.list) == FALSE & 
      all(outcome.names == names(lambda.seq.list)) == FALSE) {
    stop("Elements of lambda.seq.list must be named and the names must match outcome.names")
  }
  
  boot_cov_data <- covariate.data[bootstrap.indices,] %>% as.matrix
  
  boot_resp_data <- resp.mat[bootstrap.indices,]
  
  if(is_null(lambda.seq.list) == TRUE) {
    glmnet_fit_list <- map(.x = boot_resp_data, 
                           function(y)
                             glmnet::glmnet(x = boot_cov_data, 
                                            y = y,
                                            family = "binomial", 
                                            type.measure = "deviance",
                                            alpha = 1,
                                            standardize = TRUE,
                                            intercept = TRUE,
                                            ...))
  }
  
  if(is_null(lambda.seq.list) == FALSE) {
    glmnet_fit_list <- map2(.x = boot_resp_data, 
                            .y = lambda.seq.list,
                            function(y, lambda.seq, ...)
                              glmnet::glmnet(x = boot_cov_data, 
                                             y = y,
                                             family = "binomial", 
                                             type.measure = "deviance",
                                             standardize = TRUE,
                                             intercept = TRUE,
                                             lambda = lambda.seq,
                                             ...))
  }
  
  return(glmnet_fit_list)
}

#' This is built for running the HospOrDeath, IcuVentDeath, and death models 
#' on the same bootstrap sample. 
#' @param mice.data mice data object; gets transformed into a list of imputed data sets by the function
#' @param n.iter Number of bootstrap iterations. Total models fit = n.iter * number of imputed data sets
#' @param outcome.names Character vector whose entries are the variable names for the outcomes in the data set
#' @param parallel Logical indicating 
#' @param lambda.seq.list List whose entries are numeric vectors with lambda values to use for model fitting
BootMiceSimultaneous <- function(mice.data, 
                                 n.iter, 
                                 formula.list,
                                 outcome.names = c("HospOrDeath", "IcuVentDeath", "death"),
                                 lambda.seq.list = NULL,
                                 parallel = FALSE, ...){
  stopifnot(is_null(lambda.seq.list) | (length(lambda.seq.list) == length(outcome.names)))
  stopifnot(length(formula.list) == length(outcome.names))
  
  imputed_data_list <- mice.data %>% mice::complete("all")
  
  # Check that the covariates are the same for every model
  data_for_checking_cols <- map(formula.list, ~ConstructLassoDataMatrix(in.formula = .x, in.data = imputed_data_list[[1]]))
  
  if (suppressMessages(all_equal(data_for_checking_cols[[1]], data_for_checking_cols[[2]]) == FALSE |
                       all_equal(data_for_checking_cols[[1]], data_for_checking_cols[[3]]) == FALSE)) stop("This function requires the design matrices to be the same for all models currently")
  
  # Get a list of covariate matrices, one entry for each imputed data set
  # and a list of response matrices, one for each imputed data set, to feed into map2
  data_mat_list <- map(.x = imputed_data_list,
                       ~ConstructLassoDataMatrix(in.formula = formula.list[[1]], 
                                                 in.data = .x,
                                                 create.intercept = FALSE))
  
  resp_mat_list <- map(imputed_data_list, ~.x %>% select(all_of(outcome.names)))
  
  # Indices should be shared across all imputed data sets for a single bootstrap replication
  row_indices <- 1:nrow(imputed_data_list[[1]])
  
  bootstrap_index_sets <- replicate(n = n.iter, 
                                    sample(row_indices, 
                                           size = length(row_indices), 
                                           replace = TRUE),
                                    simplify = FALSE)
  if (parallel == FALSE) {
    data_mods_std <- 
      map(.x = bootstrap_index_sets,
          function(cur.indices) 
            map2(.x = data_mat_list, .y = resp_mat_list, 
                 function(x, y)
                   BootGlmnetSimultaneousIter(bootstrap.indices = cur.indices, 
                                              covariate.data = x, 
                                              resp.mat = y,
                                              outcome.names = outcome.names,
                                              lambda.seq.list = lambda.seq.list, ...)))
  }
  
  if (parallel == TRUE) {
    data_mods_std <-
      furrr::future_map(.x = bootstrap_index_sets,
                        function(cur.indices) 
                          map2(.x = data_mat_list, .y = resp_mat_list, 
                               function(x, y)
                                 BootGlmnetSimultaneousIter(bootstrap.indices = cur.indices, 
                                                            covariate.data = x, 
                                                            resp.mat = y,
                                                            outcome.names = outcome.names,
                                                            lambda.seq.list = lambda.seq.list, ...)))
    
  }
  
  return(data_mods_std)
}

#' Takes a list created by \code{BootMiceSimultaneous} which is [B,M,3 (# models)]
#' and returns a list which is [3, B, M]
TransformBootstrapListToModelList <- function(boot.mod.run){
  b_bootstrap_reps <- length(boot.mod.run)
  m_imputed_dfs <- length(boot.mod.run[[1]])
  
  #Pre-allocate the lists
  hosp_mod_list <- vector("list", length = b_bootstrap_reps)
  icu_mod_list <- vector("list", length = b_bootstrap_reps)
  death_mod_list <- vector("list", length = b_bootstrap_reps)
  
  for(i in 1:b_bootstrap_reps){
    hosp_mod_list[[i]] <- vector("list", length = m_imputed_dfs)
    icu_mod_list[[i]] <- vector("list", length = m_imputed_dfs)
    death_mod_list[[i]] <- vector("list", length = m_imputed_dfs)
  }
  
  # Nasty nested loops but this is only assignment so I don't think it's a big issue
  for(i in 1:b_bootstrap_reps){
    for(j in 1:m_imputed_dfs){
      hosp_mod_list[[i]][[j]] <- boot.mod.run[[i]][[j]][[1]]
      icu_mod_list[[i]][[j]] <- boot.mod.run[[i]][[j]][[2]]
      death_mod_list[[i]][[j]] <- boot.mod.run[[i]][[j]][[3]]
    }
  }
  
  reordered_list <- list(hosp_mod_list, icu_mod_list, death_mod_list)
  return(reordered_list)
}
################################################################################
################################################################################
## Functions for prediction based off of fitted bootstrap models
##
################################################################################


#' Creates a \eqn{p \times b} matrix of coefficients from a list of bootstrap replicates
#' where p is the number of parameters considered in the model 
#' (not the number of nonzero coefficients in the final model)
#' and b is the number of bootstrap replicates
#' @param glmnet.boot.list
#' @param lambda.to.use Lambda to use for selecting the model coefficients
CreateBootCoefMatrix <- function(glmnet.boot.list, lambda.to.use){
  # Data are row vectors and the coefficient vectors are column vectors
  glmnet_coef_mat <- map(glmnet.boot.list, ~as.matrix(coef(., s = lambda.to.use))) %>% 
    unlist %>% 
    matrix(., ncol = length(glmnet.boot.list), byrow = FALSE)
  
  rownames(glmnet_coef_mat) <- coef(glmnet.boot.list[[1]], lambda.to.use) %>% 
    as.matrix %>% 
    rownames
  
  return(glmnet_coef_mat)
}

#' Creates an \eqn{n \times b} matrix of coefficients from a list of bootstrap replicates
#' where p is the number of parameters considered in the model 
#' (not the number of nonzero coefficients in the final model)
#' and b is the number of bootstrap replicates
#' @param boot_coef_mat  \eqn{p \times b} matrix of coefficients as created by \code{CreateBootCoefMatrix}
#' @param mod_formula Lambda to use for selecting the model coefficients
#' @param new_obs a vector or matrix of new observations. Should only include the original data, 
#' the transformations and interactions are handled within the function
BootPredictMatrix <- function(boot_coef_mat, mod_formula, new_obs, ...){
  
  if(is.factor(new_obs$race___0) == FALSE) {
    new_obs <- new_obs %>% 
      mutate_at(., .vars = paste0("race___", c(0, 1, 3)), ~factor(., levels = c("FALSE", "TRUE")))
  }
  
  new_obs_matrix <- ConstructLassoDataMatrix(in.formula = mod_formula, 
                                                in.data = new_obs,
                                                create.intercept = TRUE)
  
  if( all(colnames(new_obs_matrix) == rownames(boot_coef_mat))  == FALSE){
    nonmatching_new_obs <- colnames(new_obs_matrix)[!colnames(new_obs_matrix) %in% rownames(boot_coef_mat)]
    nonmatching_coef_mat <- rownames(boot_coef_mat)[!rownames(boot_coef_mat) %in% colnames(new_obs_matrix)]
    
    if (is_empty(nonmatching_new_obs) & is_empty(nonmatching_coef_mat)) stop(
      "The names in the data matrix and the predictor matrix do not match, but the variables are the same - check order")
    
    nonmatching_new_obs_string <- paste(nonmatching_new_obs, collapse = "\n") %>% 
      paste("Variables in the new observation matrix not in the coefficient matrix:", ., sep = "\n")
    nonmatching_coef_mat_string <- paste(nonmatching_coef_mat, collapse = "\n") %>% 
      paste("Variables in the new observation matrix not in the coefficient matrix:", ., sep = "\n")
    
    error_msg <- paste("The names in the data matrix and the predictor matrix do not match.", 
                       nonmatching_new_obs_string, nonmatching_coef_mat_string, sep =  "\n")
    stop(error_msg)
  }
  
  new_pred_log_scale <- new_obs_matrix %*% boot_coef_mat
  new_pred <- boot::inv.logit(new_pred_log_scale)
  return(new_pred)
}

#' Average over the predicted probabilities from all the imputed data sets for 
#' each of the bootstrap replications. This assumes that the predictions came from 
#' an ordered set of coefficients so that sequential groups of columns of size \code{m.imputations}
#' correspond to a single bootstrap sample of row indices fit to each of the imputed data sets. 
#' @param pred.mat predictions created by \code{BootPredictMatrix}
#' @param m.imputations integer equal to the number of imputed data sets used in the model fitting
#' @return a n times b matrix whose entries are the predicted probabilities per bootstrap replication
#' averaged over the imputed data sets
BootPredictMatrixAvgOverImp <- function(pred.mat, m.imputations = 30){
  stopifnot(is.numeric(m.imputations))
  # Want a matrix not a vector even if there is only a single observation
  if (is.null(nrow(pred.mat))) pred.mat <- matrix(pred.mat, nrow = 1)
  
  index_begin_seq <- seq(from = 1, 
                         to = ncol(pred.mat),
                         by = m.imputations)
  
  flat_pred_avg_over_imps <- matrix(nrow = nrow(pred.mat), 
                                    ncol = length(index_begin_seq))
  
  # This is way, way faster than map
  # With 400 bootstrap replications and 30 imputed data sets 
  # map took 42 seconds, this took 2.2 seconds
  
  # rowMeans doesn't work on a matrix with a single row
  # There is definitely a more elegant way to deal with this problem, but this'll
  # do for now
  
  if (nrow(pred.mat) == 1){
    for(i in 1:length(index_begin_seq)) {
      flat_pred_avg_over_imps[, i] <- mean(
        pred.mat[, index_begin_seq[[i]]:(index_begin_seq[[i]]+m.imputations-1)])
    }
  } else {
    for(i in 1:length(index_begin_seq)) {
      flat_pred_avg_over_imps[, i] <- rowMeans(
        pred.mat[, index_begin_seq[[i]]:(index_begin_seq[[i]]+m.imputations-1)])
    }
  }
  
  return(flat_pred_avg_over_imps)
}
  