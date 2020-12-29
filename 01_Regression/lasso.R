################################################################################
## Functions for constructing a formula for use with glmnetUtils
##
################################################################################


#' Constructs an R formula for use with glmnetUtils given character vectors 
#' containing the outcome name, covariate names, and interaction terms. 
#' Optionally saves the resulting formula as an .rds file
#' @param outcome.name string containing the outcome variable name
#' @param covariate.names character vector where each entry is a covaraite to include in the model
#' @param interaction.terms character vector where each entry is a string specifying an interaction terms like 'var1:var2'. Use colon notation and not *
#' @param save.file.name It specified, the formula will be written out to the specified file. 
#' If unspecified the function simply returns the formula
#' @return an R formula constructed from the outcome, covariate, and interaction strings
ConstructLassoFormula <- function(outcome.name, covariate.names, interaction.terms, 
                                  save.file.name = NULL,
                                  terms.to.remove = NULL) {
  fit_string <- paste(covariate.names, collapse = " + ") %>% 
    paste(outcome.name, ., sep = " ~ ")
  
  if (is.null(interaction.terms) == FALSE) {
    fit_string <- paste(interaction.terms, collapse = " + ") %>% 
      paste(fit_string, ., sep = " + ")
  }
  
  if (is.null(terms.to.remove) == FALSE) {
    fit_string <- paste(terms.to.remove, collapse = " - ") %>% 
      paste(fit_string, ., sep = " - ")
  }
  
  fit_formula <- as.formula(fit_string)
  
  if (is_empty(save.file.name) == FALSE) saveRDS(fit_formula, file= save.file.name)
  
  return(fit_formula)
}

#' This function takes \code{glmnetUtil}'s \code{makeModelComponents} function
#' and modifies it to remove redundant "other country" factor levels present in our data
ConstructLassoDataMatrix <- function(in.formula, in.data, 
                                     names.to.remove = "stateCatMixWOtherOther country",
                                     create.intercept = TRUE, ...){
  return(ConstructLassoDataList(in.formula = in.formula,
                                in.data = in.data,
                                names.to.remove = names.to.remove,
                                create.intercept = create.intercept)$x)
}

#' This function takes \code{glmnetUtil}'s \code{makeModelComponents} function
#' and modifies it to remove redundant "other country" factor levels present in our data.
#' Returns a list with the data matrix and the response vector
ConstructLassoDataList <- function(in.formula, in.data, 
                                   names.to.remove = "stateCatMixWOtherOther country",
                                   create.intercept = FALSE,
                                   ...){
  
  model.comp <- glmnetUtils:::makeModelComponents(formula = in.formula,
                                                  data = in.data)
  
  x_mat <- model.comp$x[, !colnames(model.comp$x) %in% names.to.remove]
  
  # If there is only a single row model.comp$x is a vector not a matrix
  if (nrow(in.data) == 1) {
    x_mat <- matrix(data = x_mat, nrow = 1)
    colnames(x_mat) <- colnames(model.comp$x)[!colnames(model.comp$x) %in% names.to.remove]
  } 
  
  y_vec <- model.comp$y
  
  if (create.intercept == TRUE){
    x_mat <- cbind(rep(1, times = nrow(x_mat)), x_mat)
    colnames(x_mat)[1] <- "(Intercept)"
  } 
  
  if (any(duplicated(colnames(x_mat))) == TRUE) stop("Formula resulted in duplicate column names, check that there aren't repeated terms in the formula")
  
  lasso_data_list <- list(x = x_mat,
                          y = y_vec)
  
  return(lasso_data_list)
}

################################################################################
## Wrapper functions which construct a data matrix then call the relevant
## glmnetUtils or glmnet function
##
################################################################################


#' Wrapper around cva.glmnet that uses our own data matrix construction to ensure
#' the redundant other country levels are removed
CvaGlmnet<- function(in.formula, in.data, ...){
  data_list <- ConstructLassoDataList(in.formula = in.formula,
                                in.data = in.data,
                                ...)
  
  cva_fit <- glmnetUtils::cva.glmnet(x = data_list$x,
                                     y = data_list$y,
                                     ...)
  return(cva_fit)
  
}

#' Wrapper around cv.glmnet that uses our own data matrix construction to ensure
#' the redundant other country levels are removed
CvGlmnet <- function(in.formula, in.data, ...){
  data_list <- ConstructLassoDataList(in.formula = in.formula,
                                      in.data = in.data,
                                      ...)
  
  cva_fit <- glmnet::cv.glmnet(x = data_list$x,
                                     y = data_list$y,
                                     ...)
  return(cva_fit)
}

#' Wrapper around glmnet that uses our own data matrix construction to ensure
#' the redundant other country levels are removed
Glmnet <- function(in.formula, in.data, ...){
  data_list <- ConstructLassoDataList(in.formula = in.formula,
                                      in.data = in.data,
                                      ...)
  
  cva_fit <- glmnet::glmnet(x = data_list$x,
                               y = data_list$y,
                               ...)
  return(cva_fit)
}

################################################################################
##
## Lasso-related Utility functions
##
################################################################################

#' Create \math{d \times 2} data frame where each row contains the variable name and the coefficient value
#' @param lambda.to.use a string specifying which criteria to use for selecting the final model. 
#' Must be either 'lambda.min', in which case the final model is the one minimizing the cv-error, or 'lambda.1se'
#' in which case the final model is the sparsest model with cv-error within one standard deviation of the minimal cv-error
GetCVLassoCoefs <- function(cv_mod, lambda.to.use = "lambda.min") {
  if (lambda.to.use != "lambda.min" & lambda.to.use != "lambda.1se") stop("lambda.to.use must be one of 'lambda.min' or 'lambda.1se'")
  
  main_coef <- coef(cv_mod, s = lambda.to.use) %>% as.matrix()
  
  coef_df <- data.frame(name = names(main_coef[main_coef[,1] != 0,]), 
                        main_coef[main_coef[,1] != 0,], row.names = NULL,
                        stringsAsFactors = FALSE)
  colnames(coef_df) <- c("name", lambda.to.use)
  
  
  return(coef_df)
}