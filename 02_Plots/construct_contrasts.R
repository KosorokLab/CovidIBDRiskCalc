################################################################################
##  Functions for using contrasts
##  - Read all the csvs in a directory and create a single contrast matrix
## 
################################################################################

#' Given a directory containing a number of csvs with contrast matrices, reads 
#' all the csvs in the directory and creates a single matrix. 
ConstructSingleMatrixFromContrastDir <- function(contrast.dir){
  contrast_filenames <- list.files(path = contrast.dir, 
                                   pattern = ".csv", 
                                   full.names = TRUE, 
                                   recursive = FALSE)
  
  # Source all of the outcome definitions without printing them
  contrast_mat_list <- map(contrast_filenames, readr::read_csv)
  
  combined_contrast_mat <- MergeContrastMatricesFromList(contrast_mat_list)
  
  return(combined_contrast_mat)
}

MergeContrastMatricesFromList <- function(contrast.mat.list){
  # Check that the variables are all the same
  stopifnot(all(map_dbl(contrast.mat.list, nrow) == nrow(contrast.mat.list[[1]])))
  stopifnot(all(map_lgl(contrast.mat.list, ~all(.$Variable == contrast.mat.list[[1]]$Variable))))
  
  combined_contrast_mat <- map(contrast.mat.list, ~select(., -Variable)) %>% 
    bind_cols(.) %>% 
    as.matrix
  
  row.names(combined_contrast_mat) <- contrast.mat.list[[1]]$Variable
  return(combined_contrast_mat)
}
