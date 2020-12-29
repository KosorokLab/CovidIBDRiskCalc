
################################################################################
##
## Functions for Determining Variables to Include
##
################################################################################

#' Takes a vector of variables to consider and a prevalence cut off, returns a character vector 
#' where binary variables below the prevalence cut off have been removed. 
#' Categorical (with 3 or more levels) and continuous variables are always included. 
#' @param ibd.dat data set
#' @param candidate.vars character vector of all variables being considered
#' @param prevalence.cutoff a number in (0,1)
#' @param count.cutoff minimum count
DetermineVariablesForConsideration <- function(ibd.dat, candidate.vars, prevalence.cutoff = NULL, count.cutoff = NULL){
  if (is_null(prevalence.cutoff) == FALSE) stopifnot(prevalence.cutoff >=0 & prevalence.cutoff <= 1)
  if (is_null(count.cutoff) == FALSE) stopifnot(is.numeric(count.cutoff))
  
  
  iph_lgl_data <- ibd.dat %>% select(all_of(candidate_vars)) %>% select_if(is.logical) 
  
  iph_lgl_count <- colSums(iph_lgl_data == TRUE, na.rm = TRUE)
  
  if (is_null(prevalence.cutoff) == TRUE) prevalence.cutoff <- 0
  if (is_null(count.cutoff) == TRUE) count.cutoff <- 0
  
  low_prevalence_vars <- names(iph_lgl_count)[iph_lgl_count/nrow(iph_lgl_data) <= prevalence.cutoff]
  
  low_count_vars <- names(iph_lgl_count)[iph_lgl_count <= count.cutoff]
  
  vars_to_exclude <- c(low_prevalence_vars, low_count_vars)
  
  candidate_vars_with_min_prev <- candidate_vars[!candidate_vars %in% vars_to_exclude]
  
  return(candidate_vars_with_min_prev)
}