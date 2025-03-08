#
# prob_dists.r
#
# Andrew Roberts
#

library(LaplacesDemon)

# TODO: need a function to set default names. Currently only sets default names
# for scalar-valued parameters.
get_param_names <- function(dist_list, flatten=TRUE) {
  
  if(flatten) {
    names_list <- lapply(dist_list, function(x) x$scalar_names)
    for(i in seq_along(names_list)) {
      if(is.null(names_list[[i]])) {
        names_list[[i]] <- dist_list[[i]]$param_name
      }
    }
    
    param_names <- do.call(c, names_list)
    
  } else {
    param_names <- sapply(dist_list, function(x) x$param_name)
  }
  
  return(param_names)
}


get_sampling_function <- function(dist_list) {
  
  param_names <- get_param_names(dist_list, flatten=TRUE)
  sampling_func_list <- lapply(dist_list, get_dist_sampling_func)
  
  sampling_func <- function(n) {
    samp_list <- lapply(sampling_func_list, function(f) f(n))
    samp <- do.call(cbind, samp_list)
    colnames(samp) <- param_names
    return(samp)
  }
 
  return(sampling_func) 
}


get_dist_sampling_func <- function(dist_info) {
  # `dist_info` must have arguments "dist_name" and "dist_params".
  
  dist_name <- dist_info$dist_name
  dist_params <- dist_info$dist_params
  
  if(dist_name=="Uniform") {
    func_name <- "runif"
  } else if(dist_name=="Gamma") {
    func_name <- "rgamma"
  } else if(dist_name=="Beta") {
    func_name <- "rbeta"
  } else if(dist_name=="Normal") {
    func_name <- "rnorm"
  } else if(dist_name=="Dirichlet") {
    func_name <- "rdirichlet"
  } else {
    stop("Distribution ", dist_name, " not currently supported.")
  }
  
  sampling_func <- function(n) {
    dist_params$n <- n
    do.call(func_name, dist_params)
  }
  
  return(sampling_func)
}







