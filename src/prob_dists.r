#
# prob_dists.r
#
# Andrew Roberts
#

library(LaplacesDemon)


get_param_names <- function(dist_list, flatten=TRUE) {
  
  if(flatten) {
    names_list <- lapply(dist_list, get_scalar_param_names)
    param_names <- do.call(c, names_list)
  } else {
    param_names <- sapply(dist_list, function(x) x$param_name)
  }
  
  return(param_names)
}


get_scalar_param_names <- function(dist_info) {
  # Extracts the scalar parameter names from a single multivariate parameter.
  # Creates default scalar names if they are not explicitly provided in 
  # `dist_info`.
  
  if(!is.null(dist_info$scalar_names)) return(dist_info$scalar_names)
  if(dist_info$len==1L) return(dist_info$param_name)
  
  paste(dist_info$param_name, seq_len(dist_info$len), sep="_")
}

get_sampling_func <- function(dist_list) {
  
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


get_log_density_func <- function(dist_list) {
  
  param_names <- get_param_names(dist_list, flatten=TRUE)
  ldens_func_list <- lapply(dist_list, get_dist_log_density)
  
  
  lprior <- function(par) {
    if(is.null(dim(par))) par <- matrix(par, nrow=1L)
    vals <- rep(0, nrow(par))
    
    for(i in seq_along(dist_list)) {
      scalar_names <- get_scalar_param_names(dist_list[[i]])
      vals <- vals + ldens_func_list[[i]](par[,scalar_names])
    }
    
    return(vals)
  }
  
  return(lprior)
}


get_dist_log_density <- function(dist_info) {
  # We define the log density function such that:
  # 1.) It can accept a vector (indicating a single input value), or 
  #     a matrix (one row per input value). Returns a vector of log-density
  #     evaluations of length equal to the number of input values.
  # 2.) Returns -Inf for any values outside of the distribution support.

  dist_name <- dist_info$dist_name
  dist_params <- dist_info$dist_params
  
  if(dist_name=="Uniform") {
    func <- "dunif"
  } else if(dist_name=="Gamma") {
    func <- "dgamma"
  } else if(dist_name=="Beta") {
    func <- "dbeta"
  } else if(dist_name=="Normal") {
    func <- "dnorm"
  } else if(dist_name=="Dirichlet") {
    func <- get_dirichlet_ldens_func(dist_params)
  } else {
    stop("Distribution ", dist_name, " not currently supported.")
  }
  
  ldens_func <- function(x) {
    dist_params$x <- x
    dist_params$log <- TRUE
    do.call(func, dist_params)
  }
}


get_dirichlet_ldens_func <- function(dist_params) {
  # The LaplacesDemon::ddirichlet function seems to have some unfortunate
  # behavior. Instead of returning 0 when a negative value is input, it 
  # throws an error. Yet it doesn't seem to throw an error when given an input 
  # that doesn't sum to 1 (and also doesn't return zero). When given an input
  # on the boundary (e.g., [1,0,0]), the function sometimes returns Inf and 
  # sometimes returns NaN, depending on the value of `alpha`.
  # Not sure what is going on here, but we instead try to handle these issues 
  # manually in this function. Instead of passing a matrix to 
  # `ddirichlet`, this function is called for each output individually to 
  # account for the fact that it throws an error if even one of the inputs
  # has a negative value.
  
  alpha <- dist_params$alpha
  if(is.null(alpha)) {
    stop("`alpha` parameter required for Dirichlet distribution.")
  }
  
  ldens <- function(x, alpha, log=FALSE) {
    if(is.null(dim(x))) x <- matrix(x, nrow=1L)
    vals <- vector(mode="numeric", length=nrow(x))
    tol <- sqrt(.Machine$double.eps)

    for(j in 1:nrow(x)) {
      if(abs(sum(x[j,]) - 1) > tol) val <- -Inf
      else val <- try(LaplacesDemon::ddirichlet(x[j,], alpha, log=TRUE), silent=TRUE)
      
      if(class(val)=="try-error") val <- -Inf
      vals[j] <- val
    }
    
    vals[is.na(vals) | is.infinite(vals)] <- -Inf
    return(vals)
  }
  
  return(ldens)
}


# ------------------------------------------------------------------------------
# Parameter transformations:
#  Invertible maps for transforming constrained parameters to unbounded 
#  space.
# ------------------------------------------------------------------------------

get_dist_default_par_map <- function() {
  
}


logit_map <- function(par, a=0, b=1) {
  # A log-odds (logit) transform that is translated to have domain (a,b).
  # This can be derived as the composition of the translation 
  # phi(x) = (x-a)/(b-a) with the usual logit transform, which accepts 
  # values in (0,1). If `par` is a vector of length > 1, then the transformation
  # is applied elementwise. In this case, the bounds `a` and `b` can either be:
  # 1.) scalars, in which the same bounds are used for every entry of `par`.
  # 2.) vectors of length equal to the length of `par`, in which case potentially
  #     different bounds are applied to each entry.
  
  if((length(a) != 1) & (length(a) != length(par))) {
    stop("`a` must be length 1 or length equal to `par`.")
  }
  
  if((length(b) != 1) & (length(b) != length(par))) {
    stop("`b` must be length 1 or length equal to `par`.")
  }
  
  if(any(a >= b)) {
    stop("`a` < `b` must hold (elementwise).")
  }
  
  log((par - a) / (b - par))
}


inv_logit_map <- function(phi, a=0, b=1) {
  # The inverse map of `logit_map`. See comments for this function for 
  # details. 
  
  if((length(a) != 1) & (length(a) != length(par))) {
    stop("`a` must be length 1 or length equal to `par`.")
  }
  
  if((length(b) != 1) & (length(b) != length(par))) {
    stop("`b` must be length 1 or length equal to `par`.")
  }
  
  if(any(a >= b)) {
    stop("`a` < `b` must hold (elementwise).")
  }
  
  a + (b - a) / (1 + exp(-phi))
}


log_lower_map <- function(par, a=0) {
  if((length(a) != 1) & (length(a) != length(par))) {
    stop("`a` must be length 1 or length equal to `par`.")
  }
  
  log(par-a)
}


inv_log_lower_map <- function(phi, a=0) {
  if((length(a) != 1) & (length(a) != length(par))) {
    stop("`a` must be length 1 or length equal to `par`.")
  }
  
  a + exp(phi)
}

log_upper_map <- function(par, b) {
  if((length(b) != 1) & (length(b) != length(par))) {
    stop("`b` must be length 1 or length equal to `par`.")
  }
  
  log(b-par)
}


inv_log_upper_map <- function(phi, b) {
  if((length(b) != 1) & (length(b) != length(par))) {
    stop("`b` must be length 1 or length equal to `par`.")
  }
  
  b - exp(phi)
}


