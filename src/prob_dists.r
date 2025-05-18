#
# prob_dists.r
#
# Andrew Roberts
#
# Depends:
#  general_helper_functions.r

library(LaplacesDemon)


# ------------------------------------------------------------------------------
# General helper functions.
# ------------------------------------------------------------------------------


get_param_names <- function(dist_list, flatten=TRUE) {
  # `flatten=FALSE` returns vector of multivariate parameter names. 
  # `flatten=TRUE` returns vector of scalar parameter names, which will be at 
  # least as long as the vector of multivariate parameter names.
  
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


# ------------------------------------------------------------------------------
# Log Density Functions
# ------------------------------------------------------------------------------

get_log_density_func <- function(dist_list) {
  # Assembles a log-density function for a list of independent (multivariate)
  # parameters by adding the log-density functions for each parameter.
  
  param_names <- get_param_names(dist_list, flatten=TRUE)
  ldens_func_list <- lapply(dist_list, get_dist_log_density)
  
  
  ldens <- function(par) {
    if(is.null(dim(par))) par <- matrix(par, nrow=1L)
    vals <- rep(0, nrow(par))
    
    for(i in seq_along(dist_list)) {
      scalar_names <- get_scalar_param_names(dist_list[[i]])
      vals <- vals + ldens_func_list[[i]](par[,scalar_names])
    }
    
    return(vals)
  }
  
  return(ldens)
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
# Sampling functions
# ------------------------------------------------------------------------------

get_sampling_func <- function(dist_list) {
  # Generates a sampling function by concatenating the samples from a 
  # set of (multivariate) parameters.
  
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


convert_par_info_to_list <- function(par_info) {
  # A temporary function to convert the data.frame representation of prior
  # distributions for independent scalar parameters to the list format
  # used by the new PEcAn functions.
  
  prior_list <- list()
  
  for(i in 1:nrow(par_info)) {
    par_name <- par_info$par_name[i]
    dist_name <- par_info$dist[i]
    bounds <- c(par_info$bound_lower[i], par_info$bound_upper[i])
    param1 <- par_info$param1[i]
    param2 <- par_info$param2[i]
    
    if(dist_name == "Gaussian") dist_params <- list(mean=param1, sd=param2)
    else if(dist_name == "Uniform") dist_params <- list(min=param1, max=param2)
    else if(dist_name == "Gamma") dist_params <- list(shape=param1, rate=param2)
    else if(dist_name == "Beta") dist_params <- list(shape1=param1, shape2=param2)
    else {
      stop("Distribution `", dist_name, "` not supported by `convert_par_info_to_list()`.")
    }
    
    prior_list[[par_name]] <- list(param_name=par_name, dist_name=dist_name,
                                   constraint=bounds, dist_params=dist_params,
                                   len=1L)
  }
  
  return(prior_list)
}


# ------------------------------------------------------------------------------
# Parameter transformations:
#  Invertible maps for transforming constrained parameters to unbounded 
#  space.
# ------------------------------------------------------------------------------

get_par_map_funcs <- function(dist_list) {
  
  param_names <- get_param_names(dist_list, flatten=TRUE)
  param_group_names <- get_param_names(dist_list, flatten=FALSE)
  par_map_list <- lapply(dist_list, get_dist_par_maps)
  
  # Map constrained to unconstrained.
  par_map <- function(par) {
    if(is.null(dim(par))) par <- matrix(par, nrow=1L, dimnames=list(NULL, names(par)))
    phi_list <- vector(mode="list", length=length(dist_list))

    for(i in seq_along(dist_list)) {
      scalar_names <- get_scalar_param_names(dist_list[[i]])
      phi_list[[i]] <- par_map_list[[i]]$fwd(par[,scalar_names, drop=FALSE])
    }
    
    do.call(cbind, phi_list)
  }
  
  # Map unconstrained to constrained.
  inv_par_map <- function(phi) {
    if(is.null(dim(phi))) phi <- matrix(phi, nrow=1L, dimnames=list(NULL, names(phi)))
    par_list <- vector(mode="list", length=length(dist_list))
    log_det_J <- rep(0, nrow(phi))

    for(i in seq_along(dist_list)) {
      pattern <- paste0("^_", param_group_names[i], "_[0-9]+_$")
      phi_names <- grep(pattern, colnames(phi), value=TRUE)
      par_list[[i]] <- par_map_list[[i]]$inv(phi[,phi_names, drop=FALSE])
      log_det_J <- log_det_J + drop(attr(par_list[[i]], "log_det_J"))
    }

    par <- do.call(cbind, par_list)
    attr(par, "log_det_J") <- log_det_J
    return(par)
  }
  
  return(list(fwd=par_map, inv=inv_par_map))
}


get_dist_par_maps <- function(dist_info) {
  
  param_name <- dist_info$param_name
  scalar_names <- get_scalar_param_names(dist_info)
  constraint <- dist_info$constraint
  
  if(isTRUE(constraint == "simplex")) {
    map_func <- simplex_map
    inv_map_func <- inv_simplex_map
  } else if(length(constraint)==2L) {
    map_funcs <- get_bound_constraint_map_funcs(constraint)
    map_func <- map_funcs$fwd
    inv_map_func <- map_funcs$inv
  } else {
    stop("Constraint ", constraint, " not currently supported.")
  }
  
  fwd <- function(par) {
    phi <- map_func(par)
    colnames(phi) <- paste0("_", param_name, "_", 1:ncol(phi), "_")
    return(phi)
  }
  
  inv <- function(phi) {
    par <- inv_map_func(phi)
    colnames(par) <- scalar_names
    return(par)
  }
  
  
  return(list(fwd=fwd, inv=inv))
}


get_bound_constraint_map_funcs <- function(bounds) {
  # Fetches the default parameter transformation for parameters with bound
  # constraints (lower, upper, or both).
  
  if(is.infinite(bounds[1]) && is.infinite(bounds[2])) {
    fwd <- id_map
    inv <- id_map
  } else if(is.infinite(bounds[1])) {
    fwd <- function(par) log_upper_map(par, b=bounds[2])
    inv <- function(phi) inv_log_upper_map(phi, b=bounds[2])
  } else if(is.infinite(bounds[2])) {
    fwd <- function(par) log_lower_map(par, a=bounds[1])
    inv <- function(phi) inv_log_lower_map(phi, a=bounds[1])
  } else {
    fwd <- function(par) logit_map(par, a=bounds[1], b=bounds[2])
    inv <- function(phi) inv_logit_map(phi, a=bounds[1], b=bounds[2])
  }
  
  return(list(fwd=fwd, inv=inv))
}


logit_map <- function(par, a=0, b=1) {
  # A log-odds (logit) transform that is translated to have domain (a,b).
  # This can be derived as the composition of the translation 
  # phi(x) = (x-a)/(b-a) with the usual logit transform, which accepts 
  # values in (0,1). If `par` is a vector of length > 1, then the transformation
  # is applied elementwise. The function is NOT vectorized over the arguments
  # `a` and `b`; these are scalar bounds.
  
  if(length(a) != 1L) {
    stop("`a` must be a scalar.")
  }
  
  if(length(b) != 1L) {
    stop("`b` must be a scalar.")
  }
  
  if(a >= b) {
    stop("`a` < `b` must hold.")
  }

  log((par - a) / (b - par))
}



inv_logit_map <- function(phi, a=0, b=1) {
  # The inverse map of `logit_map`. See comments for this function for 
  # details. 

  if(length(a) != 1L) {
    stop("`a` must be a scalar.")
  }
  
  if(length(b) != 1L) {
    stop("`b` must be a scalar.")
  }
  
  if(a >= b) {
    stop("`a` < `b` must hold.")
  }
  
  # This creates a deep copy for all supported input types for `phi` (vector,
  # matrix).
  inverse_logit <- phi
  
  # For positive values, use form that avoids numerical overflow.
  sel_geq_0 <- (phi >= 0)
  inverse_logit[sel_geq_0] <- 1/(1 + exp(-phi[sel_geq_0]))
  
  # For negative values, use form that avoids numerical overflow. When 
  # input is very small, e^x/(1+e^x) is approx e^x.
  sel_l_0 <- !sel_geq_0
  inverse_logit[sel_l_0] <- exp(phi[sel_l_0]) / (1 + exp(phi[sel_l_0]))

  # Map to (a,b).
  par <- a + inverse_logit * (b-a)
  attr(par, "log_det_J") <- log(b-a) + log(inverse_logit) + log(1-inverse_logit)
  
  return(par)
}


log_lower_map <- function(par, a=0) {
  if(length(a) != 1) {
    stop("`a` must be a scalar.")
  }
  
  log(par-a)
}


inv_log_lower_map <- function(phi, a=0) {
  if(length(a) != 1L) {
    stop("`a` must be a scalar.")
  }
  
  par <- a + exp(phi)
  attr(par, "log_det_J") <- phi
  
  return(par)
}


log_upper_map <- function(par, b) {
  if(length(b) != 1L) {
    stop("`b` must be a scalar.")
  }
  
  log(b-par)
}


inv_log_upper_map <- function(phi, b) {
  if(length(b) != 1L) {
    stop("`b` must be a scalar.")
  }
  
  par <- b - exp(phi)
  attr(par, "log_det_J") <- phi
  
  return(par)
}


simplex_map <- function(par) {
  # Implements the map used by Stan for unit simplex-valued parameters.
  # Vectorized so that `par` can be a matrix with one row per input. 
  # Since the d-simplex can be represented by d-1 values (due to the 
  # sum-to-one constraint), this map ignores the final value and thus maps 
  # to R^{d-1}, one dimension lower. The inverse map (see `inv_simplex_map()`)
  # accepts values in R^{d-1} and maps back to R^d. Thus, for an input 
  # `par` of dimension (num_params, d), returns a matrix of dimension 
  # (num_params, d-1). The ith row contains the d-1 unconstrained variables 
  # constructed by applying the map to the ith row of `par`.
  
  if(is.null(dim(par))) par <- matrix(par, nrow=1L)
  d <- ncol(par)

  #
  # Map `par` to intermediate variables `z`.
  #

  # If 2-simplex, then there is only one intermediate variable z1, and it
  # equals par1.
  if(d==2L) {
    z <- par[,1L, drop=FALSE]
  } else {
    par <- par[,1:(d-1), drop=FALSE]
    
    # Cumulative sums of each row of the matrix, excluding the final dimension.
    # Row i becomes par1_i, par1_i + par2_i, ..., (par1_i+...+par[d-1]_i).
    lens <- cbind(rep(1,nrow(par)), 
                  1-matrixStats::rowCumsums(par[,1:(d-2),drop=FALSE]))
    
    z <- par/lens
  }
  
  #
  # Map `z` to unconstrained variables `phi`.
  #

  # In some cases, some of the later components account for almost none of the
  # stick length. In these cases, due to numerical rounding, some of the 
  # z variables can be slightly above 1. We round these down to a number 
  # slightly less than 1.
  z[z >= 1] <- 1 - .Machine$double.eps
  phi <- logit_map(z)
  shift <- log(d-seq(1,d-1))
  add_vec_to_mat_rows(shift, phi)
}


inv_simplex_map <- function(phi) {
  # The inverse map to `simplex_map()`, including computation of the 
  # Jacobian determinant. As discussed in the comments in `simplex_map()`, 
  # for a d-simplex, this forward map maps to d-1 unconstrained variables 
  # `phi`. This function takes these d-1 variables and maps back to d variables
  # lying on the d-simplex. The function is vectorized so that `phi` can be 
  # a matrix of shape (num_params, d-1). Returns matrix of shape 
  # (num_params, d). The returned matrix has attribute `det_J`, which is a
  # vector of length num_params storing the absolute determinant of the 
  # Jacobian evaluated at each input parameter in `phi`.
  
  if(is.null(dim(phi))) phi <- matrix(phi, nrow=1L)
  d <- ncol(phi) + 1L
  
  # Map to intermediate variables.
  shift <- -log(d-seq(1,d-1))
  z <- inv_logit_map(add_vec_to_mat_rows(shift, phi))
  
  # Map to simplex.
  par <- matrix(nrow=nrow(phi), ncol=d)
  g <- matrix(nrow=nrow(phi), ncol=d-1)
  par[,1] <- z[,1]
  g[,1] <- 1
  lens <- 1 - par[,1]
  
  for(j in seq(2,d-1)) {
    par[,j] <- lens * z[,j]
    g[,j] <- par[,j] * (1-z[,j])
    lens <- lens - par[,j]
  }
  
  # Construct final variable.
  par[,d] <- 1 - rowSums(par[,1:(d-1), drop=FALSE])
  
  # Construct vector of absolute Jacobian determinants.
  attr(par, "log_det_J") <- rowSums(log(g))
  
  return(par)
}


id_map <- function(par) {
  if(is.null(dim(par))) par <- matrix(par, nrow=1L)
  par
}

inv_id_map <- function(phi) {
  if(is.null(dim(phi))) phi <- matrix(phi, nrow=1L)
  attr(phi, "log_det_J") <- rep(0, nrow(phi))
  return(phi)
}






