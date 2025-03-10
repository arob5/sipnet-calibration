# mcmc_pecan.r
#
# Markov chain Monte Carlo functions for calibrating the parameters of 
# PEcAn models.
#
# Andrew Roberts
#

run_mcmc <- function(llik, prior_list, n_itr, n_chains=4L, prop_settings) {
  # A wrapper around `mcmc_mh()` that handles the running of multiple chains,
  # and parameter transformations for constrained parameters.
  
  # Get prior density, parameter transformation maps, and sampling function.
  par_names <- get_param_names(prior_list, flatten=TRUE)
  lprior <- get_log_density_func(prior_list)
  par_maps <- get_par_map_funcs(prior_list)
  rprior <- get_sampling_func(prior_list)
  
  # Construct density of transformed (unconstrained) variables.
  lprior_phi <- function(phi) {
    par <- par_maps$inv(phi)
    log_det_J <- attr(par, "log_det_J")
    
    lprior(par) + log_det_J
  }
  
  # Adjust llik to account for parameter transformation.
  llik_phi <- function(phi, ...) {
    llik(par_maps$inv(phi), ...)
  }
  
  # Initial condition.
  par_init <- rprior(1)
  phi_init <- par_maps$fwd(par_init)
  
  mcmc_output <- mcmc_mh(llik_phi, lprior_phi, 
                         par_names=colnames(phi_init), n_itr=n_itr,
                         prop_settings=prop_settings)
  
  return(mcmc_output)
}


mcmc_mh <- function(llik, lprior, par_names, n_itr, par_init, 
                    prop_settings=NULL, settings=NULL) {
  # Metropolis-Hastings with Gaussian proposal distribution. The proposal 
  # covariance is adapted by default.

  browser()
  
  # Dimension of parameter space. 
  par_init <- drop(par_init)
  d <- length(par_names)
  if(is.null(names(par_init))) names(par_init) <- par_names

  # Objects to store samples and other iteration-dependent information
  # (llik/lprior evaluations and proposal standard deviations).
  par_samp <- matrix(nrow=n_itr, ncol=d)
  chain_info <- matrix(nrow=n_itr, ncol=d+2L)
  colnames(par_samp) <- par_names
  colnames(chain_info) <- c("llik", "lprior", par_names)
  
  # Proposal covariance.
  prop_settings <- get_init_mh_prop_settings(prop_settings, d)
  adapt <- prop_settings$adapt_cov || prop_settings$adapt_scale
  cov_prop <- prop_settings$cov_prop
  log_scale_prop <- prop_settings$log_scale_prop
  L_cov_prop <- t(chol(cov_prop))

  # Iteration 1 is taken to be the initial condition.
  par_curr <- par_init
  lprior_curr <- lprior(par_init)
  if(is.infinite(lprior_curr)) {
    llik_curr <- -Inf
  } else {
    llik_curr <- llik(par_init, run_id="itr_1")
  }
  lpost_curr <- lprior_curr + llik_curr
  
  par_samp[1L,] <- par_init 
  chain_info[1L,] <- c(llik_curr, lprior_curr, 
                       get_cov_prop_diag_sd(cov_prop, log_scale_prop))
  
  # Variable to store error condition, if it occurs.
  err <- NULL
  
  tryCatch(
    {
      for(itr in 2L:n_itr) {
        # Random walk Gaussian proposal. 
        par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
        
        # Compute prior. 
        lprior_prop <- lprior(par_prop)
        
        # Immediately reject if proposal has prior density zero (which will 
        # often happen when the prior has been truncated to stay within the 
        # design bounds). 
        if(is.infinite(lprior_prop)) {
          par_samp[itr,] <- par_samp[itr-1,]
        } else {
          # Compute log-posterior approximation. 
          llik_prop <- llik(par_prop, run_id=paste0("itr_", itr))
          lpost_prop <- lprior_prop + llik_prop
          
          # Accept-Reject step.
          alpha <- min(1.0, exp(lpost_prop - lpost_curr))
          
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_curr <- lprior_prop
            llik_curr <- llik_prop
            lpost_curr <- lpost_prop
            prop_settings <- increment_accept_count(prop_settings)
          } else {
            par_samp[itr,] <- par_curr
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && adapt_this_itr(itr, prop_settings$adapt_interval)) {
          prop_settings <- increment_times_adapted(prop_settings)
          prop_settings <- prepare_adapt_settings(prop_settings, cov_prop,
                                                  log_scale_prop, par_samp, itr)
          prop_list <- do.call(adapt_mh_prop_cov, prop_settings)
          cov_prop <- prop_list$cov_prop
          log_scale_prop <- prop_list$log_scale_prop
          if(prop_settings$adapt_cov) L_cov_prop <- prop_list$L_cov_prop
          prop_settings$accept_count <- 0L
        }
        
        # Store chain info for current step.
        chain_info[itr,] <- c(llik_curr, lprior_curr, 
                              get_cov_prop_diag_sd(cov_prop, log_scale_prop))
      }
      
    }, error = function(cond) {
      err <<- cond
      message("mcmc_mh() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }
  )
  
  return(list(samp=list(par=par_samp),
              info=list(dens=chain_info[,1:2], 
                        prop_sd=chain_info[,3:ncol(chain_info),drop=FALSE]),
              log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, 
              par_curr=par_curr, par_prop=par_prop,
              itr_curr=itr, par_init=par_init, prop_settings=prop_settings,
              condition=err))
}


# ------------------------------------------------------------------------------
# MCMC helper functions.
# ------------------------------------------------------------------------------

adapt_mh_prop_cov <- function(cov_prop, log_scale_prop, adapt_cov, adapt_scale, 
                              times_adapted, accept_count, adapt_interval,
                              samp_history, accept_rate_target=0.24, 
                              adapt_factor_exponent=0.8, 
                              adapt_factor_numerator=10, ...) {
  # Follows the adaptive Metropolis scheme used in the Nimble probabilistic 
  # programming language. 
  # Described here: https://arob5.github.io/blog/2024/06/10/adapt-mh/
  # Assumes a proposal covariance of the form exp(2l)*C where l is 
  # `log_scale_prop` and C is `cov_prop`. Both l and C can be adapted.
  # Let `d` denote the dimension of C.
  #
  # Args:
  # cov_prop: the current value of C.
  # log_scale_prop: the current value of l.
  # adapt_cov: logical, whether or not to adapt C.
  # adapt_scale: logical, whether or not to adapt l.
  # time_adapted: the number of times the covariance has been adapted up to this
  #               point (not the number of MCMC iterations).
  # accept_rate_target: numeric, number between 0 and 1.
  # adapt_factor_exponent: "tau" in the above linked writeup.
  # adapt_factor_numeror: "eta" in the above linked writeup.
  
  return_list <- list()
  adapt_factor <- 1 / (times_adapted + 3)^adapt_factor_exponent
  accept_rate <- accept_count/adapt_interval
  
  if(adapt_cov) {
    if(accept_rate > 0) {
      sample_cov_history <- cov(samp_history)
      cov_prop <- cov_prop + adapt_factor * (sample_cov_history - cov_prop)
    } 
    
    # Handle case that `cov_prop` may not be positive definite. 
    L_cov_prop <- try(t(chol(cov_prop)))
    if(inherits(L_cov_prop, "try-error")) {
      message("Problematic proposal cov: ", cov_prop)
      cov_prop <- Matrix::nearPD(cov_prop, ensureSymmetry=TRUE)
      L_cov_prop <- t(chol(cov_prop))
    }
    
    return_list$L_cov_prop <- L_cov_prop
  }
  
  if(adapt_scale) {
    log_scale_factor <- adapt_factor_numerator * adapt_factor * (accept_rate - accept_rate_target)
    log_scale_prop <- log_scale_prop + log_scale_factor
  }
  
  return_list$cov_prop <- cov_prop
  return_list$log_scale_prop <- log_scale_prop
  
  return(return_list)
}

get_init_mh_prop_settings <- function(prop_settings, d) {
  # Initialize the list used to control proposal covariance adaptation in 
  # a Metropolis-Hastings algorithm with Gaussian proposal. Any proposal 
  # adaptation settings not found in `prop_settings` will be provided 
  # defaults. If `cov_prop` is not found in `prop_settings`, then a default
  # initialization value for the proposal covariance is set. Same for the 
  # scaling term `log_scale_prop`. `d` is the dimension of the parameter 
  # space (i.e., the dimension of the proposal covariance).
  
  if(is.null(prop_settings)) prop_settings <- list()
  
  # Initialize counters to zero.
  prop_settings[c("times_adapted", "accept_count")] <- c(0L, 0L)
  
  # Default initial values for the proposal covariance and scaling.
  if(is.null(prop_settings$cov_prop)) {
    prop_settings$cov_prop <- diag(rep(1,d))
  }
  
  if(is.null(prop_settings$log_scale_prop)) {
    prop_settings$log_scale_prop <- log(2.38) - 0.5*log(d)
  }
  
  # Default adaptation settings.
  if(is.null(prop_settings$adapt_cov)) {
    prop_settings$adapt_cov <- TRUE
  }
  if(is.null(prop_settings$adapt_scale)) {
    prop_settings$adapt_scale <- TRUE
  }
  if(is.null(prop_settings$adapt_interval)) {
    prop_settings$adapt_interval <- 200L
  }
  if(is.null(prop_settings$accept_rate_target)) {
    prop_settings$accept_rate_target <- 0.24
  }
  if(is.null(prop_settings$adapt_factor_exponent)) {
    prop_settings$adapt_factor_exponent <- 0.8
  }
  if(is.null(prop_settings$adapt_factor_numerator)) {
    prop_settings$adapt_factor_numerator <- 10
  }
  
  return(prop_settings)
}


prepare_adapt_settings <- function(prop_settings, cov_prop, log_scale_prop, 
                                   par_samp, itr) {
  # Assembles a list of arguments to be passed to the function 
  # `adapt_mh_prop_cov`. Intended to be called right before a call to this 
  # function. In particular, this function adds the arguments "cov_prop", 
  # "log_scale_prop", and "samp_history" to the `prop_settings` list.
  
  # Determines the length of the sample history that will be used to update
  # the covariance estimator.
  adapt_interval <- prop_settings$adapt_interval
  
  prop_settings$cov_prop <- cov_prop
  prop_settings$log_scale_prop <- log_scale_prop
  prop_settings$samp_history <- par_samp[(itr-adapt_interval+1L):itr,,drop=FALSE]
  
  return(prop_settings)
}


increment_times_adapted <- function(prop_settings) {
  prop_settings$times_adapted <- prop_settings$times_adapted + 1L
  return(prop_settings)
}

increment_accept_count <- function(prop_settings) {
  prop_settings$accept_count <- prop_settings$accept_count + 1L
  return(prop_settings)
}


adapt_this_itr <- function(itr, adapt_interval) {
  # Proposal adaptation is run every `adapt_iterval` iterations. Check if
  # adaptation should be run at iteration `itr`.
  
  ((itr-1) %% adapt_interval) == 0
}


get_cov_prop_diag_sd <- function(cov_prop, log_scale_prop=NULL) {
  # A helper function to extract the standard deviations from the covariance
  # matrix used in the Gaussian proposal of a Metropolis-Hastings algorithm.
  # In particular, assumes a parameterization of the form 
  # Cov = exp(2l)*C; here, l is `log_scale_prop` and C is `cov_prop`.
  # The standard deviations (along the diagonal) are thus of the form 
  # exp(l)*sqrt(C_ii) where C_ii is a diagonal element of C. 
  
  prop_sd <- sqrt(diag(cov_prop))
  if(!is.null(log_scale_prop)) {
    prop_sd <- prop_sd * exp(log_scale_prop)
  }
  
  return(prop_sd)
}


# ------------------------------------------------------------------------------
# Functions for constructing log-likelihood that depends on PEcAn model run.
# ------------------------------------------------------------------------------

generate_pecan_llik <- function(settings, obs_op, llik) {

  pecan_llik <- function(par, run_id) {
    # Set up paths and write configs.
    prepare_mcmc_model_run(settings, par, run_id)
    
    # Run model at parameter `par`.
    PEcAn.workflow::start_model_runs(settings, write=FALSE)
    
    # Compute the model predicted observable quantity.
    obs_pred <- obs_op(run_id, settings)
    
    # Evaluate log-likelihood.
    llik(obs_pred)
  }
  
}


eval_pecan_llik <- function(par, llik, run_id=NULL, settings=NULL, obs_op=NULL, 
                            use_pecan=TRUE) {
  # If `use_pecan` is TRUE 
  
  # If `use_pecan=FALSE`, then `llik` interpreted directly as log-likelihood 
  # function; no PEcAn model run required.
  if(!use_pecan) return(llik(par))
  
  # Otherwise, proceed with PEcAn run model workflow.
  prepare_mcmc_model_run(settings, par, run_id)
  
  # Run model at parameter `par`.
  PEcAn.workflow::start_model_runs(settings, write=FALSE)
  
  # Compute the model predicted observable quantity.
  obs_pred <- obs_op(run_id, settings)
  
  # Return the log-likelihood evaluation.
  llik(obs_pred)
}


prepare_mcmc_model_run <- function(settings, par, run_id) {
  # TODO: need to ensure par has param names attached, or pass the names in.
  
  # Write config function expects `par` to be a data.frame, with column 
  # names set to parameter names.
  if(is.atomic(par)) {
    par <- as.data.frame.list(par)
  } else if(is.data.frame(par)) {
    # TODO: this is a hack as SIPNET's write config function expects `trait.values`
    # to be a list of data.frames (one per PFT). The write config function 
    # should probably be modified to allow either a list or a single data.frame.
    par <- list(par)
  }
  
  # Create run and output directories.
  dir.create(file.path(settings$rundir, run_id), recursive=TRUE)
  dir.create(file.path(settings$modeloutdir, run_id), recursive=TRUE)
  
  # Write config to file.
  model_write_config <- paste("write.config.", settings$model$type, sep = "")
  do.call(model_write_config, args=list(defaults = settings$pfts,
                                        trait.values = par,
                                        settings = settings,
                                        run.id = run_id))
  
  # Write the run ID to a new line in the runs text file. This will override 
  # any existing "runs.txt" file.
  cat(as.character(run_id),
      file = file.path(settings$rundir, "runs.txt"),
      sep = "\n",
      append=FALSE)
}


