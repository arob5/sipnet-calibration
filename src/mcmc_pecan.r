# mcmc_pecan.r
#
# Markov chain Monte Carlo functions for calibrating the parameters of 
# PEcAn models.
#
# Andrew Roberts
#

mcmc_mh <- function(settings, llik, lprior, n_itr, par_init, 
                    prop_settings=NULL, use_pecan=TRUE) {
  # Metropolis-Hastings with Gaussian proposal distribution. The proposal 
  # covariance is adapted by default.

  # Dimension of parameter space. 
  par_init <- drop(par_init)
  d <- length(par_init)

  # Objects to store samples and other iteration-dependent information
  # (llik/lprior evaluations and proposal standard deviations).
  par_samp <- matrix(nrow=n_itr, ncol=d)
  chain_info <- matrix(nrow=n_itr, ncol=d+2L)
  colnames(par_samp) <- llik_em$input_names
  colnames(chain_info) <- c("llik", "lprior", llik_em$input_names)
  
  # Proposal covariance.
  cov_prop <- prop_settings$cov
  log_scale_prop <- log(prop_settings$scale)
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0L
  times_adapted <- 0L
  
  # Iteration 1 is taken to be the initial condition.
  lprior_curr <- lprior(par_init)
  if(is.infinite(lprior_curr)) {
    llik_curr <- -Inf
  } else {
    llik_curr <- eval_pecan_llik(par_init, llik, run_id="itr_1", 
                                 use_pecan=use_pecan, settings=settings)
  }
  
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
          llik_prop <- eval_pecan_llik(par_prop, llik, run_id=paste0("itr_", itr), 
                                       use_pecan=use_pecan, settings=settings)
          lpost_prop <- lprior_prop + llik_prop
          
          # Accept-Reject step.
          alpha <- min(1.0, exp(lpost_prop - lpost_curr))
          
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_curr <- lprior_prop
            llik_curr <- llik_prop
            lpost_curr <- lpost_prop
            accept_count <- accept_count + 1L 
          } else {
            par_samp[itr,] <- par_curr
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && (((itr-1) %% adapt_interval) == 0)) {
          times_adapted <- times_adapted + 1L
          adapt_list <- adapt_MH_prop_cov(cov_prop, log_scale_prop, 
                                          times_adapted, prop_settings,
                                          samp_history=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE])
          cov_prop <- adapt_list$cov
          log_scale_prop <- adapt_list$log_scale
          if(prop_settings$adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
          accept_count <- 0L
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
              itr_curr=itr, par_init=par_init, condition=err))
}


generate_pecan_llik <- function(settings, obs_op, llik) {

  pecan_llik <- function(par, run_id) {
    # Set up paths and write configs.
    prepare_mcmc_model_run(settings, par, run_id)
    
    # Run model at parameter `par`.
    PEcAn.workflow::start_model_runs(settings, write=FALSE)
    
    # Compute the model predicted observable quantity.
    obs_pred <- obs_op(settings, run_id)
    
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
  obs_pred <- obs_op(settings, run_id)
  
  # Return the log-likelihood evaluation.
  llik(obs_pred)
}


prepare_mcmc_model_run <- function(settings, par, run_id) {
  # TODO: need to ensure par has param names attached, or pass the names in.
  
  # Write config function expects `par` to be a data.frame, with column 
  # names set to parameter names.
  par <- as.data.frame.list(par)
  
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


get_cov_prop_diag_sd <- function(cov_prop, log_scale_prop=NULL) {
  prop_sd <- sqrt(diag(cov_prop))
  if(!is.null(log_scale_prop)) {
    prop_sd <- prop_sd * exp(log_scale_prop)
  }
  
  return(prop_sd)
}





