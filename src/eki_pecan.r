

run_eki_pecan <- function(y, fwd, Sig, prior_list, n_itr=1L, U0=NULL, G0=NULL,
                          design_method="LHS", n_ens=NULL, settings=NULL) {
  # Runs Ensemble Kalman inversion (EKI) for `n_itr` iterations. This extends 
  # the one-step update implemented by `run_eki_step()` via a likelihood 
  # tempering approach. At present, this function assumes that the inverse 
  # problem has a Gaussian likelihood with covariance `Sig`, but this can 
  # potentially be generalized in the future. To be explicit, the current 
  # assumption is an additive Gaussian noise model, with independent noise 
  # eps ~ N(0,Sig). Each iteration of EKI requires evaluation of the forward 
  # model `fwd`. The final ensemble returned by this function provides a 
  # Monte Carlo approximation of the posterior p(u|y). For linear Gaussian 
  # inverse problems, these samples will be exactly distributed according to 
  # p(u|y), otherwise this represents a pure approximation. Since this EnKF
  # based method relies on Gaussian approximations, it is typically advisable 
  # to transform the parameter ensemble members to an unconstrained space 
  # prior to performing the EnKF updates. The arguments `transform_pars` and 
  # `par_map` allow for such transformations. See details below.
  # 
  # Args:
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on.
  #    fwd: function, representing the forward model. Must be vectorized so 
  #         that it can accept a matrix with each row representing a different 
  #         parameter vector inputs. It should return a matrix with rows 
  #         corresponding to the different inputs, and number of columns equal 
  #         to `P`, the dimension of the observation space.
  #    Sig: matrix of shape (P,P), the covariance matrix for the Gaussian 
  #         likelihood.
  #    n_itr: integer, the number of iterations the algorithm will be run.
  #    par_prior: data.frame storing the prior distributions for the parameters.
  #               This is required if `U0` is not provided, or if 
  #               `transform_pars` is TRUE but `par_map` is FALSE.
  #    U0: matrix, initial parameter ensemble of dimension (J,D), where J is 
  #        the number  of ensemble members and D the parameter space dimension. 
  #        If not  provided will be sampled from prior.
  #    G0: matrix, representing the output of `fwd(U0)`, the model outputs based
  #        on the initial ensemble. Optional.
  #    transform_pars: if TRUE, will use `par_map()` (or the default transport
  #                    map) to transform the ensemble members to an unconstrained 
  #                    space prior to the application of the EnKF analysis step. 
  #                    The inverse transformation is then applied before executing
  #                    the next round of forward model evaluations.
  #    par_map: function, representing the parameter transformation map (i.e., 
  #             transport map) and its inverse. See comments on the return value
  #             of `get_default_par_map()` for the requirements on this function.
  #    design_method: character, the sampling method used to generate the 
  #                   initial ensemble; only used if `U0` is NULL. See 
  #                   `get_batch_design()` for the different sampling options.
  #    n_ens: integer, the number of ensemble members J. Only required if 
  #           `U0` is NULL. 
  #
  # Returns:
  # list, with elements:
  #    `U`: matrix, the final JxD parameter ensemble in the untransformed 
  #         (original) space.
  #    `par_map`: function, the transport map used in the algorithm. NULL if 
  #               `transform_pars` is FALSE.
  #    `eki_list`: list of length `n_itr`. The kth element is the list returned 
  #                by the call to `run_eki_step()` at the kth iteration of the 
  #                algorithm.

  # If initial ensemble is not provided, sample from prior.
  if(is.null(U0)) {
    assert_that(!is.null(n_ens))
    rprior <- get_sampling_func(prior_list)
    U <- rprior(n_ens)
  } else {
    assert_that(is.matrix(U0))
    n_ens <- nrow(U0)
    U <- U0
  }
  
  # Map to unconstrained space.
  par_maps <- get_par_map_funcs(prior_list)
  U <- par_maps$fwd(U)
  
  # List storing intermediate outputs.
  eki_list <- vector(mode="list", length=n_itr)
  
  for(k in 1:n_itr) {
    
    # Run forward model. On first iteration, don't run if initial forward 
    # model evaluations have been explicitly passed in args.
    if((k==1L) && !is.null(G0)) G <- G0
    else G <- fwd(par_maps$inv(U), k, settings)
    
    # EnKF update with tempered likelihood.
    Sig_scaled <- n_itr * Sig
    eki_step_list <- run_eki_step(U=U, y=y, G=G, Sig=Sig_scaled)
    U <- eki_step_list$U
    eki_list[[k]] <- eki_step_list
  }
  
  # Map back to original space.
  U <- par_maps$inv(U)
  
  return(list(U=U, par_maps=par_maps, eki_list=eki_list))
}


run_eki_step <- function(U, y, G, Sig) {
  # Computes one iteration of Ensemble Kalman inversion, which involves 
  # computing the sample mean and covariance estimates using the current 
  # ensembles, and then calling `compute_enkf_update()`. This function assumes
  # that the forward model has already run at the ensemble members `U`, with 
  # the corresponding outputs provided by `G`. At present, this function 
  # assumes that the inverse problem has a Gaussian likelihood with covariance 
  # `Sig`, but this can potentially be generalized in the future. To be explicit,
  # the current assumption is an additive Gaussian noise model, with independent
  # noise eps ~ N(0,Sig). Note that no parameter transformations are performed
  # here; see `run_eki()` for such transformations.
  # 
  # Args:
  #    U: matrix of shape (J,D), where J = number of ensemble members and 
  #       D = dimension of each ensemble member.
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on.
  #    G: matrix of shape (J,P), storing the results of the forward model runs 
  #       evaluated at the ensemble members `U`. 
  #    Sig: matrix of shape (P,P), the covariance matrix for the Gaussian 
  #         likelihood. 
  #
  # Returns:
  #  list, with elements:
  #     `U`: the updated "u" ensemble, stored in a (J,D) matrix.
  #     `G`: the argument `G`.
  #     `m_u`: the estimated mean for the "u" part of the joint Gaussian 
  #            approximation.
  #     `m_y`: the estimated mean for the "y" part. 
  #     `C_u`: the estimated covariance for the "u" part. 
  #     `C_y`: the estimated covariance matrix for the "y" part.
  #     `L_y`: the lower Cholesky factor of `C_y`. 
  #     `C_uy`: the estimated cross-covariance matrix.
  #
  # The means and covariances define the joint Gaussian approximation implicit 
  # in the EnKF update.
  
  # Estimate means.
  m_u <- colMeans(U)
  m_y <- colMeans(G)
  
  # Estimate covariances.
  C_u <- cov(U)
  C_y <- cov(G) + Sig
  C_uy <- cov(U,G)
  L_y <- t(chol(C_y))
  
  # Generate "simulated observations".
  P <- ncol(Sig)
  J <- nrow(U)
  eps <- matrix(rnorm(P*J),nrow=J,ncol=P) %*% chol(Sig)
  Y <- G + eps
  
  # Compute EnKF update.
  U_updated <- compute_enkf_update(U, y, Y, C_uy, L_y=L_y)
  
  # Return list.
  list(U=U_updated, G=G, m_u=m_u, m_y=m_y, C_u=C_u, C_y=C_y,
       L_y=L_y, C_uy=C_uy)
}


compute_enkf_update <- function(U, y, Y, C_uy, C_y=NULL, L_y=NULL) {
  # Applies the standard EnKF update (analysis step) to an ensemble of particles 
  # `U`. No sample means or covariances are estimated in this function; this 
  # function simply evaluates the update equation given the requisite covariances
  # as arguments. 
  #
  # Args:
  #    U: matrix of shape (J,D), where J = number of ensemble members and 
  #       D = dimension of each ensemble member.
  #    y: numeric vector of length `P`, where `P` is the observation dimension. 
  #       This vector is the observation being conditioned on. 
  #    Y: matrix of shape (J,P), the ensemble of "simulated observations" that 
  #       are subtracted from the true observation `y` in the update equation.
  #    C_uy: matrix of shape (D,P), the cross-covariance matrix used in the 
  #          update equation. 
  #    C_y: matrix of shape (P,P), the "y" part of the covariance. Optional if 
  #         its Cholesky factor `L_y` is provided. 
  #    L_y: matrix of shape (P,P), the lower Cholesky factor of `C_y`. If NULL, 
  #         will be computed using `C_y`. 
  #
  # Returns:
  #    matrix of shape (J,D), the updated ensemble. 
  
  assert_that(!(is.null(C_y) && is.null(L_y)))
  
  # Compute lower Cholesky factor.
  if(is.null(L_y)) L_y <- t(chol(C_y))
  
  # Update ensemble.
  Y_diff <- add_vec_to_mat_rows(y,-Y)
  U + t(C_uy %*% backsolve(t(L_y), forwardsolve(L_y, t(Y_diff))))
}


calc_cond_Gaussian_moments <- function(y_obs, m_u, m_y, C_u, C_uy, C_y=NULL, L_y=NULL) {
  # Given the means and covariances to a joint Gaussian vector (u,y), returns
  # the mean and covariance of the conditional Gaussian u | [y=y_obs]. The 
  # covariance of y can be provided via `C_y`, or alternatively the lower
  # Cholesky factor of the covariance can be given by `L_y`.
  
  if(is.null(C_y) && is.null(L_y)) {
    stop("Either `C_y` or `L_y` required to compute conditional moments.")
  }

  y_obs <- drop(y_obs)
  m_y <- drop(m_y)
  
  # Compute Cholesky factor of cov[y].
  if(is.null(L_y)) L_y <- t(chol(C_y))
  
  # Compute conditional moments.
  L_inv_C_yu <- forwardsolve(L_y, t(C_uy))
  m_cond <- m_u + C_uy %*% backsolve(t(L_y), forwardsolve(L_y, y_obs-m_y))
  C_cond <- C_u - crossprod(L_inv_C_yu)
  
  list(mean=m_cond, cov=C_cond)
}


get_Gaussian_ldens_approx <- function(eki_results, y_obs, itr=NULL, normalize=TRUE,
                                      transform=FALSE) {
  # A convenience function that constructs the conditional log Gaussian density  
  # implied by the joint Gaussian ansatz at iteration `itr` of an EKI run.
  # This is a post-processing function that takes in the list returned by
  # `run_eki_pecan()` as an argument. By default constructs the density
  # for the final iteration. By default returns a normalized (log) Gaussian
  # density, representing an approximation to log p(u|y=y_obs). 
  # If `normalize = FALSE` then subtracts off the normalizing constant p(y_obs), 
  # which returns an unnormalized density that can be viewed as an approximation 
  # to log p(u,y_obs). By default, this density will be defined with respect
  # to the transformed/unconstrained space (i.e., density over phi). If
  # `transform = TRUE`, the density will be transformed back to the original
  # input parameter, so in this case the returned density may not be Gaussian.
  #
  # Note that this density will be defined with respect to the 
  # transformed/unconstrained space.
  #
  # TODO: currently the EKI results don't contain the observed data. This should
  # be updated, in which case the observations will not need to be provided
  # as an argument here. Storing the obs in the EkI list is especially important
  # when observations have been batched over iterations.
  
  if(is.null(itr)) itr <- length(eki_results$eki_list)
  cond_moments <- get_Gaussian_cond_moments(eki_results, y_obs, itr=itr) 
  
  # If transforming, need to add log det[DT(u)] to the density, where 
  # phi = T(u) is the parameter transformation. This is the same as 
  # subtracting log det[D T_inv(phi)], which is the quantity computed here.
  if(transform) {
    change_vars <- function(U) {
      U <- par_maps$inv(par_maps$fwd(U))
      -1 * attr(U, "log_det_J")
    }
  } else {
    change_vars <- function(...) 0
  }
  
  # If not normalizing, need to shift by p(y_obs) to obtain joint, indead
  # of conditional, density.
  if(normalize) {
    prob_y_obs <- 0
  } else {
    moments <- eki_results$eki_list[[itr]]
    prob_y_obs <- mvtnorm::dmvnorm(drop(inv_prob$y), mean=moments$m_y, 
                                   sigma=moments$C_y, log=TRUE)
  }
  
  # Construct the log-density function.
  if(transform) {
    
    cond_dens <- function(input) {
      # The input to this function is "par", the un-transformed parameter.
      mvtnorm::dmvnorm(par_maps$fwd(input), 
                       mean=drop(cond_moments$mean), 
                       sigma=cond_moments$cov, log=TRUE) + prob_y_obs + change_vars(input)
    }
    
  } else {
    
    cond_dens <- function(input) {
      # The input to this function is "phi", the transformed parameter.
      mvtnorm::dmvnorm(input, mean=drop(cond_moments$mean), 
                       sigma=cond_moments$cov, log=TRUE) + prob_y_obs + change_vars(input)
    }
    
  }

  return(cond_dens)
}


get_Gaussian_cond_moments <- function(eki_results, y_obs, itr=NULL, 
                                      reverse_conditional=FALSE) {
  # Same as `get_Gaussian_ldens_approx` but returns the conditional mean and
  # covariance, rather than the density defined by these quantities. If
  # `reverse_conditional = TRUE`, then the moments for y|u are returned as
  # opposed to u|y.
  #
  # Note that this density will be defined with respect to the 
  # transformed/unconstrained space.
  #
  # TODO: currently the EKI results don't contain the observed data. This should
  # be updated, in which case the observations will not need to be provided
  # as an argument here. Storing the obs in the EKI list is especially important
  # when observations have been batched over iterations.
  #
  # TODO: add support for `reverse_conditional` option, which returns 
  # the moments for y|u instead of u|y.
  
  if(reverse_conditional) {
    .NotYetImplemented()
  }
  
  if(is.null(itr)) itr <- length(eki_results$eki_list)
  
  info <- eki_results$eki_list[[itr]]
  
  cond_moments <- calc_cond_Gaussian_moments(y_obs, info$m_u, info$m_y, info$C_u, 
                                             info$C_uy, info$C_y, info$L_y)
  
  return(cond_moments)
}








