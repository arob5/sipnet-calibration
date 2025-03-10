#
# run_exact_mcmc.r
# Initial test of the brute-force MCMC functionality for PEcAn. This script
# sets up a statistical model for calibrating a subset of SIPNET's parameters
# and initial conditions, constrained by annual averages of state variables.
#
# Andrew Roberts
# 
# Working directory: calibrating_ics
#

library(PEcAn.settings)
library(PEcAn.workflow)
library(PEcAn.logger)
library(PEcAn.utils)
library(PEcAn.remote)
library(PEcAn.uncertainty)
library(dplyr)
library(ggplot2)
library(data.table)
library(assertthat)
library(lubridate)

# ------------------------------------------------------------------------------
# User-specified settings/filepaths.  
# ------------------------------------------------------------------------------

# Directories. `base_out_dir` will be used to overwrite the settings `outdir`,
# and will be used at the base output directory for the MCMC run.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration")
pecan_base_dir <- file.path(base_dir, "runs", "calibrating_ics")
src_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration", "src")
pecan_src_dir <- file.path(base_dir, "src")
base_out_dir <- file.path(pecan_base_dir, "mcmc")

# PEcAn XML path
pecan_settings_path <- file.path(pecan_base_dir, "sipnet_calibration.xml")

# Paths to constraint data, used to define observation vector "y" used in 
# calibration.
obs_mean_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites",
                           "Obs", "Rdata", "obs.mean.Rdata")
obs_cov_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites", 
                          "Obs", "Rdata", "obs.cov.Rdata")

# Sourcing helper functions.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "basis_function_emulation.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(pecan_src_dir, "param_calibration_functions.r"))
source(file.path(pecan_src_dir, "mcmc_pecan.r"))
source(file.path(pecan_src_dir, "prob_dists.r"))

# ------------------------------------------------------------------------------
# General setup.  
# ------------------------------------------------------------------------------

# Read settings. 
settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path)

# Update output directories.
settings$outdir <- base_out_dir
settings$modeloutdir <- file.path(base_out_dir, "out")
settings$rundir <- file.path(base_out_dir, "run")

# Check/validate settings. 
settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)

# Create directories, if not already created. 
dir.create(settings$outdir)
dir.create(settings$rundir)
dir.create(settings$modeloutdir)

# Write validated settings. 
PEcAn.settings::write.settings(settings, outputfile="pecan.CHECKED.xml")

# Specify the number of MCMC iterations.
n_mcmc_itr <- 20L

# Set so that all model runs are run on the same job, no new qsub jobs are
# dispatched.
settings$host <- list(name="localhost")


# ------------------------------------------------------------------------------
# Specify calibration parameters and their prior distributions.
# ------------------------------------------------------------------------------

prior_list <- list()

# Dirichlet prior on allocation parameters with sum-to-one constraint.
# Note that `root_allocation_fraction` maps to the SIPNET 
# fine root allocation parameter. The parameter `coarse_root_allocation_fraction`
# does not actually exist; it is implicitly defined as 1 minus the sum of the 
# other allocation parameters. We introduce it here for the purpose of 
# defining a Dirichlet prior on the allocation parameters. It will be
# removed before being passed to the write config function.
prior_list$alloc_pars <- list(
  param_name = "alloc_pars",
  dist_name = "Dirichlet",
  constraint = "simplex",
  length = 4L,
  scalar_names = c("root_allocation_fraction", "wood_allocation_fraction",
                   "leaf_allocation_fraction", "coarse_root_allocation_fraction"),
  dist_params = list(alpha=c(0.1, 0.5, 0.2, 0.2))
)

# Fractions controlling initial conditions for various wood pools. Also 
# subject to sum-to-one constraint.
prior_list$wood_frac_pars <- list(
  param_name = "wood_frac_pars",
  dist_name = "Dirichlet",
  constraint = "simplex",
  length = 3L,
  scalar_names = c("abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac"),
  dist_params = list(alpha=c(0.8, 0.1, 0.1))
)

# Aboveground wood initial condition. Note that this represents the total wood
# IC (which is not explicitly specified) times `abvGrndWoodFrac`.
prior_list$AbvGrndWood <- list(
  param_name = "AbvGrndWood",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length = 1L,
  dist_params = list(shape=10, rate=6.25e-05)
)

# Soil IC.
prior_list$soil <- list(
  param_name = "soil",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length = 1L,
  dist_params = list(shape=10, rate=5e-04)
)

# Root turnover rate constrained to (0,1).
prior_list$root_turnover_rate <- list(
  param_name = "root_turnover_rate",
  dist_name = "Beta",
  constraint = c(0,1),
  length = 1L,
  dist_params = list(shape1=2, shape2=4)
)

# Amax constrained to be nonnegative.
prior_list$Amax <- list(
  param_name = "Amax",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length = 1L,
  dist_params = list(shape=224, rate=2)
)


# ------------------------------------------------------------------------------
# Sample parameter values from prior to initialize MCMC settings.
# ------------------------------------------------------------------------------

# Number of design points to sample from prior.
n_prior_samp <- 10000L

# Sample from prior, and convert to unconstrained space. We use "phi" when 
# referring to the transformed/unconstrained parameters.
lprior <- get_sampling_func(prior_list)
par_maps <- get_par_map_funcs(prior_list)
par_samp <- lprior(n_prior_samp)
phi_samp <- par_maps$fwd(par_samp)

# Estimate prior mean and covariance for use in initializing MCMC.
prior_mean <- colMeans(phi_samp)
prior_cov <- cov(phi_samp)

# ------------------------------------------------------------------------------
# Set up likelihood function.
# ------------------------------------------------------------------------------

load(obs_mean_path)
load(obs_cov_path)

constraint_vars <- c("LAI", "AbvGrndWood")

# Specify single site.
site_id <- settings$run$site$id

# Data vector, determines observation order convention.
y_obs <- create_y_vec_annual(obs.mean, site_id=site_id, 
                             output_vars=constraint_vars, 
                             order_first_by_year=TRUE)
obs_order <- names(y_obs)
dim_obs <- length(drop(y_obs))

# Observation covariance matrix.
Sig <- create_y_noise_cov_annual(obs.cov, obs_mean_list=obs.mean, 
                                 site_id=site_id, obs_order=obs_order)

# Assuming Gaussian likelihood with diagonal covariance.
sig2 <- diag(Sig)

# Create a function that maps from the observable quantity to the 
# log-likelihood. We're treating the observation covariance as fixed here, 
# so the log-det normalizing constant term is excluded.
obsvb_to_llik_map <- function(obsvb) {
  -0.5 * sum((y_obs-obsvb)^2/sig2)
}

# Create observation operator, defined as a map from forward model outputs
# (states) to observable. In this context, we define it as a map from the run 
# ID to the observable. 

obs_op <- function(run_id, settings) {
  
  # Read predicted model trajectories of constraint variables.
  run_id_path <- file.path(settings$modeloutdir, run_id)
  model_output <- PEcAn.utils::read.output(run_id, outdir=run_id_path, 
                                           variables=constraint_vars, 
                                           dataframe=TRUE)
  model_output <- as.data.table(model_output)
  model_output[, posix := NULL]
  model_output[, year := as.character(year)]
  
  # Compute annual averages.
  model_output <- melt.data.table(model_output, id.vars="year", 
                                  variable.name="output_var", 
                                  variable.factor=FALSE)
  model_output <- model_output[, .(value=mean(value)), by=.(output_var, year)]
  
  # Extract predicted observables in the correct order. Note that `sort=FALSE`
  # ensures that the order of `dt_order` is maintained.
  dt_order <- get_dt_obs_order_from_vec(obs_order)
  obsvb <- merge.data.table(dt_order, model_output, by=c("output_var", "year"),
                            all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  return(obsvb$value)
}


prepare_mcmc_model_run <- function(settings, par, run_id) {
  # `par` must be a names attribute set.
  
  # Extract model initial conditions. Write config function expects this 
  # to be a named list.
  ic <- as.list(par[ic_names])
  
  # Extract model parameters. Write config function expects `par` to be a 
  # list of data.frames, with column names set to parameter names.
  par <- par[setdiff(param_names, exclude_par)]
  par <- as.data.frame.list(par)
  
  # Create run and output directories.
  dir.create(file.path(settings$rundir, run_id), recursive=TRUE)
  dir.create(file.path(settings$modeloutdir, run_id), recursive=TRUE)
  
  # Write config to file.
  model_write_config <- paste("write.config.", settings$model$type, sep = "")
  do.call(model_write_config, args=list(defaults = settings$pfts,
                                        trait.values = par,
                                        IC = ic,
                                        settings = settings,
                                        run.id = run_id))
  
  # Write the run ID to a new line in the runs text file. This will override 
  # any existing "runs.txt" file.
  cat(as.character(run_id),
      file = file.path(settings$rundir, "runs.txt"),
      sep = "\n",
      append=FALSE)
}


# Log-likelihood function.
llik <- function(par, run_id) {
  # Set up paths and write configs.
  prepare_mcmc_model_run(settings, par, run_id)
  
  # Run model at parameter `par`.
  PEcAn.workflow::start_model_runs(settings, write=FALSE)
  
  # Compute the model predicted observable quantity.
  obs_pred <- obs_op(run_id, settings)
  
  # Evaluate log-likelihood.
  obsvb_to_llik_map(obs_pred)
}


# ------------------------------------------------------------------------------
# Run MCMC
# ------------------------------------------------------------------------------

d <- length(all_par_names)
cov_prop_init <- prior_cov + diag(1e-8, nrow=d, ncol=d)

mcmc_output <- mcmc_mh(llik, lprior, par_names=all_par_names, n_itr=n_mcmc_itr, 
                       par_init=prior_mean, prop_settings=list(cov_prop=cov_prop_init))









