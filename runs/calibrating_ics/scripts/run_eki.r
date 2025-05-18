#
# run_eki.r
# Initial test of the Ensemble Kalman Inversion for PEcAn.
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
library(tictoc)
library(docopt)

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --n_eki_itr=<n_mcmc_itr>                  Number of EKI iterations (integer)
  --n_ens=<chain_idx>                       EKI ensemble size (integer)
" -> doc

# Read command line arguments.
cmd_args <- docopt(doc)
n_eki_itr <- as.integer(cmd_args$n_eki_itr)
n_ens <- as.integer(cmd_args$n_ens)

# ------------------------------------------------------------------------------
# User-specified settings/filepaths.  
# ------------------------------------------------------------------------------

# Directories. `chain_dir` will be used to overwrite the settings `outdir`,
# and will be used at the base output directory for the MCMC run.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration")
pecan_base_dir <- file.path(base_dir, "runs", "calibrating_ics")
src_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration", "src")
pecan_src_dir <- file.path(base_dir, "src")
eki_dir <- file.path(pecan_base_dir, "eki")
stat_model_dir <- file.path(pecan_base_dir, "statistical_model")

# PEcAn XML path
pecan_settings_path <- file.path(pecan_base_dir, "sipnet_calibration.xml")

# Paths to constraint data, used to define observation vector "y" used in 
# calibration.
obs_mean_path <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/Obs/obs.mean.Rdata"
obs_cov_path <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/Obs/obs.cov.Rdata"

# Sourcing helper functions.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "basis_function_emulation.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "ens_Kalman_inversion.r"))
source(file.path(pecan_src_dir, "param_calibration_functions.r"))
source(file.path(pecan_src_dir, "eki_pecan.r"))
source(file.path(pecan_src_dir, "prob_dists.r"))

# ------------------------------------------------------------------------------
# General setup.  
# ------------------------------------------------------------------------------

# Read settings. 
settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path)

# Update output directories. These will be appended to by the EKI code. Only 
# `outdir` needs to be specified here as the "base" output directory. The EKI 
# code will create sub-directories for each EKI iteration, setting `rundir`
# and `modeloutdir` as needed.
settings$outdir <- eki_dir

# Check/validate settings. 
settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)

# Create base output directory, if not already created. 
dir.create(settings$outdir, recursive=TRUE)

# Write validated settings. 
PEcAn.settings::write.settings(settings, outputfile="pecan.CHECKED.xml")

# Load SIPNET package.
model_write_config <- paste("write.config.", settings$model$type, sep = "")
PEcAn.utils::load.modelpkg(settings$model$type)

# ------------------------------------------------------------------------------
# Specify calibration parameters and their prior distributions.
# ------------------------------------------------------------------------------

prior_list <- list()
pecan_ic_names <- c("abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac", 
                    "AbvGrndWood", "soil")
pecan_par_names <- c("root_turnover_rate", "Amax", "root_allocation_fraction", 
                     "wood_allocation_fraction", "leaf_allocation_fraction")
ic_vs_par_list <- list(ic_names=pecan_ic_names, par_names=pecan_par_names)



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

saveRDS(prior_list, file.path(stat_model_dir, "prior_list.rds"))
saveRDS(ic_vs_par_list, file.path(stat_model_dir, "ic_vs_par_list.rds"))

# ------------------------------------------------------------------------------
# Set up forward model and likelihood function.
# ------------------------------------------------------------------------------

load(obs_mean_path)
load(obs_cov_path)

constraint_vars <- c("LAI", "AbvGrndWood", "TotSoilCarb")

# TODO: TEMP - manually setting this for now.
# Specify single site.
# site_id <- settings$run$site$id
site_id <- 3380

# Data vector, determines observation order convention. Soil carbon obs are 
# constant across years so we're only keeping the first one. Also cap 
# at 2021 for now.
y_obs <- create_y_vec_annual(obs.mean, site_id=site_id, 
                             output_vars=constraint_vars, 
                             order_first_by_year=TRUE)
drop_vars <- paste0("TotSoilCarb_", 2013:2021)
obs_order <- setdiff(names(y_obs), drop_vars)
dt_order <- get_dt_obs_order_from_vec(obs_order)
dt_order <- dt_order[year <= 2021]
obs_order <- paste(dt_order$output_var, dt_order$year, sep="_")
y_obs <- y_obs[obs_order]
dim_obs <- length(drop(y_obs))

# Observation covariance matrix.
Sig <- create_y_noise_cov_annual(obs.cov, obs_mean_list=obs.mean, 
                                 site_id=site_id, obs_order=obs_order)

# Create observation operator, defined as a map from forward model outputs
# (states) to observable. In this context, we define it as a map from the run 
# ID to the observable. 

obs_op <- function(run_ids, settings) {
  
  # List of matrices storing model output, one element per ensemble run (i.e., per run ID). 
  output_list <- list()
  
  for(run_id in run_ids) {
    run_id_path <- file.path(settings$modeloutdir, run_id)
    model_run_output <- PEcAn.utils::read.output(run_id, outdir=run_id_path, 
                                                 variables=constraint_vars, 
                                                 dataframe=TRUE)
    model_run_output <- as.data.table(model_run_output)
    model_run_output[, run_id := run_id]
    output_list[[run_id]] <- model_run_output
  }
  
  # Combine into single data.table. 
  model_output <- rbindlist(output_list, use.names=TRUE)
  
  # Compute annual averages.
  model_output[, posix := NULL]
  model_output[, year := as.character(year)]
  
  model_output <- melt.data.table(model_output, id.vars=c("year", "run_id"), 
                                  variable.name="output_var", 
                                  variable.factor=FALSE)
  model_output <- model_output[, .(value=mean(value)), by=.(output_var, year, run_id)]
  
  # Extract predicted observables in the correct order. Note that `sort=FALSE`
  # ensures that the order of `dt_order` is maintained. Format into matrix 
  # of shape (num runs, num obs).
  model_output[, obs_id := paste(output_var, year, sep="_")]
  model_output[, `:=`(output_var=NULL, year=NULL)]
  model_output <- data.table::dcast(data=model_output, formula=run_id~obs_id, 
                                    value.var="value")
  
  # Ensure rows are ordered by run ID and columns align with order dictated
  # by `obs_order`. Convert to matrix. Not actually checking row order here; 
  # I think it should be ordered properly but need to improve this later.
  as.matrix(model_output[, ..obs_order])
}


prepare_eki_model_run <- function(settings, par) {

  # Setup.
  model_write_config <- paste("write.config.", settings$model$type, sep = "")

  # `par` should be a matrix of shape (num runs, param dim), where param 
  # dimension includes both model parameters and initial conditions. colnames 
  # must be set. Extract ICs and model parameters as separate matrices here.
  n_ens <- nrow(par)
  ic <- par[,pecan_ic_names, drop=FALSE]
  par <- par[, pecan_par_names, drop=FALSE]
  
  # One run ID per ensemble member.
  run_ids <- paste0("ens_", 1:n_ens)
  
  # Write a model config for each ensemble member.
  for(i in seq_along(run_ids)) {
    # Create sub-directories for each run.
    run_id <- run_ids[i]
    dir.create(file.path(settings$rundir, run_id), recursive=TRUE)
    dir.create(file.path(settings$modeloutdir, run_id), recursive=TRUE)
    
    # Write config function expects IC to be a named list.
    ic_ens <- as.list(ic[i,])
    
    # Write config function expects `par` to be a list of data.frames, with 
    # column names set to parameter names.
    par_ens <- list(as.data.frame.list(par[i,]))
    
    # Write config to file.
    do.call(model_write_config, args=list(defaults = settings$pfts,
                                          trait.values = par_ens,
                                          IC = ic_ens,
                                          settings = settings,
                                          run.id = run_id))
    
    # Write the run ID to a new line in the runs text file. Append run ID to
    # "runs.txt" file.
    cat(run_id,
        file = file.path(settings$rundir, "runs.txt"),
        sep = "\n",
        append = (i != 1L))
  }
  
  return(run_ids)
}


# Forward model, i.e., param-to-observable map ("G" function). Vectorized so 
# that `par` can be a matrix with one row per parameter value.
fwd <- function(par, itr, settings) {

  # Set up paths and write configs. Each EKI itr gets its own modeloutdir and
  # rundir.
  settings <- update_settings_path_for_eki_itr(settings, itr)
  
  run_ids <- prepare_eki_model_run(settings, par)
  
  # Run model at parameter `par`.
  PEcAn.workflow::start_model_runs(settings, write=FALSE)
  
  # Compute the model predicted observable quantity.
  obs_op(run_ids, settings)
}


update_settings_path_for_eki_itr <- function(settings, itr) {
  itr_lbl <- paste0("itr_", itr)
  
  # Update paths. Each EKI itr gets its own base output directory.
  settings$outdir <- file.path(settings$outdir, itr_lbl)
  settings$rundir <- file.path(settings$outdir, "run")
  settings$modeloutdir <- file.path(settings$outdir, "out")
  settings$host$rundir <- settings$rundir 
  settings$host$outdir <- settings$modeloutdir
  settings$host$folder <- settings$modeloutdir
  
  # Create directories.
  dir.create(settings$rundir, recursive=TRUE)
  dir.create(settings$modeloutdir, recursive=TRUE)
  
  # Save updated settings. Note that these are written to `settings$outdir`.
  PEcAn.settings::write.settings(settings, outputfile="pecan.CHECKED.xml")
  
  return(settings)
}

lik_list <- list(y=y_obs, Sig=Sig, obs_op=obs_op, fwd=fwd,
                 prepare_eki_model_run=prepare_eki_model_run,
                 update_settings_path_for_eki_itr=update_settings_path_for_eki_itr)

saveRDS(lik_list, file.path(stat_model_dir, "lik_list.rds"))


# ------------------------------------------------------------------------------
# Run EKI
# ------------------------------------------------------------------------------

tic()
eki_output <- run_eki_pecan(y_obs, fwd, Sig, prior_list, 
                            n_itr=n_eki_itr, n_ens=n_ens, settings=settings)
print("EKI runtime:")
toc()

saveRDS(eki_output, file.path(eki_dir, "eki_output.rds"))


