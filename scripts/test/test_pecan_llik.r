# test_pecan_llik.r
#
# Testing the ability to generate a log-likelihood function that depends on 
# a PEcAn model run. The idea is to have a function with an interface like
# a typical log-likelihood function, which hides the complexity of the 
# underlying model run from the end user. This script assumes that a PEcAn
# XML object has already been set up, and also assumes access to some 
# data constraints to use in defining the log-likelihood function. 
#
# Andrew Roberts
#

library(dplyr)
library(ggplot2)
library(lubridate)
library(data.table)
library(assertthat)
library(PEcAn.settings)

# Directories. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration")
pecan_base_dir <- file.path(base_dir, "runs", "sipnet_calibration_2024")
calibration_dir <- file.path(pecan_base_dir, "calibration")
src_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration", "src")
pecan_src_dir <- file.path(base_dir, "src")

# Sourcing helper functions.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "basis_function_emulation.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(pecan_src_dir, "param_calibration_functions.r"))
source(file.path(pecan_src_dir, "mcmc_pecan.r"))

# Paths to files.
obs_mean_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites",
                           "Obs", "Rdata", "obs.mean.Rdata")
obs_cov_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites", 
                          "Obs", "Rdata", "obs.cov.Rdata")
param_path_plant_pft <- file.path(pecan_base_dir, "param", "param_design", 
                                  "param_design_plant_pft.csv")
param_path_soil_pft <- file.path(pecan_base_dir, "param", "param_design", 
                                 "param_design_soil_pft.csv")
pecan_settings_path <- file.path(pecan_base_dir, "sipnet_calibration.xml")
# obs_path <- file.path(calibration_dir, "obs_forecast_aligned.csv")


# ------------------------------------------------------------------------------
# Setup.
# ------------------------------------------------------------------------------

# Read PEcAn settings object.
settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path)

# Specify single site.
site_id <- settings$run$site$id

# Load SIPNET package. 
model_write_config <- paste("write.config.", settings$model$type, sep = "")
PEcAn.utils::load.modelpkg(settings$model$type)

# Read parameter settings.
param_design_plant_pft <- fread(param_path_plant_pft)
param_design_soil_pft <- fread(param_path_soil_pft)

# Store in a list over PFTs. 
param_design_list_by_pft <- list(param_design_plant_pft, param_design_soil_pft) 

# ------------------------------------------------------------------------------
# Set up data vector and likelihood information. 
# ------------------------------------------------------------------------------

load(obs_mean_path)
load(obs_cov_path)

constraint_vars <- c("LAI", "AbvGrndWood")

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

# Generate log-likelihood function.
llik <- generate_pecan_llik(settings, obs_op, obsvb_to_llik_map)


# ------------------------------------------------------------------------------
# Evaluate llik at parameter setting.
# ------------------------------------------------------------------------------

run_id <- "test_run"
par_idx <- 1L
par <- lapply(param_design_list_by_pft, 
              function(pft_params) as.data.frame(pft_params[par_idx,,drop=FALSE]))


# TODO: add settings as another argument?
llik_eval_test <- llik(par, run_id)






