#
# param_calibration_functions.r
#
# This file contains functions that process PEcAn forward model outputs and 
# observations into the the quantities required for setting up the inverse 
# problem of calibrating model parameters. Some functions are aimed at setting 
# up the likelihood, and others assemble the model outputs into the "design
# matrix" format which can then be used to fit emulators of the forward model. 
# Crucial to all of these steps is the concept of the observation operator, 
# which maps the raw model output (trajectories of the state variables) to 
# the observation space. These functions are designed to handle model runs 
# at a single site. In the future, wrappers can be written to leverage these 
# functions in extending to multi-site calibration.
#
# These functions are currently geared towards the setting of calibrating with
# respect to constraint data that comes in the form of annual means. We will 
# want to think about generalizing these functions to handle more varied 
# situations in the future. 
#
# Andrew Roberts
#

library(data.table)


create_y_obs_annual <- function(obs_mean_list, site_id) {
  # This function takes the "obs.mean" list in the format required by the SDA
  # code (see details below), which is assumed to store annual observations of 
  # various output variables for various sites. This function extracts the 
  # available observations for the single site with site ID `site_id` and 
  # assembles them into a vector. This vector is intended to represent the "y"
  # vector in a typical inverse problem setup: y = G(u) + eps. Note that some 
  # output variables may be available for some years but not others. The 
  # order of the variables in "y" is critical for aligning with the model
  # outputs and for other uses. At present, the ordering is encoded via the 
  # names attribute of the vector. Names are defined as `<output_var>_<year>`;
  # e.g., "LAI_2017".
  #
  # Args:
  # obs_mean_list:
  # site_id: 
  #
  # Returns:
  #
  
  # Extract only data for the specified site. Create data.table with columns
  # "", "year". Missing data results in NAs.
  dt_obs_list <- lapply(seq_along(obs_mean_list), extract_site_obs, 
                        site_id=site_id)
  dt_obs <- rbindlist(dt_obs_list, use.names=TRUE, fill=TRUE)
  dt_obs[, posix := lubridate::ymd_hms(paste(posix, "23:59:59"))]
  dt_obs[, year := lubridate::year(posix)]
  dt_obs[, posix := NULL]
  
  return(dt_obs)
}


get_obs_cov_annual <- function(obs_cov_list, site_id, obs_order) {
  # This function constructs an observation covariance matrix for a single 
  # site based on the values found in the "obs.cov" list in the format required
  # by the SDA code (see details below). The observations must be at an annual
  # time step. The order of the covariance matrix is determined by `obs_order`
  # which is a charater vector with entries of the form `<output_var>_<year>`;
  # e.g., "LAI_2017".
  #
  #
  #
  #
  
}


extract_site_obs <- function(time_idx, site_id) {
  site_dt <- as.data.table(obs.mean[[time_idx]][[site_id]])
  site_dt[, posix := names(obs.mean)[time_idx]]
  return(site_dt)
}









