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


create_y_vec_annual <- function(obs_mean_list, site_id, output_vars=NULL,
                                order_first_by_year=TRUE) {
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
  # obs_mean_list: list, each element corresponding to a year. The names of the 
  #                list are currently assumed to be in ymd format (ultimately
  #                only the year is extracted). Each element is itself a list
  #                over sites. Names of this sublist must correspond to the
  #                site IDs. Elements of the sublist are single-row data.frames 
  #                with column names set to output variables.
  # site_id: the single site ID used to extract the relevant site from 
  #          `obs_mean_list`.
  # output_vars: character, optional vector of output variables to include. If 
  #              NULL, uses all output variables found in `obs_mean_list`.
  # order_first_by_year: logical, if TRUE orders the observations first by
  #                      year (ascending) then output variable 
  #                      (alphabetically ascending). If FALSE, orders by 
  #                      output variable then year.
  #
  # Returns:
  # vector, with names of the form `<output_var>_<year>`, containing data
  # extracted from `obs_mean_list` for the specified site. No NA values are 
  # included in the returned vector. The order of the vector is determined by
  # the `order_first_by_year` argument.
  
  # Assemble outputs in table.
  dt_obs <- create_y_dt_annual(obs_mean_list, site_id, output_vars=output_vars)
  output_vars <- setdiff(colnames(dt_obs), "year")
  
  # Create "y" vector.
  dt_obs_long <- melt.data.table(dt_obs, id.vars="year", 
                                 variable.name="output_var", na.rm=TRUE,
                                 variable.factor=FALSE)
  
  # Set observation order.
  if(order_first_by_year) setorder(dt_obs_long, year, output_var)
  else setorder(dt_obs_long, output_var, year)
  
  # Construct observation vector.
  y_vec <- dt_obs_long$value
  names(y_vec) <- paste(dt_obs_long$output_var, dt_obs_long$year, sep="_")  
  
  return(y_vec)
}


create_y_dt_annual <- function(obs_mean_list, site_id, output_vars=NULL) {
  # This function takes the "obs.mean" list in the format required by the SDA
  # code (see details below), which is assumed to store annual observations of 
  # various output variables for various sites. This function extracts the 
  # available observations for the single site with site ID `site_id` and 
  # assembles them into a table (see `Returns` below). See also
  # `create_y_vec_annual()`, which takes the table returned by this function
  # and converts it into a vector.
  #
  # Args:
  # obs_mean_list: list, each element corresponding to a year. The names of the 
  #                list are currently assumed to be in ymd format (ultimatlely
  #                only the year is extracted). Each element is itself a list
  #                over sites. Names of this sublist must correspond to the
  #                site IDs. Elements of the sublist are single-row data.frames 
  #                with column names set to output variables.
  # site_id: the single site ID used to extract the relevant site from 
  #          `obs_mean_list`.
  # output_vars: character, optional vector of output variables to include. If 
  #              NULL, uses all output variables found in `obs_mean_list`.
  #
  # Returns:
  # data.table, with one row per year. Column names include "year" as well 
  # as all output variables found in "obs_mean_list". NA values indicate that
  # the site did not have data for a certain output variable for a certain 
  # year.

  # Extract only data for the specified site. Create data.table with columns
  # "", "year". Missing data results in NAs.
  dt_obs_list <- lapply(seq_along(obs_mean_list), extract_site_obs, 
                        site_id=site_id, obs_mean_list=obs_mean_list)
  dt_obs <- rbindlist(dt_obs_list, use.names=TRUE, fill=TRUE)
  dt_obs[, posix := lubridate::ymd_hms(paste(posix, "23:59:59"))]
  dt_obs[, year := lubridate::year(posix)]
  dt_obs[, posix := NULL]
  setorder(dt_obs, year)
  output_vars_dt <- setdiff(colnames(dt_obs), "year")
  
  # Restrict to specified output variables.
  if(!is.null(output_vars)) {
    missing_vars <- setdiff(output_vars, output_vars_dt)
    output_vars <- intersect(output_vars, output_vars_dt)
    if(length(output_vars)==0) {
      stop("No output variables in `output_vars` found in `obs_mean_list`.")
    }
    
    if(length(missing_vars) > 0) {
      message("Specified variables not found in `obs_mean_list`: ",
              paste(missing_vars, sep=", "))
    }
    
    keep_cols <- c(output_vars, "year")
    dt_obs <- dt_obs[,..keep_cols]
  } else {
    output_vars <- output_vars_dt
  }
  
  # Order columns.
  setcolorder(dt_obs, output_vars)
  
  return(dt_obs)
}


get_y_noise_cov_annual <- function(obs_cov_list, obs_mean_list, site_id, 
                                   obs_order) {
  # This function constructs an observation covariance matrix for a single 
  # site based on the values found in the "obs.cov" list in the format required
  # by the SDA code (see details below). The observations must be at an annual
  # time step. The order of the covariance matrix is determined by `obs_order`
  # which is a character vector with entries of the form `<output_var>_<year>`;
  # e.g., "LAI_2017" (the convention used to set the column names in the 
  # vector returned by `create_y_vec_annual()`).
  #
  # Args:
  # obs_cov_list: list, each element corresponding to a year. The names of the 
  #               list are currently assumed to be in ymd format (ultimatlely
  #               only the year is extracted). Each element is itself a list
  #               over sites. Names of this sublist must correspond to the
  #               site IDs. Elements of the sublist are covariance matrices 
  #               across output variables for that specific year/site. The
  #               ordering of this covariance matrix is determined by the 
  #               order of the variables in the corresponding element of 
  #               "obs_mean_list".
  # obs_mean_list: list, see `create_y_vec_annual()` for a more detailed 
  #                description. Only required here in order to extract the 
  #                variable names/order in order to interpret the covariances
  #                in `obs_cov_list`.
  # site_id: the single site ID used to extract the relevant site from 
  #          `obs_mean_list` and `obs_cov_list`.
  # obs_order: character, vector with elements of the form `<output_var>_<year>`
  #            used to determine the ordering of the returned covariance matrix.
  #            See details below.
  #
  # Returns:
  # Covariance matrix of dimension `length(obs_order)`. Row and col names will 
  # be set to `obs_order`. The covariance matrix stores covariances across 
  # years and output variables for the single site with site ID `site_id`.

  # Extract the output variable and year.
  dt_order <- as.data.table(data.table::tstrsplit(obs_order, split="_", 
                                                  fixed=TRUE))
  colnames(dt_order) <- c("output_var", "year")
  
  # Convert names of "obs.mean" and "obs.cov" lists to years.
  names(obs_mean_list) <- year(lubridate::ymd(names(obs_mean_list)))
  names(obs_cov_list) <- year(lubridate::ymd(names(obs_cov_list)))

  # Empty covariance matrix to be filled in.
  obs_dim <- length(obs_order)
  cov_mat <- matrix(nrow=obs_dim, ncol=obs_dim, 
                    dimnames=list(obs_order, obs_order))
  
  # Fill in covariance matrix year by year.
  years <- unique(dt_order$year)
  for(yr in years) {
    
    # The "obs.mean" list is required as it determines the order of the
    # variables in "obs.cov".
    output_vars <- names(obs_mean_list[[yr]][[site_id]])
    
    # Covariance matrix across output variables for the year.
    cov_yr <- as.matrix(obs_cov_list[[yr]][[site_id]])
    colnames(cov_yr) <- rownames(cov_yr) <- output_vars
    
    # Select the observation indices in the correct order. This may extract
    # a subset of the variables in `output_vars`.
    idcs <- dt_order[year==yr, which=TRUE]
    vars_include <- dt_order[idcs, output_var]
    
    # Fill in covariances.
    cov_mat[idcs, idcs] <- cov_yr[vars_include, vars_include]
  }

  return(cov_mat)
}


extract_site_obs <- function(time_idx, site_id, obs_mean_list) {
  # Helper function for `create_y_obs_annual()`. Extracts the data.frame
  # of output variables from `obs_mean_list` for the specified time index
  # and site. Converts to a data.table and adds a "posix" column storing
  # the date. Note that `time_idx` here is simply the integer index for the 
  # top level of the list `obs_mean_list`.
  
  site_dt <- as.data.table(obs_mean_list[[time_idx]][[site_id]])
  site_dt[, posix := names(obs_mean_list)[time_idx]]
  return(site_dt)
}









