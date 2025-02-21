#
# prepare_obs_constraints.r
# Formatting the constraint data to be used for the parameter calibration. The 
# goal of the parameter calibration is to characterize the posterior 
# p(par|Y) where `par` is the calibration parameters. By "constraint data"
# we are referring to "Y". 
#
# Andrew Roberts
#

# # pecan standard unit: kg C m-2
# model_agw <- model_output[, .(posix,year,run_id,AbvGrndWood)]
# model_agw_mean <- model_output[, .(agw_mean=mean(AbvGrndWood)), by=.(year,run_id)]
# obs_agw <- dt_obs[, .(posix, AbvGrndWood)]
# obs_agw[, year := lubridate::year(posix)]
# 
# model_agw_plt <- ggplot(model_agw_mean, aes(x=as.factor(year),y=agw_mean)) + 
#   geom_boxplot() + geom_point(aes(x=as.factor(year),y=AbvGrndWood/10), obs_abw, color="red")
# plot(model_agw_plt)

library(PEcAn.settings)
library(PEcAn.workflow)
library(PEcAn.logger)
library(PEcAn.utils)
library(PEcAn.remote)
library(PEcAn.uncertainty)
library(PEcAn.benchmark)
library(dplyr)
library(ggplot2)
library(data.table)
library(assertthat)
library(lubridate)

# ------------------------------------------------------------------------------
# User-specified settings/filepaths.  
# ------------------------------------------------------------------------------

# General file paths. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration", 
                      "runs", "sipnet_calibration_2024")
calibration_dir <- file.path(base_dir, "calibration")
local_edits_path <- file.path(base_dir, "..", "..", "local_edits", "local_edits.r")
pecan_settings_path <- file.path(base_dir, "sipnet_calibration.xml")

# Observation/constraint data paths.
obs_mean_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites",
                           "Obs", "Rdata", "obs.mean.Rdata")
obs_cov_path <- file.path("/projectnb", "dietzelab", "dongchen", "anchorSites", 
                          "Obs", "Rdata", "obs.cov.Rdata")
obs_out_path <- file.path(calibration_dir, "obs_forecast_aligned.csv")

source(local_edits_path)

# ------------------------------------------------------------------------------
# General setup.  
# ------------------------------------------------------------------------------

# Read settings. 
settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path)

# Check/validate settings. 
settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)

# Create directories, if not already created. 
dir.create(settings$outdir)
dir.create(settings$rundir)
dir.create(settings$modeloutdir)
dir.create(calibration_dir)

# ------------------------------------------------------------------------------
# Read Observation Data.   
# ------------------------------------------------------------------------------

# These are lists over time periods (with names attributes set to the time 
# stamps). Each element is itself a list over spatial locations (with names 
# set to site ID). Each element of these sub-lists in `obs.mean` is a data.frame 
# with 1 row and colnames "AbvGrndWood", "LAI", and "TotSoilCarb". The respective 
# elements in `obs.cov` are 3-by-3 matrices (the covariance matrix). 
load(obs_mean_path)
load(obs_cov_path)

# Restrict to the single spatial location. 
site_id <- settings$run$site$id
create_site_dt <- function(time_idx) {
  site_dt <- as.data.table(obs.mean[[time_idx]][[site_id]])
  site_dt[, posix := names(obs.mean)[time_idx]]
  return(site_dt)
}

dt_obs_list <- lapply(seq_along(obs.mean), create_site_dt)
dt_obs <- rbindlist(dt_obs_list, use.names=TRUE, fill=TRUE)
dt_obs[, posix := lubridate::ymd_hms(paste(posix, "23:59:59"))]
dt_obs[, year := lubridate::year(posix)]
dt_obs[, posix := NULL]

# ------------------------------------------------------------------------------
# Read model output.     
# ------------------------------------------------------------------------------

# Output variables to read. 
output_vars <- c("time", "NEE", "LAI", "GPP", "TotSoilCarb", "AbvGrndWood", 
                 "leaf_carbon_content", "fine_root_carbon_content", "AGB", 
                 "coarse_root_carbon_content", "wood_carbon_content", 
                 "litter_carbon_content")
output_vars <- union(output_vars, setdiff(colnames(dt_obs), "datetime"))

# Used to specify path to read model outputs from. 
run_ids <- readLines(file.path(settings$rundir, "runs.txt"))

# List of matrices storing model output, one element per ensemble run (i.e., per run ID). 
output_list <- list()

for(run_id in run_ids) {
  run_id_path <- file.path(settings$modeloutdir, run_id)
  model_run_output <- PEcAn.utils::read.output(run_id, outdir=run_id_path, 
                                               variables=output_vars, dataframe=TRUE)
  model_run_output <- as.data.table(model_run_output)
  model_run_output[, run_id := run_id]
  output_list[[run_id]] <- model_run_output
}

# Combine into single data.table. 
model_output <- rbindlist(output_list, use.names=TRUE)
fwrite(model_output, file.path(calibration_dir, "model_forecasts_design.csv"))

# ------------------------------------------------------------------------------
# Format model output to align temporally with observation data.    
# ------------------------------------------------------------------------------

# Compute number of seconds elapsed relative to the earliest timestamp. 
# model_output[, seconds := PEcAn.utils::ud_convert(time, "days" ,"seconds")]

# Temporal alignment. 
# Previously using the PEcAn helper, but in this case I'm simply taking 
# annual averages, so now changing to doing this manually. 
# common_colnames <- intersect(colnames(model_output), colnames(dt_obs))
# run_ids <- unique(model_output$run_id)
# align_model_run <- function(run_id) {
#   model_run_aligned <- PEcAn.benchmark::align_data(model.calc=model_output[run_id==run_id, ..common_colnames], 
#                                                    obvs.calc=dt_obs[, ..common_colnames], 
#                                                    var=setdiff(common_colnames, "posix"), 
#                                                    align_method="mean_over_larger_timestep")
#   model_run_aligned[, run_id := run_id]
#   return(model_run_aligned)
# }
# 
# dt_aligned_list <- lapply(run_ids, align_model_run)
# dt_aligned <- rbindlist(dt_aligned_list, use.names=TRUE)

# Manually computing annual averages of model outputs.
model_output_annual <- model_output[, .(year, run_id, LAI, TotSoilCarb, AbvGrndWood)]
model_output_annual <- model_output_annual[, .(LAI=mean(LAI), 
                                               TotSoilCarb=mean(TotSoilCarb),
                                               AbvGrndWood=mean(AbvGrndWood)), by=.(year, run_id)]
model_output_annual[, type := "model"]

# Merge model outputs and observations.
dt_obs[, type := "obs"]
dt_obs[, run_id := NA_character_]
dt_obs[, SoilMoistFrac := NULL]

dt_aligned <- rbindlist(list(dt_obs, model_output_annual), use.names=TRUE)

fwrite(dt_aligned, obs_out_path)





