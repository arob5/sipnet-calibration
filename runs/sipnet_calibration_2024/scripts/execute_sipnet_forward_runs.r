#
# execute_sipnet_forward_runs.r
# Forward evaluations of SIPNET using different parameter values. A single 
# initial condition and model driver trajectory are used. These are read from 
# file but note that the initial condition is manually adjusted in this script.
#
# Andrew Roberts
# 
# Working directory: sipnet_test_run
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

# General file paths. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration", 
                      "runs", "sipnet_calibration_2024")
local_edits_path <- file.path(base_dir, "..", "..", "local_edits", "local_edits.r")
helper_funcs_path <- file.path(base_dir, "..", "..", "local_edits", "helper_functions_temp.r")
pecan_settings_path <- file.path(base_dir, "sipnet_calibration.xml")

# Paths to files containing the parameter design at which to run the forward model. 
param_path_plant_pft <- file.path(base_dir, "param", "param_design", 
                                  "param_design_plant_pft.csv")
param_path_soil_pft <- file.path(base_dir, "param", "param_design", 
                                 "param_design_soil_pft.csv")

source(local_edits_path)
source(helper_funcs_path)

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

# Write validated settings. 
PEcAn.settings::write.settings(settings, outputfile="pecan.CHECKED.xml")

# ------------------------------------------------------------------------------
# Load parameter design. 
#     - Parameter input design created and saved by `prepare_model_parameters.r`
#       is just read here. 
# ------------------------------------------------------------------------------

# Number of design points. 
N_design <- as.integer(settings$ensemble$size)

# Read files. 
param_design_plant_pft <- fread(param_path_plant_pft)
param_design_soil_pft <- fread(param_path_soil_pft)

# Store in a list over PFTs. 
param_design_list_by_pft <- list(param_design_plant_pft, param_design_soil_pft) 

# ------------------------------------------------------------------------------
# Prepare for forward runs: write one config per design point. 
# ------------------------------------------------------------------------------

# `pda.init.run` requires these settings or will throw error. 
settings$assim.batch$method <- "emulator"  
settings$assim.batch$chain <- 1

# Load SIPNET package. 
model_write_config <- paste("write.config.", settings$model$type, sep = "")
PEcAn.utils::load.modelpkg(settings$model$type)

# Save write configs: replacing call to pda.init.run here to be more explicit
# about what is going on.

run_ids <- paste("run", 1:N_design, sep=".")
for(i in seq_along(run_ids)) {
  # Create sub-directories for each run.
  dir.create(file.path(settings$rundir, run_ids[i]), recursive=TRUE)
  dir.create(file.path(settings$modeloutdir, run_ids[i]), recursive=TRUE)
  
  # Extract parameter value for current run.
  param_curr <- lapply(param_design_list_by_pft, 
                       function(pft_params) as.data.frame(pft_params[i,,drop=FALSE]))
  
  # Write config to file.
  do.call(model_write_config, args=list(defaults = settings$pfts,
                                        trait.values = param_curr,
                                        settings = settings,
                                        run.id = run_ids[i]))
  
  # Write the run ID to a new line in the runs text file.
  cat(as.character(run_ids[i]),
      file = file.path(settings$rundir, "runs.txt"),
      sep = "\n",
      append = (i != 1L))
}

# run_ids <- PEcAn.assim.batch::pda.init.run(settings, con=NULL, model_write_config, 
#                                            workflow.id=-1, params=param_design_list_by_pft, 
#                                            n=N_design)


# ------------------------------------------------------------------------------
# Set initial conditions. 
#
# Soil IC:
# Fixing initial conditions at the values for a single ensemble member in the 
# IC ensemble, except for the initial soil pool, which is manually set to 
# align with an observed value or around 20. Ideally this would be done  
# in the call to `pda.init.run()` to avoid having to rewrite all of these files,
# but manually setting the `IC` parameter in the write config function seems 
# to effect the other non-soil initial conditions due to the SLA (specific 
# leaf area) calculation that is performed in the write config file. Making 
# the code more modular should make it easier to modify a single IC in this way.
#
# Wood/Roots ICs:
# To investigate the weird model behavior with respect to roots, I am also 
# manually varying the plantWoodInit initial condition. The value I have been 
# using is 156520.871478648 grams, so I create a grid of values around this 
# baseline value.
# ------------------------------------------------------------------------------

# Create grid of wood initial conditions.
wood_ics <- seq(10, 200000, length.out=N_design)
names(wood_ics) <- run_ids

# Randomly sample fineRootFrac and coarseRootFrac initial conditions. The sum 
# of these two fractions should be less than 1. 
fine_root_frac_ics <- runif(N_design, 0, 1)
coarse_root_frac_ics <- runif(N_design, rep(0, N_design), 1-fine_root_frac_ics)
names(fine_root_frac_ics) <- run_ids
names(coarse_root_frac_ics) <- run_ids

# Note that the IC file stores the soil IC in kg, while the value in 
# sipnet.param is converted to grams. The value `IC_soil` is in kg here, and 
# will be converted to grams when modifying the sipnet.param table.
IC_soil <- 1.6
IC_lai <- 3.4

for(run_id in run_ids) {
  # Read sipnet.param file.
  sipnet_param_path <- file.path(settings$rundir, run_id, "sipnet.param")
  sipnet_param <- fread(sipnet_param_path)
  
  # Update soil IC, converting from kilograms to grams.
  sipnet_param[V1=="soilInit", V2 := IC_soil*1000]
  
  # Update the wood initial condition.
  sipnet_param[V1=="plantWoodInit", V2 := wood_ics[[run_id]]]
  
  # Update LAI initial condition.
  sipnet_param[V1=="laiInit", V2 := IC_lai]
  
  # Update fine and coarse root fraction initial conditions.
  sipnet_param[V1=="fineRootFrac", V2 := fine_root_frac_ics[[run_id]]]
  sipnet_param[V1=="coarseRootFrac", V2 := coarse_root_frac_ics[[run_id]]]
  
  #
  # TEMP: trying to force zero fluxes out of roots.
  #
  
  sipnet_param[V1=="microbePulseEff", V2 := 1.0]
  sipnet_param[V1=="coarseRootTurnoverRate", V2 := 0.0]
  sipnet_param[V1=="fineRootTurnoverRate", V2 := 0.0]
  sipnet_param[V1=="baseCoarseRootResp", V2 := 0.0]
  sipnet_param[V1=="baseFineRootResp", V2 := 0.0]
  
  # Overwrite the param file with the modification.
  utils::write.table(as.data.frame(sipnet_param), sipnet_param_path, 
                     row.names=FALSE, col.names=FALSE, quote=FALSE)
}


# ------------------------------------------------------------------------------
# Execute forward model runs. 
# ------------------------------------------------------------------------------

PEcAn.workflow::start_model_runs(settings, write=FALSE)

