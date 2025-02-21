#
# pecan_first_sipnet_test.r
# Running a single forward evaluation of SIPNET at a single site. 
#
# Andrew Roberts
# 
# Working directory: sipnet_test_run
#

# TODO: 
# - Figure out how to access model outputs after running models (pda.get.model.output?, runModule.get.results?)
# - Figure out how to specify parameter values without sub-sampling. 
# - Figure out how to set up priors. 
# - Figure out how to incorporate calibration data (load.pda.data). 

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


#
# Setup. 
#

base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration", 
                      "runs", "sipnet_test_run")
local_edits_path <- file.path(base_dir, "..", "..", "local_edits", "local_edits.r")
pecan_settings_path <- file.path(base_dir, "pecan_first_sipnet_test.xml")
param_samples_path <- file.path("/projectnb", "dietzelab", "dongchen", 
                                "All_NEON_SDA", "NEON42", "SDA", "samples.Rdata")

source(local_edits_path)

#
# Read settings file. 
#

settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path)

#
# Prepare for model run. 
#

# Prepare/validate settings. 
settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)

# Create directories. 
dir.create(settings$outdir)
dir.create(settings$rundir)
dir.create(settings$modeloutdir)

# Write validated settings. 
PEcAn.settings::write.settings(settings, outputfile="pecan.CHECKED.xml")

# Write run files to prepare for model execution. 
# This function reads a set of parameter samples (either from the db or file), 
# then sub-samples them using a method specified by `ens.sample.method`. 
# Using my edited function in place of `PEcAn.workflow::RunModule.run.write.configs`. 
settings <- run.write.configs_test_edit(settings, overwrite=TRUE)
                            
#
# Run model. 
#

PEcAn.workflow::start_model_runs(settings, write=FALSE, stop.on.error=TRUE)

#
# Read model output. 
#

# Output variables to read. 
output_vars <- c("NEE", "LAI", "GPP", "TotSoilCarb", 
                 "leaf_carbon_content", "fine_root_carbon_content", "AGB")

# Used to specify path to read model outputs from. 
run_ids <- readLines(file.path(settings$rundir, "runs.txt"))

# List of matrices storing model output, one element per ensemble run (i.e., per run ID). 
output_list <- list()

for(run_id in run_ids) {
  run_id_path <- file.path(settings$modeloutdir, run_id)
  model_output <- PEcAn.utils::read.output(run_id, outdir=run_id_path, 
                                           variables=output_vars, dataframe=TRUE)
  model_output <- as.data.table(model_output)
  model_output[, run_id := run_id]
  output_list[[run_id]] <- model_output
}

# Combine into single data.table. 
model_output <- rbindlist(output_list, use.names=TRUE)


#
# Visualize model output.  
#

time_steps <- unique(model_output$posix)
time_steps <- time_steps[order(time_steps)]
N_time_steps <- length(time_steps)

ggplot(data=model_output) + 
  geom_line(aes(x=posix, y=NEE, color=run_id)) + 
  xlab("time")

ggplot(data=model_output[year==2021]) + 
  geom_line(aes(x=posix, y=NEE, color=run_id)) + 
  xlab("time")

ggplot(data=model_output[posix %in% tail(posix, 56)]) + 
  geom_line(aes(x=posix, y=NEE, color=run_id)) + 
  xlab("time")

ggplot(data=model_output[year==2021]) + 
  geom_line(aes(x=posix, y=leaf_carbon_content, color=run_id)) + 
  xlab("time")




