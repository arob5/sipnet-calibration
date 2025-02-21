#
# construct_looser_priors.r
#
# This file is a temporary script that reads the output saved by 
# `get_sipnet_default_pars.r` and widens the prior bounds on parameters 
# related to wood and root dynamics. This is done to test whether the default
# bounds are much too restrictive, resulting in the observed model behavior 
# in which root carbon pools are plummetting. Not modifying the priors for the 
# allocation factors here, as these are jointly assigned a Dirichlet prior
# elsewhere.
#
# Andrew Roberts
#

library(data.table)

# General file paths. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration", 
                      "runs", "sipnet_calibration_2024")

# Paths to files containing the parameter design at which to run the forward model. 
prior_dir <- file.path(base_dir, "param", "priors")
plant_pft_path <- file.path(prior_dir, "sipnet_params_default_priors_plant.csv")
soil_pft_path <- file.path(prior_dir, "sipnet_params_default_priors_soil.csv")

# Update plant parameter priors.
param_plant <- fread(plant_pft_path)

param_plant[par_name=="veg_respiration_Q10", `:=`(parama=1.05, paramb=4.0)]
param_plant[par_name=="root_turnover_rate", `:=`(parama=0.00001, paramb=1.0)]
param_plant[par_name=="fine_root_respiration_Q10", `:=`(parama=1.05, paramb=7.0)]
param_plant[par_name=="coarse_root_respiration_Q10", `:=`(parama=1.05, paramb=7.0)]
param_plant[par_name=="coarse_root_respiration_Q10", `:=`(parama=1.05, paramb=7.0)]

fwrite(param_plant, file.path(prior_dir, "sipnet_params_default_priors_plant_loose.csv"))

# Update soil parameter priors. No modifications for now.
param_soil <- fread(soil_pft_path)
fwrite(param_soil, file.path(prior_dir, "sipnet_params_default_priors_soil_loose.csv"))






