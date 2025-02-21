#
# get_par_design_from_post.r
#
# The job of this script is to save CSV files storing a set of parameter
# values that will be read by `execute_sipnet_forward_runs.r`. One CSV 
# file is saved per PFT, which for SIPNET implies 2 files: one plant PFT and 
# one soil PFT. For each PFT, the CSV stores a matrix of dimension 
# (number design points, number PFT parameters), where the number of PFT 
# parameters is the number of parameters whose SIPNET default parameter 
# values will be (potentially) overwritten. This means the parameter set here
# includes both the calibration parameters and fixed parameters, where the 
# latter typically means the parameters are fixed at their prior median. 
# The only parameters whose values should vary across the rows of the CSV 
# are the calibration parameters. The fixed parameters will be constant 
# across the rows. The format of the saved matrices should follow the format 
# of the matrices in the list returned by `PEcAn.assim.batch::pda.generate.knots`. 
#
# For this test forward model run, the values for the calibration parameters 
# are sampled from previous parameter posteriors - hence the "from_post" in the 
# file name for this script. The specifics are discussed in below comments. 
#
# Andrew Roberts
#
# Dependencies: 
#  dietzelab/arober/sipnet_calibration/local_edits/local_edits.r
#  dietzelab/arober/sipnet_calibration/local_edits/helper_functions_temp.r
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
out_dir <- file.path(base_dir, "param", "param_design")

# Paths to previous parameter posteriors that will be sub-sampled to generate 
# the design points. Note the below comments on the difference between the 
# plant and soil PFT files. The plant PFT path was copied from:  
# /projectnb/dietzelab/dongchen/All_NEON_SDA/NEON42/SDA/samples.Rdata
# while the soil PFT path was copied from: 
# /projectnb/dietzelab/dietze/hf_landscape_SDA/test03/pfts/trait.mcmc.pda.soil.HPDA_som_respiration_rate_correction.Rdata
plant_pft_path <- file.path(base_dir, "param", "prev_posteriors", "param_samples_plant_pfts.Rdata")
soil_pft_path <- file.path(base_dir, "param", "prev_posteriors", "param_samples_soil.HPDA.Rdata")

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
# Site, PFT, and design point settings. 
# ------------------------------------------------------------------------------

# Look up plant PFT for the specified site. 
site_pft_lookup <- fread(settings$run$inputs$pft.site$path)
site_id <- settings$run$site$id
plant_pft_name <- site_pft_lookup[site==site_id, pft]
if(length(plant_pft_name)==0) stop("No PFT found for site ID ", site_id)

# Soil PFT specified in settings. 
soil_pft_name <- grep("soil", names(settings$assim.batch$param.names), value=TRUE)
if(length(soil_pft_name) != 1L) {
  stop("Must have exactly one soil PFT; found: ", paste(soil_pft_name, sep=", "))
}

# Number design points. 
N_design <- as.integer(settings$ensemble$size)

# Calibration parameter names, for each PFT. 
param_cal_names <- lapply(settings$assim.batch$param.names, unlist)

# ------------------------------------------------------------------------------
# Generate parameter design: Plant PFT
#   The plant PFT path specified points to parameter samples that have already 
#   been sub-sampled from previous MCMC output (down to 100 samples). Thus, 
#   here we sub-sample further if `N_design` < 100. 
# ------------------------------------------------------------------------------

# Load samples, which include various PFTs, only one of which should be 
# selected below. 
plant_pft_env <- new.env()
load(plant_pft_path, env=plant_pft_env)
plant_pft_samp <- plant_pft_env[[names(plant_pft_env)]]

# Select plant PFT. 
plant_pft_samp <- plant_pft_samp[[plant_pft_name]]

# Select only columns corresponding to calibration parameters. 
plant_pft_samp <- plant_pft_samp[, param_cal_names[[plant_pft_name]]]

# Sub-sample further to obtain the number of specified design points. 
N_prev_samp <- nrow(plant_pft_samp)
if(N_design > N_prev_samp) {
  stop("Number of design points ", N_design, 
       " exceeds number of available samples ", N_prev_samp)
} else if (N_design < N_prev_samp) {
  plant_pft_samp <- plant_pft_samp[sample(seq_len(N_prev_samp), size=N_design, replace=FALSE),]
}

# Save CSV. 
fwrite(plant_pft_samp, file.path(out_dir, "param_design_plant_pft.csv"))

# ------------------------------------------------------------------------------
# Generate parameter design: Soil PFT
#    In contrast to the plant PFT, the soil PFT parameter samples file contains 
#    direct MCMC output. Therefore, relevant considerations should be taken 
#    into account when sub-sampling (warm-up, autocorrelation, etc.). I am 
#    not being very format about this in this test; I checked the trace plots 
#    and the MCMC output looked reasonable. I am therefore simply sub-sampling 
#    from the second half of each chain. 
# ------------------------------------------------------------------------------

# Load soil parameter samples. 
soil_pft_env <- new.env()
load(soil_pft_path, env=soil_pft_env)
soil_pft_samp <- soil_pft_env[[names(soil_pft_env)]]

# Select soil calibration parameters. 
soil_pft_samp <- soil_pft_samp[param_cal_names[[soil_pft_name]]]
assert_that(setequal(names(soil_pft_samp), param_cal_names[[soil_pft_name]]), 
            msg="Soil PFT calibration parameter names in settings not aligning with names in samples file.")

# Convert to list of matrices (from coda MCMC list), and then to single matrix. 
soil_cal_param_names <- names(soil_pft_samp)
soil_pft_samp <- lapply(soil_pft_samp, function(x) as.matrix(coda::as.mcmc(x)))
soil_pft_samp <- do.call(cbind, soil_pft_samp)
colnames(soil_pft_samp) <- soil_cal_param_names

# Sub-sample. 
subsample_idx <- sample(1:nrow(soil_pft_samp), size=N_design, replace=FALSE)
soil_pft_samp <- as.data.table(soil_pft_samp[subsample_idx,])

# Save CSV. 
fwrite(soil_pft_samp, file.path(out_dir, "param_design_soil_pft.csv"))

