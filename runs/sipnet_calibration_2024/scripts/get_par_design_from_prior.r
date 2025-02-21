#
# get_par_design_from_prior.r
#
# This function is the analog of `get_par_design_from_post.r`, but it generates
# the initial design from prior distributions, rather than by sub-sampling 
# previous MCMC output. 
# 
# The job of this script is to save CSV files storing a set of parameter
# values that will be read by `execute_sipnet_forward_runs.r`. One CSV 
# file is saved per PFT, which for SIPNET implies 2 files: one plant PFT and 
# one soil PFT. For each PFT, the CSV stores a matrix of dimension 
# (number design points, number PFT parameters), where the number of PFT 
# parameters is the number of parameters whose SIPNET default parameter 
# values will be (potentially) overwritten. This means the parameter set here
# includes both the calibration parameters and fixed parameters. The latter 
# will be fixed at their prior median if the prior is provided; otherwise 
# they will not be included in the saved files and the SIPNET code will fall 
# back on the default parameter value. The only parameters whose values should 
# vary across the rows of the CSV are the calibration parameters. The fixed 
# parameters will be constant across the rows. The format of the saved matrices 
# should follow the format of the matrices in the list returned by 
# `PEcAn.assim.batch::pda.generate.knots`. 
#
# NOTE: for now saving files only containing the calibration parameters; 
#       ignoring priors for fixed parameters for simplicity. 
# NOTE: My function `get_batch_design()` currently only supports sampling from 
#       from distributions that are the product of independent marginals. Below 
#       we consider a Dirichlet prior for the allocation parameters (which 
#       must sum to 1), so we manually sample these values and attach them to 
#       the Latin Hypercube design post-hoc. This implies that the prior 
#       distributions in the prior CSV files for the allocation parameters are 
#       ignored, and are instead specified here. These parameters should still
#       be included in the XML and the prior files though.
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
local_edits_path <- file.path(base_dir, "..", "..", "local_edits", 
                              "local_edits.r")
helper_funcs_path <- file.path(base_dir, "..", "..", "local_edits", 
                               "helper_functions_temp.r")
pecan_settings_path <- file.path(base_dir, "sipnet_calibration.xml")
design_func_path <- file.path(base_dir, "..", "..", "..", "gp-calibration", 
                              "src", "seq_design.r")
out_dir <- file.path(base_dir, "param", "param_design")

# Paths to the CSV files containing the information defining the prior 
# distributions. The prior distributions for the calibration parameters must 
# be present - prior distributions for the fixed parameters can optionally be 
# provided, in which case the parameter values will be fixed at the prior 
# median. 
plant_pft_path <- file.path(base_dir, "param", "priors",
                            "sipnet_params_default_priors_plant_loose.csv")
soil_pft_path <- file.path(base_dir, "param", "priors", 
                           "sipnet_params_default_priors_soil_loose.csv")
                           
source(local_edits_path)
source(helper_funcs_path)
source(design_func_path)

# ------------------------------------------------------------------------------
# Specify joint prior over allocation parameters.
# This is currently a bit fragmented, since the priors which are independent 
# are read from the prior CSV files, while these parameters are instead set 
# here. Moreover, "coarse_root_allocation_fraction" is not actually a SIPNET
# parameter, since it is determined immediately once the other 3 allocation 
# parameters are fixed. The samples from this variable are thus excluded 
# from the design. "root_allocation_fraction" maps to the SIPNET parameter
# "fineRootAllocation".
# ------------------------------------------------------------------------------

# Define Dirichlet distribution with concentration parameter alpha.
alloc_pars <- c("wood_allocation_fraction", "root_allocation_fraction", 
                "leaf_allocation_fraction", "coarse_root_allocation_fraction")
exclude_par <- "coarse_root_allocation_fraction"
alpha <- setNames(c(0.25, 0.25, 0.25, 0.25), alloc_pars)

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
# ------------------------------------------------------------------------------

# Load prior for plant PFT parameters.
plant_pft_prior <- fread(plant_pft_path)

# Changing column names to align with the requirement for `get_batch_design()`.
setnames(plant_pft_prior, c("distn","parama", "paramb"), 
                          c("dist", "param1", "param2"))
plant_pft_prior[dist=="unif", dist := "Uniform"]
plant_pft_prior <- as.data.frame(plant_pft_prior)
rownames(plant_pft_prior) <- plant_pft_prior$par_name

# Remove allocation parameters, as these will be sampled separately.
non_alloc_pars <- setdiff(plant_pft_prior$par_name, alloc_pars)
plant_pft_prior <- plant_pft_prior[non_alloc_pars,]

# Generate Latin Hypercube Design for plant PFT parameters, other than 
# allocation parameters.
plant_pft_design <- get_batch_design("LHS", N_batch=N_design, 
                                     prior_params=plant_pft_prior)

# Draw independent samples from Dirichlet prior over allocation parameters.
alloc_par_design <- LaplacesDemon::rdirichlet(n=N_design, alpha)
colnames(alloc_par_design) <- alloc_pars
alloc_par_design <- alloc_par_design[, setdiff(alloc_pars, exclude_par)]

# Construct complete design.
plant_pft_design <- cbind(plant_pft_design, alloc_par_design)

# Save CSV. 
fwrite(plant_pft_design, file.path(out_dir, "param_design_plant_pft.csv"))

# ------------------------------------------------------------------------------
# Generate parameter design: Soil PFT
#    In contrast to the plant PFT, the soil PFT parameter samples file contains 
#    direct MCMC output. Therefore, relevant considerations should be taken 
#    into account when sub-sampling (warm-up, autocorrelation, etc.). I am 
#    not being very format about this in this test; I checked the trace plots 
#    and the MCMC output looked reasonable. I am therefore simply sub-sampling 
#    from the second half of each chain. 
# ------------------------------------------------------------------------------

# Load prior for plant PFT parameters.
soil_pft_prior <- fread(soil_pft_path)

# Changing column names to align with the requirement for `get_batch_design()`.
setnames(soil_pft_prior, c("distn","parama", "paramb"), 
                         c("dist", "param1", "param2"))
soil_pft_prior[dist=="unif", dist := "Uniform"]
soil_pft_prior <- as.data.frame(soil_pft_prior)
rownames(soil_pft_prior) <- soil_pft_prior$par_name

# Generate Latin Hypercube Design for plant PFT parameters.
soil_pft_design <- get_batch_design("LHS", N_batch=N_design, 
                                    prior_params=soil_pft_prior)

# Save CSV. 
fwrite(soil_pft_design, file.path(out_dir, "param_design_soil_pft.csv"))

