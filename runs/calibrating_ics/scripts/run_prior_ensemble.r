#
# run_prior_ensemble.r
# Forward evaluations of SIPNET using different parameter and initial condition
# (IC) values. A single model driver trajectory is used. The primary purpose is
# is to select a set of calibration parameters (treating ICs as parameters
# here) and prior distributions that put probability mass in regions of 
# parameter space that lead to a reasonable spread in model outputs.
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

# ------------------------------------------------------------------------------
# User-specified settings/filepaths.  
# ------------------------------------------------------------------------------

# Directories. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration")
pecan_base_dir <- file.path(base_dir, "runs", "calibrating_ics")
src_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration", "src")
pecan_src_dir <- file.path(base_dir, "src")

# PEcAn XML path
pecan_settings_path <- file.path(pecan_base_dir, "sipnet_calibration.xml")

# Sourcing helper functions.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "basis_function_emulation.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(pecan_src_dir, "param_calibration_functions.r"))
source(file.path(pecan_src_dir, "mcmc_pecan.r"))


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
# Specify parameters and initial conditions to vary.
#  Note that the names of the params/ICs must follow the PEcAn naming standards;
#  they will be converted to SIPNET parameters in `write.config.SIPNET()`.
#  Parameter names can be passed to this function by PFT; we make no PFT 
#  distinctions here and lump all parameters together.
# ------------------------------------------------------------------------------

# Parameter names. Note that `root_allocation_fraction` maps to the SIPNET 
# fine root allocation parameter. The parameter `coarse_root_allocation_fraction`
# does not actually exist; it is implicitly defined as 1 minus the sum of the 
# other allocation parameters. We introduce it here for the purpose of 
# defining a Dirichlet prior on the allocation parameters. It will be
# removed before being passed to the write config function.
alloc_pars <- c("root_allocation_fraction", "wood_allocation_fraction",
                "leaf_allocation_fraction", "coarse_root_allocation_fraction")
param_names <- c(alloc_pars, "root_turnover_rate", "Amax")
exclude_par <- "coarse_root_allocation_fraction"

# Initial conditions.
frac_pars <- c("abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac")
ic_names <- c(frac_pars, "AbvGrndWood", "soil")

all_par_names <- setdiff(c(param_names, ic_names), exclude_par)


# ------------------------------------------------------------------------------
# Specify prior distributions.
# ------------------------------------------------------------------------------

rprior_list <- list()

# Allocation parameters
rprior_list$alloc_pars <- function(n) {
  n_alloc_pars <- length(alloc_pars)
  
  alphas <- setNames(rep(1/n_alloc_pars, n_alloc_pars), alloc_pars)
  alloc_par_design <- LaplacesDemon::rdirichlet(n, alphas)
  colnames(alloc_par_design) <- alloc_pars
  alloc_par_design[, setdiff(alloc_pars, exclude_par), drop=FALSE]
}

# Root turnover rate
rprior_list$root_turnover_rate <- function(n) {
  matrix(rbeta(n, shape1=2, shape2=4), ncol=1, dimnames=list(NULL, "root_turnover_rate"))
}

# Amax
rprior_list$Amax <- function(n) {
  matrix(rgamma(n, shape=224, rate=2), ncol=1, dimnames=list(NULL, "Amax"))
}

# Fraction initial condition parameters.
rprior_list$frac_pars <- function(n) {
  n_frac_pars <- length(frac_pars)
  alphas <- setNames(rep(1/n_frac_pars, n_frac_pars), frac_pars)
  frac_par_design <- LaplacesDemon::rdirichlet(n, alphas)
  colnames(frac_par_design) <- frac_pars
  return(frac_par_design)
}

# Aboveground wood initial condition (in SIPNET not PEcAn units).
rprior_list$AbvGrndWood <- function(n) {
  matrix(rgamma(n, shape=10, rate=6.25e-05), ncol=1, dimnames=list(NULL, "AbvGrndWood"))
}

# Soil pool initial condition (in SIPNET not PEcAn units).
rprior_list$soil <- function(n) {
  matrix(rgamma(n, shape=10, rate=5e-04), ncol=1, dimnames=list(NULL, "soil"))
}



# ------------------------------------------------------------------------------
# Sample parameter values from prior. 
# ------------------------------------------------------------------------------

# Function to generate parameter design.
sample_prior <- function(n) {
  samp_list <- lapply(rprior_list, function(f) f(n))
  do.call(cbind, samp_list)
}
  
# Number of design points to sample from prior (i.e., the number of forward
# model runs).
n_design <- as.integer(settings$ensemble$size)

# Assemble design and ensure that all parameters are present.
par_design <- sample_prior(n_design)

if(!setequal(colnames(par_design), all_par_names)) {
  stop("Set of parameter names in `par_design` is not equal to `all_par_names`.")
}

# Split into parameters and ICs.
ic_design <- par_design[,ic_names, drop=FALSE]
par_design <- par_design[,setdiff(param_names, exclude_par), drop=FALSE]

# ------------------------------------------------------------------------------
# Prepare for forward runs: write one config per design point. 
# ------------------------------------------------------------------------------

# Load SIPNET package. 
model_write_config <- paste("write.config.", settings$model$type, sep = "")
PEcAn.utils::load.modelpkg(settings$model$type)

# Save write configs.
run_ids <- paste("run", 1:n_design, sep=".")
for(i in seq_along(run_ids)) {
  # Create sub-directories for each run.
  dir.create(file.path(settings$rundir, run_ids[i]), recursive=TRUE)
  dir.create(file.path(settings$modeloutdir, run_ids[i]), recursive=TRUE)
  
  # Extract parameter value for current run. Write config currently expects
  # a list of data.frames as an argument.
  par_curr <- list(par_design[i,,drop=FALSE])
  
  # Extract initial condition for current run. Must be named list.
  ic_curr <- as.list(ic_design[i,])
  
  # Write config to file.
  do.call(model_write_config, args=list(defaults = settings$pfts,
                                        trait.values = par_curr,
                                        IC = ic_curr,
                                        settings = settings,
                                        run.id = run_ids[i]))
  
  # Write the run ID to a new line in the runs text file.
  cat(as.character(run_ids[i]),
      file = file.path(settings$rundir, "runs.txt"),
      sep = "\n",
      append = (i != 1L))
}


# ------------------------------------------------------------------------------
# Execute forward model runs. 
# ------------------------------------------------------------------------------

# Run locally.
settings$host <- list(name="localhost")

PEcAn.workflow::start_model_runs(settings, write=FALSE)


# ------------------------------------------------------------------------------
# Read and save table of model outputs. 
# ------------------------------------------------------------------------------

# Output variables to read. 
output_vars <- c("time", "NEE", "LAI", "GPP", "TotSoilCarb", "AbvGrndWood", 
                 "leaf_carbon_content", "fine_root_carbon_content", "AGB", 
                 "coarse_root_carbon_content", "wood_carbon_content", 
                 "litter_carbon_content")

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

fwrite(model_output, file.path(settings$modeloutdir, "model_output.csv"))
















