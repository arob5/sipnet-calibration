#!/bin/bash -l
#$ -N run_prior_ensemble
#$ -l h_rt=12:00:00             # Hard time limit.
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.

# Runs the R script `run_prior_ensemble.r`.
 
# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/sipnet_calibration/runs/calibrating_ics/scripts/run_prior_ensemble.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH}