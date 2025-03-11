#!/bin/bash -l
#$ -N run_exact_mcmc
#$ -l h_rt=72:00:00             # Hard time limit.
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.

# Runs the R script `run_exact_mcmc.r` remotely via a qsub job submission.
# Passes chain index and number of iterations to the script. This script is 
# intended to be called by the wrapper `run_mcmc_chains.sh`.

# Read commandline arguments.
N_MCMC_ITR=$1
CHAIN_IDX=$2

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/sipnet_calibration/runs/calibrating_ics/scripts/run_exact_mcmc.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --n_mcmc_itr=${N_MCMC_ITR} \
    --chain_idx=${CHAIN_IDX}
    
    