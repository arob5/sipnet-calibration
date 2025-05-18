#!/bin/bash -l
#$ -N run_eki
#$ -l h_rt=72:00:00             # Hard time limit.
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.

# Runs the R script `run_eki.r` remotely via a qsub job submission.
# Passes ensemble size and number of iterations to the script.

N_EKI_ITR=3
N_ENS=200

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/sipnet_calibration/runs/calibrating_ics/scripts/run_eki.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --n_eki_itr=${N_EKI_ITR} \
    --n_ens=${N_ENS}
    