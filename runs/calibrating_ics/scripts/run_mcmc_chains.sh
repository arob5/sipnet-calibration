#!/bin/bash -l
#$ -N run_mcmc_chains
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# This script is a wrapper around `remote_exec_run_exact_mcmc.sh`, 
# intended to call this script multiple times in order to run different 
# MCMC chains remotely in parallel.
#
# Settings passed to `remote_exec_run_exact_mcmc.sh`:
# N_MCMC_ITR: Number of MCMC iterations.
# CHAIN_IDX: Integer index of MCMC chain.

# Settings for this script.
N_MCMC_ITR=50000
N_CHAINS=4

# Path to bash script to execute.
SCRIPT_PATH="/projectnb/dietzelab/arober/sipnet_calibration/runs/calibrating_ics/scripts/remote_exec_run_exact_mcmc.sh"

# Print settings.
echo "Commandline arguments:"
echo "N_MCMC_ITR: ${N_MCMC_ITR}"
echo "N_CHAINS: ${N_CHAINS}"

# Run run_init_emulator.sh` `N_CHAINS` times.
for i in $( seq 1 ${N_CHAINS} )
do
  qsub ${SCRIPT_PATH} ${N_MCMC_ITR} ${i}
done


