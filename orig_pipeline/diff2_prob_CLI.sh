#!/bin/bash

#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=10   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=24gb   # memory per CPU core
#SBATCH -J "afqP"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab


# Notes:
#
# This script runs pyAFQ via the reference file
#   config.toml. Various files (templates, masks, etc.)
#   are required in $HOME/AFQ_data.

# make sure that AFQ_data exists in $HOME
if [ ! -d ${HOME}/AFQ_data ]; then
    cp -r /home/data/madlab/atlases/AFQ_data $HOME
fi

# Submit pyAFQ
pyAFQ config_prob.toml --notrack