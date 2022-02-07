#!/bin/bash

#SBATCH --time=60:00:00   # walltime
#SBATCH --cpus-per-task=10
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
config_file=$1
pyAFQ $config_file --notrack
