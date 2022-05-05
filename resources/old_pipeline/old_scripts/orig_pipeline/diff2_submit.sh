#!/bin/bash

# submit jobs
time=`date '+%Y_%m_%d-%H:%M'`
out_dir=/scratch/madlab/emu_UNC/derivatives/Slurm_out/pyafq_${time}
mkdir -p $out_dir

sbatch \
    -e ${out_dir}/err.txt \
    -o ${out_dir}/out.txt \
    diff2_prob_CLI.sh
