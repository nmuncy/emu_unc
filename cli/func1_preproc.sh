#!/bin/bash

# Notes ---
#
# sbatch command for running pre-processing via
# github.com/emu-project/func_preproc.git cli/afni_task_subj.py.
#
# Usage:
#   ./func1_preproc.sh

proj_dir="$(dirname "$(pwd)")"
code_dir=/home/nmuncy/compute/func_processing/cli
sbatch \
    --job-name=runAfniTask \
    --output=${code_dir}/logs/runAfniTask_log \
    --mem-per-cpu=4000 \
    --partition=IB_44C_512G \
    --account=iacc_madlab \
    --qos=pq_madlab \
    ${code_dir}/afni_task_subj.py \
    --batch-num 15 \
    --json-dir ${proj_dir}/data/timing_files \
    --out-dir /home/data/madlab/McMakin_EMUR01/derivatives/emu_unc \
    -c $code_dir \
    -s ses-S1 \
    -t task-study
