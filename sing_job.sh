#!/bin/bash

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH --nodes 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 16000
#SBATCH --job-name nate_sing

sing_dir=/home/nmuncy/bin/singularities
proj_dir=/scratch/madlab/test_sing
dset_dir=${proj_dir}/dset
deriv_dir=${proj_dir}/derivatives

label=4009
work_dir=${deriv_dir}/sub-${label}/tmp_work
mkdir -p $work_dir

export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/data/madlab/singularity-images/templateflow
cd /
singularity run --cleanenv ${sing_dir}/nipreps_fmriprep_20.2.3.simg \
  $dset_dir $deriv_dir \
  participant \
  --participant-label $label \
  --work-dir $work_dir \
  --skull-strip-template MNIPediatricAsym:cohort-5 \
  --output-spaces MNIPediatricAsym:cohort-5:res-1.8 \
  --nthreads 16 \
  --omp-nthreads 16 \
  --fs-license-file $FS_LICENSE