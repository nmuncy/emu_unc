#!/bin/bash

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -p IB_44C_512G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 16000
#SBATCH --job-name nm_fmriprep


module load singularity-3.8.2


# receive args
label=$1
sing_img=$2
proj_dir=$3

# set paths
dset_dir=${proj_dir}/dset
deriv_dir=${proj_dir}/derivatives
work_dir=${deriv_dir}/tmp_work/sub-${label}
mkdir -p $work_dir

# reference template fow, fs, avoid root issues
export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/data/madlab/singularity-images/templateflow
export FS_LICENSE=~/bin/licenses/fs_license.txt
cd /

# do job
singularity run --cleanenv $sing_img \
  $dset_dir \
  $deriv_dir \
  participant \
  --participant-label $label \
  --work-dir $work_dir \
  --skull-strip-template MNIPediatricAsym:cohort-5 \
  --output-spaces emur MNIPediatricAsym:cohort-5:res-2 \
  --nthreads 10 \
  --omp-nthreads 10 \
  --fs-license-file $FS_LICENSE \
  --stop-on-first-crash
