#!/bin/bash

#SBATCH --time=00:30:00   # walltime
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "emuM"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab

module load c3d-1.0.0-gcc-8.2.0

scratch_dir=/scratch/madlab/emu_unc/derivatives/kmeans_warp
deriv_dir=/home/data/madlab/McMakin_EMUR01/derivatives/afni
sess=ses-S2

# find subjects with mask
subj_all=($(ls $deriv_dir | grep "sub-*"))
file_list=()
for subj in ${subj_all[@]}; do
    check_file=${deriv_dir}/${subj}/${sess}/dwi/${subj}_${sess}_space-MNIPediatricAsym_cohort-5_res-2_desc-bla_mask.nii.gz
    if [ -f $check_file ]; then
        file_list+=($check_file)
    fi
done

# combine masks
c3d \
    ${file_list[@]} \
    -accum -add -endaccum \
    -o ${scratch_dir}/bla_mask_sum.nii.gz
