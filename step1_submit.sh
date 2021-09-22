#!/bin/bash

# Set up
sing_img=/home/nmuncy/bin/singularities/nipreps_fmriprep_20.2.3.simg

proj_dir=/scratch/madlab/emu_UNC
dset_dir=${proj_dir}/dset
deriv_dir=${proj_dir}/derivatives
subj_all=(`ls $dset_dir | grep "sub-*"`)

# Make a list of 8 subjects who don't have fMRIprep
#   output in deriv_dir
max_num=8
sess=ses-S2
task=task-test
space=space-MNIPediatricAsym_cohort-5_res-2

unset subj_list
for subj in ${subj_all[@]}; do
    check_file=${deriv_dir}/fmriprep/${subj}/${sess}/func/${subj}_${sess}_${task}_run-1_${space}_desc-preproc_bold.nii.gz
    if [ ! -f $check_file ]; then
        subj_list+=(${subj})
    fi
    if [ ${#subj_list[@]} == $max_num ]; then
        break
    fi
done

# submit jobs
time=`date '+%Y_%m_%d-%H_%M'`
out_dir=${deriv_dir}/Slurm_out/fmriprep_${time}
mkdir -p $out_dir

for subj in ${subj_list[@]}; do

    # clean
    for check in freesurfer tmp_work fmriprep; do
        if [ -d ${deriv_dir}/${check}/${subj} ]; then
            rm -r ${deriv_dir}/${check}/${subj}*
        fi
    done

    # submit jobs
    sbatch \
        -e ${out_dir}/err_${subj}.txt \
        -o ${out_dir}/out_${subj}.txt \
        step1_fmriprep.sh \
            ${subj#*-} \
            $sing_img \
            $proj_dir
    sleep 1
done
