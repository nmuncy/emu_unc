#!/bin/bash

sing_img=/home/nmuncy/bin/singularities/nipreps_fmriprep_20.2.3.simg

proj_dir=/scratch/madlab/emu_UNC
dset_dir=${proj_dir}/dset
deriv_dir=${proj_dir}/derivatives
subj_all=(`ls $dset_dir | grep "sub-*"`)

grow_list=true
unset subj_list
c=0; while [ $c -lt ${#subj_all[@]} ] && [ "$grow_list" = true ]; do

    check_file=${deriv_dir}/fmriprep/${subj_all[$c]}.html

    if [ ! -f $check_file ]; then
        subj_list+=(${subj_all[$c]})
    fi

    if [ ${#subj_list[@]} == 11 ]; then
        grow_list=false
    fi

    let c+=1
done


sbatch -e err.txt -o out.txt sing_job.sh
