#!/bin/bash

module load r-3.5.1-gcc-8.2.0-djzshna

proj_dir=/scratch/madlab/emu_UNC
afni_dir=${proj_dir}/derivatives/afni
sess=ses-S2
task=test

subj_list=(`ls $afni_dir | grep "sub-*"`)
for subj in ${subj_list[@]}; do
    echo -e "\t Starting R script for $subj ..."
    write_dir=${afni_dir}/${subj}/${sess}/timing_files
    mkdir -p $write_dir
    Rscript func3_timing_files.R $proj_dir $subj $sess $task $write_dir
done
