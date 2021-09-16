#!/bin/bash

# TODO ranem anat files, copy new event files

# orient to existing data
source_dir=/home/data/madlab/McMakin_EMUR01
source_dset=${source_dir}/dset
source_dwi=${source_dir}/derivatives/dwi_preproc

# start project dir
proj_dir=/scratch/madlab/emu_UNC
proj_dset=${proj_dir}/dset
proj_deriv=${proj_dir}/derivatives
mkdir -p $proj_deriv $proj_dset

# extras vars for ease
sess=ses-S2

# get bids helper files
cp ${source_dset}/dataset_description.json $proj_dset
cp ${source_dset}/participants.tsv $proj_dset
cp ${source_dset}/task-test_bold.json $proj_dset

# only pull data when pre-processed DWI exists
cd $source_dwi
for subj in sub-*; do

    # files to be moved
    dwi_file=${subj}/${sess}/dwi/${subj}_${sess}_run-1_desc-eddyCorrected_dwi.nii.gz
    dwi_bvec=${subj}/${sess}/dwi/${subj}_${sess}_run-1_desc-eddyCorrected_dwi.bvec
    dwi_bval=${source_dset}/${subj}/${sess}/dwi/${subj}_${sess}_run-1_dwi.bval
    anat_file=${source_dset}/${subj}/ses-S1/anat/${subj}_ses-S1_run-2_T1w.nii.gz
    func_list=(`ls ${source_dset}/${subj}/${sess}/func/${subj}_${sess}_task-test_run-*`)
    fmap_list=(`ls ${source_dset}/${subj}/${sess}/fmap/${subj}_${sess}_acq-func_dir-*_run-?_epi.*`)

    # don't repeat work
    dwi_check=${proj_deriv}/dwi_preproc/${subj}/${sess}/dwi/${subj}_${sess}_run-1_desc-eddyCorrected_dwi.nii.gz
    if [ -f $dwi_check ]; then
        continue
    fi

    # make sure relevant files exist
    if [ ! -f $dwi_file ] || [ ! -f $anat_file ] || [ ${#func_list[@]} != 9 ]; then
        continue
    fi

    # print progress
    echo -e "\t Copying data for $subj ..."

    # set file destinations in new proj dir
    proj_dwi=${proj_deriv}/dwi_preproc/${subj}/${sess}/dwi
    proj_anat=${proj_dset}/${subj}/${sess}/anat
    proj_func=${proj_dset}/${subj}/${sess}/func
    proj_fmap=${proj_dset}/${subj}/${sess}/fmap
    mkdir -p $proj_dwi $proj_anat $proj_func $proj_fmap

    # copy data to proj dir, rename some
    cp $dwi_file ${proj_dwi}/${subj}_${sess}_dwi.nii.gz
    cp $dwi_bvec ${proj_dwi}/${subj}_${sess}_dwi.bvec
    cp $dwi_bval ${proj_dwi}/${subj}_${sess}_dwi.bval

    cp ${anat_file%%.*}* $proj_anat

    for tmp_file in ${func_list[@]}; do
        cp $tmp_file $proj_func
    done

    for tmp_file in ${fmap_list[@]}; do
        cp $tmp_file $proj_fmap
    done
done
