#!/bin/bash

# Construct and run a Joint Label Fusion command to create
# amygdala priors in template space with participant FreeSufer
# posteriors.
#
# Receives three positional argument from make_amgPriors_submit.sh:
#   [1] <out_dir> = location of final output
#   [2] <fs_dir> = location of FreeSurer posteriors
#   [3] <atlas> = template/atlas

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -p IB_44C_512G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 16000
#SBATCH --job-name amgJLF

# load relevant modules
module load ants-2.3.5
module load c3d-1.0.0-gcc-8.2.0
module load freesurfer-7.1

# receive args
out_dir=$1
fs_dir=$2
atlas=$3

# roi values
arrLabel=(18 54)

# start command - build as an array instead of using eval due to user input
jlfCmd=(antsJointLabelFusion.sh
    -d 3
    -t $atlas
    -o ${out_dir}/JLF_
    -p ${out_dir}/priors_amg/label_%04d.nii.gz
    -c 5
    -j 10
)

# orient in fs_dir
cd $fs_dir
for subj in sub-*; do

    # set output dir
    subj_out=${out_dir}/$subj
    mkdir $subj_out

    # set fs, out files
    fs_seg=${subj}/mri/aparc+aseg.mgz
    fs_anat=${subj}/mri/orig.mgz
    out_seg=${subj_out}/${subj}_desc-amg_mask.nii.gz
    out_anat=${subj_out}/${subj}_T1w.nii.gz

    # check for fs output
    if [ ! -f $fs_seg ] || [ ! -f $fs_anat ]; then
        continue
    fi
    echo -e "\nAdding $subj to JLF command ...\n"

    # get t1w anat
    if [ ! -f $out_anat ]; then
        mri_convert $fs_anat $out_anat
    fi

    # make out_seg
    if [ ! -f $out_seg ]; then

        # extract fs amg labels
        mri_convert $fs_seg ${subj_out}/aparc_aseg.nii.gz
        for h_label in ${arrLabel[@]}; do
            out_file=${subj_out}/label_${h_label}.nii.gz
            if [ ! -f $out_file ]; then
                c3d \
                    ${subj_out}/aparc_aseg.nii.gz \
                    -thresh $h_label $h_label 1 0 \
                    -o $out_file
            fi
        done

        # stitch together
        c3d \
            ${subj_out}/label_*.nii.gz \
            -accum -add -endaccum \
            -o $out_seg
    fi

    # update jlfCmd with files
    jlfCmd+=(-g $out_anat -l $out_seg)
done

# print for review, run JLF
echo -e "\nStarting:\n\t${jlfCmd[@]}\n"
"${jlfCmd[@]}"

# split L/R priors
atlas_str=tpl-MNIPediatricAsym_cohort-5_res-2
cp ${out_dir}/JLF_Labels.nii.gz ${out_dir}/${atlas_str}_desc-amg_mask.nii.gz
c3d \
    ${out_dir}/${atlas_str}_desc-amg_mask.nii.gz \
    -as SEG -cmv -pop -pop -thresh 50% inf 1 0 \
    -as MASK -push SEG -times \
    -o ${out_dir}/${atlas_str}_desc-amgR_mask.nii.gz \
    -push MASK -replace 1 0 0 1 -push SEG -times \
    -o ${out_dir}/${atlas_str}_desc-amgL_mask.nii.gz
