#!/bin/bash

#SBATCH --time=00:30:00   # walltime
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "emuM"   # job name
#SBATCH -p IB_44C_512G   # partition name
#SBATCH --account iacc_madlab  # account
#SBATCH --qos pq_madlab

function Usage {
    cat <<USAGE

    Use BLA masks produced by func0_masks.py to make left and right bla masks.
    This involves reslicing each mask into reference space, summing the masks,
    thresholding the summation at some level (20), and then splitting the resulting
    binarized mask into left/right by the origin.

    Final masks are written to -w <scratch working dir>.

    Visual inspection revealed registration problems for a number of subjects
    (possibly due to issues with distortion correction), these are excluded
    from the group-level masks.

    Required Arguments:
        -d <project_derivatives> = path to project derivatives location
        -w <scratch_directory> = path to working scratch/derivatives/kmeans_warp directory
        -s <session> = BIDS session string

    Example Usage:
        sbatch \\
            -e err.log \\
            -o out.log \\
            func1_combine.sh \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives \\
            -w /scratch/madlab/emu_unc/derivatives/kmeans_warp \\
            -s ses-S2
USAGE
}

# capture arguments
while getopts ":d:s:w:h" OPT; do
    case $OPT in
    d)
        deriv_dir=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    w)
        scratch_dir=${OPTARG}
        ;;
    h)
        Usage
        exit 0
        ;;
    \?)
        echo -e "\n\t ERROR: invalid option." >&2
        Usage
        exit 1
        ;;
    esac
done

# print help if no arg
if [ $OPTIND == 1 ]; then
    Usage
    exit 0
fi

# check args
if [ ! -d $deriv_dir ]; then
    echo -e "\n\t ERROR: -d directory not found.\n" >&2
    Usage
    exit 1
fi

if [ -z $scratch_dir ]; then
    echo -e "\n\t ERROR: please specify -w directory.\n" >&2
    Usage
    exit 1
fi

if [ -z $sess ]; then
    echo -e "\n\t ERROR: please specify -s session.\n" >&2
    Usage
    exit 1
fi

# load required modules (c3d)
module load c3d-1.0.0-gcc-8.2.0

# set up paths
afni_dir=${deriv_dir}/afni
out_dir=${deriv_dir}/emu_unc

# find subjects with func0_masks.py mask
subj_all=($(ls $afni_dir | grep "sub-*"))
file_list=()
subj_list=()
for subj in ${subj_all[@]}; do
    check_file=${afni_dir}/${subj}/${sess}/dwi/${subj}_${sess}_space-MNIPediatricAsym_cohort-5_res-2_desc-bla_mask.nii.gz
    if [ -f $check_file ]; then
        echo -e "\t Found file for $subj"
        file_list+=($check_file)
        subj_list+=($subj)
    fi
done

# find index of subjects w/registration problems
problem_list=(sub-4{008,012,019,038,057,065,082,083,091,108,110,114,118,119,125,142,144})
ind_problem=()
for ind in ${!subj_list[@]}; do
    for prob in ${problem_list[@]}; do
        if [[ "${subj_list[$ind]}" == "${prob}" ]]; then
            ind_problem+=($ind)
        fi
    done
done

# remove problematic subjects from arrays
for ind in "${ind_problem[@]}"; do
    unset "subj_list[$ind]"
    unset "file_list[$ind]"
done
declare -a subj_list=(${subj_list[@]})
declare -a file_list=(${file_list[@]})

# forcibly solve space/origin issue - use first subject as reference
ref_file=${file_list[0]}
mask_list=($ref_file)
c=1
while [ $c -lt ${#file_list[@]} ]; do
    subj=${subj_list[$c]}
    mask_file=${file_list[$c]}
    reslice_file=${scratch_dir}/${subj}/${sess}/dwi/tmp_bla_resliced.nii.gz
    if [ ! -f $reslice_file ]; then
        echo -e "\t Reslicing $subj"
        c3d \
            $ref_file \
            $mask_file \
            -reslice-identity \
            -o $reslice_file
    fi
    mask_list+=($reslice_file)
    let c+=1
done

# combine masks
if [ ! -f ${scratch_dir}/bla_mask_sum.nii.gz ]; then
    echo -e "\n\t Stitching together ${scratch_dir}/bla_mask_sum.nii.gz"
    c3d \
        ${mask_list[@]} \
        -accum -add -endaccum \
        -o ${scratch_dir}/bla_mask_sum.nii.gz
fi

# thresh, split LR
echo -e "\n\t Making binary ${scratch_dir}/bla_mask_<left/right>.nii.gz"
c3d \
    ${scratch_dir}/bla_mask_sum.nii.gz \
    -thresh 20 inf 1 0 \
    -o ${scratch_dir}/bla_mask_bin.nii.gz

c3d \
    ${scratch_dir}/bla_mask_bin.nii.gz \
    -as SEG -cmv -pop -pop -thresh 50% inf 1 0 -as MASK \
    -push SEG -times \
    -o ${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-blaL_mask.nii.gz \
    -push MASK -replace 1 0 0 1 \
    -push SEG -times \
    -o ${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-blaR_mask.nii.gz
