#!/bin/bash

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -p IB_44C_512G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16000
#SBATCH --job-name ppiROI

# load relevant modules
module load afni-20.2.06
module load c3d-1.0.0-gcc-8.2.0

# help
function Usage {
    cat <<USAGE
    Extract PPI coefficient for supplied behaviors.

    First, find subjects who have PPI output (ref func2_ppi.py).
    Second, make a group intersection mask referencing *desc-intersect_mask.nii.gz
        of each subject (ref resources.afni.masks.make_intersect_mask from
        github.com/emu-project/func_processing.git).
    Third, multiply ROI <mask_name> by intersection mask to make clean mask.
    Fourth, use clean mask to extract coefs for each subject.

    Assumes example strings/locations described in required args.

    Required Arguments:
        -d </path/to/dir> = location of project derivatives directory, should contain both
            afni and emu_unc sub-directories.
        -m <mask_name> = identifying mask name
            (e.g. NSdmpfcL to find <mask_dir>/tpl-MNIPediatricAsym_cohort-5_res-2_desc-NSdmpfcL_mask.nii.gz)
        -p <ppi_seed> = identifying PPI seed name
            (e.g. amgL to find <data_dir>/<subj>/<sess>/func/decon_task-test_UniqueBehs_PPI-amgL_stats_REML+tlrc.HEAD)
        -s <session> = BIDS session string
        <behaviors> = remaining args are sub-brick behaviors to extract (ref 3dinfo -verb)

    Example Usage:
        sbatch func3_roiAnalysis.sh \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives \\
            -m NSdmpfcL \\
            -p amgL \\
            -s ses-S2 \\
            SnegLF SneuLF

USAGE
}

# receive args
while getopts ":d:m:p:s:h" OPT; do
    case $OPT in
    d)
        proj_dir=${OPTARG}
        if [ ! -d ${proj_dir}/afni ] || [ ! -d ${proj_dir}/emu_unc ]; then
            echo -e "\n\t ERROR: did not detect $proj_dir or required sub-directories." >&2
            Usage
            exit 1
        fi
        ;;
    m)
        mask_name=${OPTARG}
        ;;
    p)
        ppi_seed=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    h)
        Usage
        exit 0
        ;;
    :)
        echo -e "\n\t ERROR: option '$OPTARG' missing argument." >&2
        Usage
        exit 1
        ;;
    \?)
        echo -e "\n\t ERROR: invalid option '$OPTARG'." >&2
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

# get remaning args as behaviors
shift "$((OPTIND - 1))"
beh_list=("$@")

# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
    proj_dir)
        h_ret="-d"
        ;;
    mask_name)
        h_ret="-m"
        ;;
    sess)
        h_ret="-s"
        ;;
    ppi_seed)
        h_ret="-p"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in proj_dir sess mask_name ppi_seed; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for passed behaviors
if [ ${#beh_list[@]} -eq 0 ]; then
    echo -e "\n\n\t ERROR: Missing behavior inputs." >&2
    Usage
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, options captured:
        -d : $proj_dir
        -m : $mask_name
        -p : $ppi_seed
        -s : $sess
        beh : ${beh_list[@]}

EOF

# # for testing
# sess=ses-S2
# ppi_seed=amgL
# proj_dir=/home/data/madlab/McMakin_EMUR01/derivatives
# mask_name=NSdmpfcL
# beh_list=(SnegLF SneuLF)

# set up
afni_dir=${proj_dir}/afni
ppi_dir=${proj_dir}/emu_unc
mask_dir=${ppi_dir}/template
analysis_dir=${ppi_dir}/analyses

# find subjs with PPI output
subj_list_all=($(ls $ppi_dir | grep "sub-*"))
subj_list=()
for subj in ${subj_list_all[@]}; do
    check_file=${ppi_dir}/${subj}/${sess}/func/decon_task-test_UniqueBehs_PPI-${ppi_seed}_stats_REML+tlrc.HEAD
    if [ -f $check_file ]; then
        subj_list+=($subj)
    fi
done
echo -e "Subject list:\n\t${subj_list[@]}"

# make group intersection mask
group_mask=${mask_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-grpIntx_mask.nii.gz
if [ ! -f $group_mask ]; then
    echo -e "\n\tMaking: $group_mask"
    maskCmd=(3dmask_tool
        -frac 1
        -prefix $group_mask
        -input
    )
    for subj in ${subj_list[@]}; do
        mask_file=${afni_dir}/${subj}/${sess}/anat/${subj}_${sess}_space-MNIPediatricAsym_cohort-5_res-2_desc-intersect_mask.nii.gz
        if [ -f $mask_file ]; then
            maskCmd+=($mask_file)
        fi
    done
    echo -e "Starting:\n\t${maskCmd[@]}"
    "${maskCmd[@]}"
fi

# multiply ROI mask by group intx mask, deal with origin issue, binarize
mask_clean=${mask_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-${mask_name}Clean_mask.nii.gz
if [ ! -f $mask_clean ]; then
    echo -e "\n\tMaking: $mask_clean"
    mask_raw=${mask_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-${mask_name}_mask.nii.gz
    c3d \
        $group_mask $mask_raw \
        -reslice-identity \
        -o $mask_raw
    c3d \
        $mask_raw $group_mask \
        -multiply \
        -o $mask_clean
    c3d \
        $mask_clean \
        -thresh 0.5 1 1 0 \
        -o $mask_clean
fi

# extract coefs for each subj/behavior
out_file=${analysis_dir}/Coefs_${ppi_seed}-${mask_name}.txt
echo -e "\n\tWriting: $out_file ...\n"
echo -e "Mask\t$mask_name" >$out_file

for subj in ${subj_list[@]}; do
    ppi_file=${ppi_dir}/${subj}/${sess}/func/decon_task-test_UniqueBehs_PPI-${ppi_seed}_stats_REML+tlrc

    # find sub-brick for behavior coefficient
    brick_list=()
    for beh in ${beh_list[@]}; do
        h_brick=$(3dinfo -label2index "${beh}#0_Coef" $ppi_file)
        brick_list+=($h_brick)
    done

    # make sure all planned behaviors were found
    if [ ${#brick_list[@]} != ${#beh_list[@]} ]; then
        echo -e "\n\t Missing behavior for $subj, skipping ..."
        continue
    fi

    # convert array to comma-delimited list
    IFS=","
    brick_str="${brick_list[*]}"
    IFS=$' \t\n'

    # get, print coefs
    stats=$(3dROIstats -mask $mask_clean "${ppi_file}[${brick_str}]")
    echo -e "$subj\t$stats" >>$out_file
done
