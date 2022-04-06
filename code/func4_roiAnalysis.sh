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

    Find subjects who have PPI output (ref func3_ppi.py). Then multiply ROI <mask_name>
    by intersection mask (-g) to make clean mask. Finally, use clean mask to extract
    coefs for each subject.

    Required Arguments:
        -d </path/to/dir> = location of project derivatives directory, should contain template dir.
        -g <group_intx_mask> = group intersection mask
        -m <mask_name> = identifying mask name
            (e.g. NSlacc to find <mask_dir>/tpl-MNIPediatricAsym_cohort-5_res-2_desc-NSlacc_mask.nii.gz)
        -n <decon_name> = identifying deconvolution name
        -p <ppi_seed> = identifying PPI seed name
            (e.g. amgL to find <data_dir>/<subj>/<sess>/func/decon_task-test_UniqueBehs_PPI-amgL_stats_REML+tlrc.HEAD)
        -s <session> = BIDS session string
        -t <task> = BIDS task string
        <behaviors> = remaining args are sub-brick behaviors to extract (ref 3dinfo -verb),
            "foo" for foo#0_Coef

    Example Usage:
        deriv_dir=/home/data/madlab/McMakin_EMUR01/derivatives/emu_unc
        sess=ses-S1
        task=task-study
        sbatch func4_roiAnalysis.sh \\
            -d \$deriv_dir \\
            -g \${deriv_dir}/emu_unc/template/tpl-MNIPediatricAsym_cohort-5_res-2_\${sess}_\${task}_desc-grpIntx_mask.nii.gz \\
            -m NSlacc \\
            -n precTest \\
            -p amgL \\
            -s \$sess \\
            -t \$task \\
            SPnegLF SPneuLF

USAGE
}

# receive args
while getopts ":d:g:m:n:p:s:t:h" OPT; do
    case $OPT in
    d)
        deriv_dir=${OPTARG}
        if [ ! -d ${deriv_dir} ] || [ ! -d ${deriv_dir}/template ]; then
            echo -e "\n\t ERROR: did not detect -d directory or required template sub-directory." >&2
            Usage
            exit 1
        fi
        ;;
    g)
        group_mask=${OPTARG}
        ;;
    m)
        mask_name=${OPTARG}
        ;;
    n)
        decon_name=${OPTARG}
        ;;
    p)
        ppi_seed=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    t)
        task=${OPTARG}
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
    deriv_dir)
        h_ret="-d"
        ;;
    mask_name)
        h_ret="-m"
        ;;
    decon_name)
        h_ret="-n"
        ;;
    sess)
        h_ret="-s"
        ;;
    task)
        h_ret="-t"
        ;;
    ppi_seed)
        h_ret="-p"
        ;;
    group_mask)
        h_ret="-g"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in deriv_dir sess task mask_name decon_name group_mask ppi_seed; do
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
        -d : $deriv_dir
        -g : $group_mask
        -m : $mask_name
        -n : $decon_name
        -p : $ppi_seed
        -s : $sess
        -t : $task
        beh : ${beh_list[@]}

EOF

# set up
mask_dir=${deriv_dir}/template
analysis_dir=${deriv_dir}/analyses
mkdir -p $analysis_dir

# find subjs with PPI output
subj_list_all=($(ls $deriv_dir | grep "sub-*"))
subj_list=()
for subj in ${subj_list_all[@]}; do
    check_file=${deriv_dir}/${subj}/${sess}/func/decon_${task}_${decon_name}_PPI-${ppi_seed}_stats_REML+tlrc.HEAD
    if [ -f $check_file ]; then
        subj_list+=($subj)
    fi
done
echo -e "Subject list:\n\t${subj_list[@]}"

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
out_file=${analysis_dir}/Coefs_${sess}_${task}_${ppi_seed}-${mask_name}.txt
echo -e "\n\tWriting: $out_file ...\n"
echo -e "Mask\t$mask_name" >$out_file

for subj in ${subj_list[@]}; do
    ppi_file=${deriv_dir}/${subj}/${sess}/func/decon_${task}_${decon_name}_PPI-${ppi_seed}_stats_REML+tlrc

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

    # convert brick array to comma-delimited list
    IFS=","
    brick_str="${brick_list[*]}"
    IFS=$' \t\n'

    # get, print coefs
    stats=$(3dROIstats -mask $mask_clean "${ppi_file}[${brick_str}]")
    echo -e "$subj\t$stats" >>$out_file
done
