#!/bin/bash

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -p IB_44C_512G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16000
#SBATCH --job-name ppiINTX

# load relevant modules
module load afni-20.2.06
module load c3d-1.0.0-gcc-8.2.0

# help
function Usage {
    cat <<USAGE
    Group intersection mask.

    Make a group intersection mask referencing <sess>_<task>_*desc-intersect_mask.nii.gz
    of each subject (ref resources.afni.masks.make_intersect_mask from
    github.com/emu-project/func_processing.git).

    Optionally, multiply seed masks by intersection mask to generate foo_*Clean_mask.nii.gz.
    These masks will be resampled into dimension of intersection mask. BIDS filename required.

    Required Arguments:
        -d </path/to/dir> = location of project derivatives directory, should contain
            the sub-directory "template".
        -s <session> = BIDS session string
        -t <task> = BIDS task string

    Optional Arguments:
        <list of seeds> = extra args can be supplied as a list of seed files
            (output of func1_amgPriors.sh). These will be multiplied by
            group intersection mask.

    Example Usage:
        deriv_dir=/home/data/madlab/McMakin_EMUR01/derivatives/emu_unc
        seed_dir=\${deriv_dir}/template
        sbatch priors_intx_mask.sh \\
            -d \$deriv_dir \\
            -s ses-S1 \\
            -t task-study \\
            \${seed_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_desc-amgL_mask.nii.gz \\
            \${seed_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_desc-amgR_mask.nii.gz

USAGE
}

# receive args
while getopts ":d:s:t:h" OPT; do
    case $OPT in
    d)
        deriv_dir=${OPTARG}
        if [ ! -d $deriv_dir ] || [ ! -d ${deriv_dir}/template ]; then
            echo -e "\n\t ERROR: did not detect $deriv_dir." >&2
            Usage
            exit 1
        fi
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

# get remaning args as seed list
shift "$((OPTIND - 1))"
seed_list=("$@")
if [ ${#seed_list[@]} -eq 0 ]; then
    unset seed_list && seed_list=("None")
fi

# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
    deriv_dir)
        h_ret="-d"
        ;;
    sess)
        h_ret="-s"
        ;;
    task)
        h_ret="-t"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in deriv_dir sess task; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# print report
cat <<-EOF

    Checks passed, options captured:
        -d : $deriv_dir
        -s : $sess
        -t : $task
        seed : ${seed_list[@]}

EOF

# start command
intx_mask=${deriv_dir}/template/tpl-MNIPediatricAsym_cohort-5_res-2_${sess}_${task}_desc-grpIntx_mask.nii.gz
maskCmd=(3dmask_tool
    -frac 1
    -prefix $intx_mask
    -input
)

# find subjs with epi-anat intx mask
subj_list_all=($(ls $deriv_dir | grep "sub-*"))
for subj in ${subj_list_all[@]}; do
    mask_file=${deriv_dir}/${subj}/${sess}/anat/${subj}_${sess}_${task}_space-MNIPediatricAsym_cohort-5_res-2_desc-intersect_mask.nii.gz
    [ -f $mask_file ] && maskCmd+=($mask_file)
done

# make group intersection mask
if [ ! -f $intx_mask ]; then
    echo -e "Starting:\n\t${maskCmd[@]}"
    "${maskCmd[@]}"
fi

# multiply seed mask by intersection mask
[ ${seed_list[0]} == "None" ] && exit 0
for seed in ${seed_list[@]}; do

    # separate $seed into path, file
    seed_file="$(basename $seed)"
    seed_dir="$(dirname $seed)"

    # determine input, ref resolution from BIDs string
    seed_res=res-$(echo $seed | grep -o -P '(?<=res-).(?=_)')
    ref_res=res-$(echo $intx_mask | grep -o -P '(?<=res-).(?=_)')

    # setup output string
    seed_new_res="${seed_file/$seed_res/$ref_res}"
    seed_clean=${seed_dir}/"${seed_new_res/_mask.nii.gz/Clean_mask.nii.gz}"
    echo -e "\n\t Making $seed_clean"

    # resample
    tmp_res=${seed_dir}/tmp-res_${seed_new_res}
    3dresample \
        -master $intx_mask \
        -rmode NN \
        -input $seed \
        -prefix $tmp_res

    # multiply and binarize
    tmp_mult=${seed_dir}/tmp-mult_${seed_new_res}
    c3d \
        $intx_mask $tmp_res \
        -multiply \
        -o $tmp_mult

    c3d \
        $tmp_mult \
        -thresh 0.3 1 1 0 \
        -o $seed_clean && rm $tmp_res $tmp_mult
done
