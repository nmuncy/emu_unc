#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper for func0_masks.py. Checks for which subjects have
    kmeans amygdala masks and not a resampled BLA mask in template
    space. Submits sbatch jobs for N subjects meeting checks.

    Required Arguments:
        -d <project_derivatives> = path to project derivatives location
        -w <scratch_directory> = path to parent scratch directory
        -s <session> = BIDS session string
        -n <number> = number of subjects to submit jobs

    Example Usage:
        $0 \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives \\
            -w /scratch/madlab/emu_unc \\
            -s ses-S2 \\
            -n 8
USAGE
}

# capture arguments
while getopts ":d:f:n:s:w:h" OPT; do
    case $OPT in
    d)
        deriv_dir=${OPTARG}
        if [ ! -d $deriv_dir ]; then
            echo -e "\n\t ERROR: -d directory not found." >&2
            Usage
            exit 1
        fi
        ;;
    n)
        num_subj=${OPTARG}
        if [ $num_subj -lt 1 ]; then
            echo -e "\n\t ERROR: -n arg must be greater than 0." >&2
            Usage
            exit 1
        fi
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

# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
    deriv_dir)
        h_ret="-d"
        ;;
    num_subj)
        h_ret="-n"
        ;;
    sess)
        h_ret="-s"
        ;;
    scratch_dir)
        h_ret="-w"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n\t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in deriv_dir num_subj sess scratch_dir; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for conda env
which python | grep "emuR01_unc" >/dev/null 2>&1
if [ $? != 0 ]; then
    echo -e "\n\t ERROR: Please conda activate emuR01_unc_env and try again.\n" >&2
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, options captured:
        -w : $scratch_dir
        -d : $deriv_dir
        -s : $sess
        -n : $num_subj

EOF

# find subjects missing kmeans masks
echo -e "Building subject lists ...\n"
subj_list=()
kmeans_dir=${deriv_dir}/kmeans
afni_dir=${deriv_dir}/afni
subj_all=($(ls $kmeans_dir | grep "sub-*"))
for subj in ${subj_all[@]}; do
    mask_file=${afni_dir}/${subj}/${sess}/dwi/${subj}_${sess}_space-MNIPediatricAsym_cohort-5_res-2_desc-bla_mask.nii.gz
    if [ ! -f $mask_file ]; then
        subj_list+=($subj)
    fi
done

# submit N jobs
time=$(date '+%Y-%m-%d_%H:%M')
out_dir=${scratch_dir}/slurm_out/kmean_${time}
mkdir -p $out_dir

cat <<-EOF
    Submitting sbatch job

        func0_masks.py \\
            -s <subj>

    with the following subjects:

        ${subj_list[@]:0:$num_subj}

EOF

c=0
while [ $c -lt $num_subj ]; do

    subj=${subj_list[$c]}
    sbatch \
        --job-name=p${subj#*-} \
        --output=${out_dir}/${subj}.out \
        --mem-per-cpu=4000 \
        --partition=IB_44C_512G \
        --account=iacc_madlab \
        --qos=pq_madlab \
        func0_masks.py \
        -s $subj

    let c+=1
    sleep 1
done
