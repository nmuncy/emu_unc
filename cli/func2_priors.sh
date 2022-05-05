#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper for resources.afni.amg_priors.sh, which uses FreeSurfer posteriors
    to create template priors via Joint Label Fusion.

    Required Arguments:
        -c <code_dir> = path to clone of github.com/nmuncy/emu_unc.git
        -f <fs_dir> = path to FreeSurfer derivatives
        -o <out_dir> = path to output directory
        -r <template> = path to template
        -s <session> = BIDS session string
        -t <task> = BIDS task string


    Example Usage:
        code_dir="\$(dirname "\$(pwd)")"
        out_dir=/scratch/madlab/McMakin_EMUR01/derivatives/emu_unc/template
        $0 \\
            -c \$code_dir \\
            -f /home/data/madlab/McMakin_EMUR01/derivatives/freesurfer \\
            -o \$out_dir \\
            -r \${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_T1w.nii.gz \\
            -s ses-S1 \\
            -t task-study

USAGE
}

# capture arguments
while getopts ":c:f:o:r:s:t:h" OPT; do
    case $OPT in
    c)
        code_dir=${OPTARG}
        if [ ! -d $code_dir ]; then
            echo -e "\n\t ERROR: $code_dir not detected." >&2
            Usage
            exit 1
        fi
        ;;
    f)
        fs_dir=${OPTARG}
        if [ ! -d $fs_dir ]; then
            echo -e "\n\t ERROR: -f directory not found.\n" >&2
            Usage
            exit 1
        fi
        ;;
    o)
        out_dir=${OPTARG}
        ;;
    r)
        atlas=${OPTARG}
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
    \?)
        echo -e "\n\t Error: invalid option -${OPTARG}." >&2
        Usage
        exit 1
        ;;
    :)
        echo -e "\n\t Error: -${OPTARG} requires an argument." >&2
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
    code_dir)
        h_ret="-c"
        ;;
    fs_dir)
        h_ret="-f"
        ;;
    out_dir)
        h_ret="-o"
        ;;
    atlas)
        h_ret="-r"
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

for opt in code_dir fs_dir out_dir atlas sess task; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# make amg, intx scripts
cmd_amg=(
    ${code_dir}/resources/func/priors_amg.sh
    $out_dir
    $fs_dir
    $atlas
)

deriv_dir=$(dirname $out_dir)
cmd_intx=(
    ${code_dir}/resources/func/priors_intx_mask.sh
    -d $deriv_dir
    -s $sess
    -t $task
    ${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_desc-amgL_mask.nii.gz
    ${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_desc-amgR_mask.nii.gz
)

# print feedback
cat <<-EOF

    Making amygdala priors via:
        ${cmd_amg[@]}

    and, making intersection priors via:
        ${cmd_intx[@]}

EOF

# submit job
time=$(date '+%Y-%m-%d_%H:%M')
slurm_dir=${out_dir}/Slurm_out/amgJLF_$time
mkdir -p $slurm_dir
sbatch \
    -e ${slurm_dir}/err_amgJLF.txt \
    -o ${slurm_dir}/out_amgJLF.txt \
    "${cmd_amg[@]}"

# get job ID for waiting
wait_amg=$(squeue -u $(whoami) | grep "amgJLF" | awk '{print $1}')

# submit job, wait for previous job
time=$(date '+%Y-%m-%d_%H:%M')
slurm_dir=${out_dir}/Slurm_out/amgIntx_$time
mkdir -p $slurm_dir
sbatch \
    --dependency=afterany:$wait_amg \
    -e ${slurm_dir}/err_intx.txt \
    -o ${slurm_dir}/out_intx.txt \
    "${cmd_intx[@]}"
