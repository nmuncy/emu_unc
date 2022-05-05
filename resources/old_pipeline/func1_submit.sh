#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper script for func1_amgPriors.sh, which uses FreeSurfer posteriors
    to create template priors via Joint Label Fusion. Stderr/out of
    func1_amgPriors.sh captured in <out_dir>/Slurm_out/amgJLF_<timestamp>.

    Usage:
        $0 --run y

    Optional Arguments:
        --run = [y/N] whether to run or print this help
        -o <out_dir> = path to output directory
            default: /scratch/madlab/McMakin_EMUR01/derivatives/emu_unc/template
        -f <fs_dir> = path to FreeSurfer derivatives
            default: /home/data/madlab/McMakin_EMUR01/derivatives/freesurfer
        -t <template> = path to template
            default: \${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_T1w.nii.gz

USAGE
}

if [ $# == 0 ]; then
    Usage
    exit 0
fi

# set defaults
run=n
out_dir=/scratch/madlab/McMakin_EMUR01/derivatives/emu_unc/template
fs_dir=/home/data/madlab/McMakin_EMUR01/derivatives/freesurfer
atlas=${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-1_T1w.nii.gz

# Check options
while (($# >= 1)); do
    case $1 in
    -f)
        fs_dir=$2
        ;;
    -o)
        out_dir=$2
        ;;
    -t)
        atlas=$2
        ;;
    --run)
        run=$2
        ;;
    \?)
        echo -e "\n\n ERROR: invalid option." >&2
        Usage
        exit 1
        ;;
    esac
    shift 2

    # catch failure to specify option/parameter (avoid infinite shift)
    num_arg=$#
    if [ ! $(($num_arg % 2)) -eq 0 ]; then
        echo -e "\n\n ERROR: Failed to specify an option/parameter (number of args is not even)." >&2
        Usage
        exit 1
    fi
done

# check for atlas, fs_dir
if [ ! -f $atlas ]; then
    echo -e "\n\n ERROR: option \"-t\" template not detected." >&2
    Usage
    exit 1
fi

if [ ! -d $fs_dir ]; then
    echo -e "\n\n ERROR: option \"-f\" directory not detected." >&2
    Usage
    exit 1
fi

# set up
codeDir=$(pwd)
time=$(date '+%Y-%m-%d_%H:%M')
slurmDir=${out_dir}/Slurm_out/amgJLF_$time

# print feedback
cat <<-EOF

    ${codeDir}/func1_amgPriors.sh will use
    following parameters (use "--run y" to start job):
        -o <out_dir> = $out_dir
        -f <fs_dir> = $fs_dir
        -t <template> = $atlas

EOF

# submit job
if [ $run == "y" ]; then
    echo -e "Detected \$run == \"y\", starting ...\n"
    mkdir -p $slurmDir
    sbatch \
        -e ${slurmDir}/err_amgJLF.txt \
        -o ${slurmDir}/out_amgJLF.txt \
        ${codeDir}/func1_amgPriors.sh $out_dir $fs_dir $atlas
fi
