#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper script for func1_amgPriors.sh, which uses FreeSurfer posteriors
    to create template priors via Joint Label Fusion. Stderr/out of
    func1_amgPriors.sh captured in <outDir>/Slurm_out/amgJLF_<timestamp>.

    Usage:
        $0 --run y

    Optional Arguments:
        --run = [y/N] whether to run or print this help
        -o <outDir> = path to output directory
            default: /scratch/madlab/McMakin_EMUR01/derivatives/emu_unc/template
        -f <fsDir> = path to FreeSurfer derivatives
            default: /home/data/madlab/McMakin_EMUR01/derivatives/freesurfer
        -t <template> = path to template
            default: \${outDir}/tpl-MNIPediatricAsym_cohort-5_res-1_T1w.nii.gz

USAGE
}

if [ $# == 0 ]; then
    Usage
    exit 0
fi

# set defaults
run=n
outDir=/scratch/madlab/McMakin_EMUR01/derivatives/emu_unc/template
fsDir=/home/data/madlab/McMakin_EMUR01/derivatives/freesurfer
atlas=${outDir}/tpl-MNIPediatricAsym_cohort-5_res-1_T1w.nii.gz

# Check options
while (($# >= 1)); do
    case $1 in
    -f)
        fsDir=$2
        ;;
    -o)
        outDir=$2
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

# check for atlas, fsDir
if [ ! -f $atlas ]; then
    echo -e "\n\n ERROR: option \"-t\" template not detected." >&2
    Usage
    exit 1
fi

if [ ! -d $fsDir ]; then
    echo -e "\n\n ERROR: option \"-f\" directory not detected." >&2
    Usage
    exit 1
fi

# set up
codeDir=$(pwd)
time=$(date '+%Y-%m-%d_%H:%M')
slurmDir=${outDir}/Slurm_out/amgJLF_$time

# print feedback
cat <<-EOF

    ${codeDir}/func1_amgPriors.sh will use
    following parameters (use "--run y" to start job):
        -o <outDir> = $outDir
        -f <fsDir> = $fsDir
        -t <template> = $atlas

EOF

# submit job
if [ $run == "y" ]; then
    echo -e "Detected \$run == \"y\", starting ...\n"
    mkdir -p $slurmDir
    sbatch \
        -e ${slurmDir}/err_amgJLF.txt \
        -o ${slurmDir}/out_amgJLF.txt \
        ${codeDir}/func1_amgPriors.sh $outDir $fsDir $atlas
fi
