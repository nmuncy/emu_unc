#!/bin/bash

function Usage {
    cat << USAGE

    Submit pyAFQ via CLI on scheduled resources. Conda environment
    emuR01_unc is required.

    Required Arguments:
        -c <code_dir> = path to code directory, containing diff2_prob_CLI.sh
        -t <config_file> = path to config.toml
        -w <work_dir> = path to derivatives directory

    Example Usage:
        ./diff2_submit.sh \\
            -c /home/nmuncy/compute/emu_unc/code \\
            -w /scratch/madlab/emu_unc \\
            -t /home/nmuncy/compute/emu_unc/docs/config.toml

USAGE
}


# capture arguments
unset code_dir config_file work_dir
while getopts ":c:t:w:h" OPT
    do
    case $OPT in
        c) code_dir=${OPTARG}
            ;;
        t) config_file=${OPTARG}
            ;;
        w) work_dir=${OPTARG}
            ;;
        h)
            Usage
            exit 0
            ;;
        \?) echo -e "\n \t ERROR: invalid option." >&2
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

# check required directories
if [ ! -d $code_dir ] || [ -z $code_dir ]; then
    echo -e "\n\t ERROR: \"-c\" directory not detected." >&2
    Usage
    exit 1
fi

if [ ! -d $work_dir ] || [ -z $work_dir ]; then
    echo -e "\n\t ERROR: \"-w\" directory not detected." >&2
    Usage
    exit 1
fi

if [ ! -f $config_file ] || [ -z $config_file ]; then
    echo -e "\n\t ERROR: \"-t\" config.toml file not detected." >&2
    Usage
    exit 1
fi


# check for conda env
which python | grep "emuR01_unc"
if [ $? != 0 ]; then
    echo "ERROR: Please conda activate emuR01_unc_env and try again." >&2
    exit 1
fi

# submit jobs
time=`date '+%Y-%m-%d_%H:%M'`
out_dir=${work_dir}/slurm_out/pyafq_${time}
mkdir -p $out_dir

sbatch \
    -e ${out_dir}/err.txt \
    -o ${out_dir}/out.txt \
    ${code_dir}/diff2_prob_CLI.sh $config_file
