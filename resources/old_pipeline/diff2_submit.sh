#!/bin/bash

function Usage {
    cat <<USAGE

    Submit pyAFQ via CLI on scheduled resources. Conda environment
    emuR01_unc is required.

    Required Arguments:
        -t <config_file> = path to config.toml
        -w <work_dir> = path to derivatives directory

    Example Usage:
        $0 \\
            -w /scratch/madlab/emu_unc \\
            -t /home/nmuncy/compute/emu_unc/docs/config.toml

USAGE
}

# capture arguments
while getopts ":t:w:h" OPT; do
    case $OPT in
    t)
        config_file=${OPTARG}
        if [ ! -f $config_file ] || [ -z $config_file ]; then
            echo -e "\n\t ERROR: $config_file file not found or is empty." >&2
            Usage
            exit 1
        fi
        ;;
    w)
        work_dir=${OPTARG}
        if [ ! -d $work_dir ]; then
            echo -e "\n\t ERROR: $work_dir not detected." >&2
            Usage
            exit 1
        fi
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
    config_file)
        h_ret="-t"
        ;;
    work_dir)
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

for opt in work_dir config_file; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for conda env
which python | grep "emuR01_unc" >/dev/null 2>&1
if [ $? != 0 ]; then
    echo "ERROR: Please conda activate emuR01_unc_env and try again." >&2
    exit 1
fi

cat <<-EOF

    Checks passed, options captured:
        -w : $work_dir
        -t : $config_file

    Starting:
        sbatch diff2_probl_cli.sh \\
            $config_file

EOF

# submit jobs
time=$(date '+%Y-%m-%d_%H:%M')
out_dir=${work_dir}/slurm_out/pyafq_${time}
mkdir -p $out_dir

sbatch \
    -e ${out_dir}/err.txt \
    -o ${out_dir}/out.txt \
    diff2_prob_CLI.sh $config_file
