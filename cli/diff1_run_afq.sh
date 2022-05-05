#!/bin/bash

function Usage {
    cat <<USAGE

    Copy data to a working directory.

    This script will copy data from a BIDS-structured project
    directory to a BIDS-structured working/processing directory.
    This is because we store our project data in one location
    (/home/data/madlab/...) but process our data and manage
    intermediates in another (/scratch/madlab/...).

    Pre-processed dwi data will be copied from
        <data_dir>/derivatives/<deriv_dir> to
        <proc_dir>/derivatives/<deriv_dir>.

    Then, submit pyAFQ.

    Required Arguments:
        -c <code_dir> = path to clone of github.com/nmuncy/emu_unc.git
        -d <data_dir> = location of BIDS-structured stored project data
        -j <json-file> = BIDS dataset_description.json sidecar
        -p <proc_dir> = location to where data will be copied,
            and processed. Will be created if it does not exist
        -r <run-str> = BIDS run string, for naming
        -s <ses-str> = BIDS session string, for data structure and naming
        -t <config_file> = path to config.toml
        -x <deriv_dir> = directory within derivatives containing pre-processed
            diffusion weighted data

    Example Usage:
        code_dir="\$(dirname "\$(pwd)")"
        $0 \\
            -c \$code_dir
            -d /home/data/madlab/McMakin_EMUR01 \\
            -j /home/data/madlab/McMakin_EMUR01/dset/dataset_description.json \\
            -p /scratch/madlab/emu_unc \\
            -r run-1 \\
            -s ses-S2 \\
            -t /home/nmuncy/compute/emu_unc/docs/config.toml \\
            -x dwi_preproc

    Notes:
        In our file structure, bval exists in dset, while
        bvec and nii exist in derivatives/<deriv_dir>.

USAGE
}

# capture arguments
while getopts ":c:d:j:p:r:s:t:x:h" OPT; do
    case $OPT in
    c)
        code_dir=${OPTARG}
        if [ ! -d $code_dir ]; then
            echo -e "\n\t ERROR: $code_dir not detected." >&2
            Usage
            exit 1
        fi
        ;;
    d)
        data_dir=${OPTARG}
        if [ ! -d $data_dir ]; then
            echo -e "\n\t ERROR: $data_dir not found or is not a directory." >&2
            Usage
            exit 1
        fi
        ;;
    j)
        json_file=${OPTARG}
        if [ ! -f $json_file ]; then
            echo -e "\n\t ERROR: $json_file not found." >&2
            Usage
            exit 1
        fi
        ;;
    p)
        proc_dir=${OPTARG}
        ;;
    r)
        run=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    t)
        config_file=${OPTARG}
        if [ ! -f $config_file ] || [ -z $config_file ]; then
            echo -e "\n\t ERROR: $config_file file not found or is empty." >&2
            Usage
            exit 1
        fi
        ;;
    x)
        diff_dir=${OPTARG}
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
    data_dir)
        h_ret="-d"
        ;;
    proc_dir)
        h_ret="-p"
        ;;
    sess)
        h_ret="-s"
        ;;
    run)
        h_ret="-r"
        ;;
    diff_dir)
        h_ret="-x"
        ;;
    json_file)
        h_ret="-j"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in data_dir proc_dir sess run diff_dir json_file; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for derivatives dir
if [ ! -d ${data_dir}/derivatives/$diff_dir ]; then
    echo -e "\n\t ERROR: $diff_dir not found or is not a directory." >&2
    Usage
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, options captured:
        -d : $data_dir
        -p : $proc_dir
        -s : $sess
        -r : $run
        -x : $diff_dir
        -j : $json_file

EOF

# Get oriented
deriv_dir=${data_dir}/derivatives/$diff_dir
dset_dir=${data_dir}/dset
work_dir=${proc_dir}/derivatives/$diff_dir

# Get json
mkdir -p $work_dir
cp $json_file $proc_dir
cp $json_file $work_dir

# Copy, BIDs format pre-processed dwi data
unset subj_list
subj_list=($(ls $deriv_dir | grep "sub-*"))
for subj in ${subj_list[@]}; do

    source_dir=${deriv_dir}/${subj}/${sess}/dwi
    out_dir=${work_dir}/${subj}/${sess}/dwi

    if [ ! -d $out_dir ] || [ ! -f ${out_dir}/${subj}_${sess}_dwi.nii.gz ]; then

        echo -e "\t Copying data for $subj ..."
        mkdir -p $out_dir

        cp ${dset_dir}/${subj}/${sess}/dwi/${subj}_${sess}_${run}_dwi.bval \
            ${out_dir}/${subj}_${sess}_dwi.bval
        cp ${source_dir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.bvec \
            ${out_dir}/${subj}_${sess}_dwi.bvec
        cp ${source_dir}/${subj}_${sess}_${run}_desc-eddyCorrected_dwi.nii.gz \
            ${out_dir}/${subj}_${sess}_dwi.nii.gz

    fi
done

# submit afq
time=$(date '+%Y-%m-%d_%H:%M')
out_dir=${work_dir}/slurm_out/pyafq_${time}
mkdir -p $out_dir

sbatch \
    -e ${out_dir}/err.txt \
    -o ${out_dir}/out.txt \
    ${code_dir}/resources/diff/afq_cli.sh $config_file
