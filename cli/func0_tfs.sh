#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper for resources.func.tf_mk.R. Makes timing files for each
    subject in <project_dir>.

    Required Arguments:
        -c <code_dir> = path to clone of github.com/nmuncy/emu_unc.git
        -n <decon_name> = deconvolution identifier string
            "rVal" for decon-rVal_*
        -p <project_dir> = BIDS project directory
        -r <hpc_dir> = path to hpc timing file location, written
            into the json sidecar by resources.func.tf_mk_json.py for
            use with github.com/emu-project/func_processing
            resources.afni.deconvolve.write_new_decon
        -s <session> = BIDS session string of study session
        -t <session> = BIDS session string of test session
        -w <local_dir> = path to local output dir for timing files

    Example Usage:
        code_dir="\$(dirname "\$(pwd)")"
        $0 \\
            -c \$code_dir \\
            -n rVal \\
            -p /Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset \\
            -r /home/nmuncy/compute/emu_unc/data/timing_files \\
            -s ses-S1 \\
            -t ses-S2 \\
            -w \${code_dir}/data/timing_files

USAGE
}

# capture arguments
while getopts ":c:n:p:r:s:t:w:h" OPT; do
    case $OPT in
    c)
        code_dir=${OPTARG}
        if [ ! -d $code_dir ]; then
            echo -e "\n\t ERROR: $code_dir not detected." >&2
            Usage
            exit 1
        fi
        ;;
    n)
        decon_name=${OPTARG}
        ;;
    p)
        proj_dir=${OPTARG}
        if [ ! -d $proj_dir ]; then
            echo -e "\n\t ERROR: -p directory not found.\n" >&2
            Usage
            exit 1
        fi
        ;;
    r)
        hpc_path=${OPTARG}
        ;;
    s)
        study=${OPTARG}
        ;;
    t)
        test=${OPTARG}
        ;;
    w)
        timing_dir=${OPTARG}
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
    decon_name)
        h_ret="-n"
        ;;
    proj_dir)
        h_ret="-p"
        ;;
    hpc_path)
        h_ret="-r"
        ;;
    study)
        h_ret="-s"
        ;;
    test)
        h_ret="-t"
        ;;
    timing_dir)
        h_ret="-w"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in code_dir decon_name proj_dir hpc_path study test timing_dir; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check python
python --version | grep "3." >/dev/null 2>&1
if [ $? != 0 ]; then
    echo "ERROR: Please use python3." >&2
    exit 1
fi

# make timing files for each subject in $proj_dir
mkdir -p $timing_dir
subj_list=($(ls $proj_dir | grep --exclude=*.json "sub-*"))
for subj in ${subj_list[@]}; do

    # set output path
    subj_out=${timing_dir}/${subj}/$study

    # set subject paths, find tsv files for study, test
    subj_study=${proj_dir}/${subj}/${study}/func
    study_list=($(ls ${subj_study}/*_task-study_*events.tsv))
    subj_test=${proj_dir}/${subj}/${test}/func
    test_list=($(ls ${subj_test}/*_task-test_*events.tsv))

    # require 2 study and 3 test run files
    if [ ${#study_list[@]} != 2 ] || [ ${#test_list[@]} != 3 ]; then
        continue
    fi

    # submit R work
    echo "Making TFs for $subj ..."
    mkdir -p $subj_out
    Rscript \
        ${code_dir}/resources/func/tf_mk.R \
        $proj_dir \
        $subj \
        $subj_out \
        "${study_list[@]}" \
        "${test_list[@]}"
done

if [ $? != 0 ]; then
    echo -e "\n\tERROR: Issue making timing files, check reources.func.tf_mk.R" >&2
    exit 1
fi

# generate json sidecars
python \
    ${code_dir}/resources/func/tf_mk_json.py \
    -s $study \
    -w $timing_dir \
    -p $hpc_path \
    -n $decon_name
