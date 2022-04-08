#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper for foo_studyTF.R. Makes timing files for each
    subject in <project_dir>.

    Required Arguments:
        -p <project_dir> = BIDS project directory

    Example Usage:
        $0 \\
            -p /Volumes/homes/MaDLab/projects/McMakin_EMUR01/dset

USAGE
}

# capture arguments
while getopts ":p:h" OPT; do
    case $OPT in
    p)
        proj_dir=${OPTARG}
        if [ ! -d $proj_dir ]; then
            echo -e "\n\t ERROR: -p directory not found.\n" >&2
            Usage
            exit 1
        fi
        ;;
    h)
        Usage
        exit 0
        ;;
    \?)
        echo -e "\n\t Error: invalid option -${OPTARG}."
        Usage
        exit 1
        ;;
    :)
        echo -e "\n\t Error: -${OPTARG} requires an argument."
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

# set up - resolve path to data_dir, make timing_files dir
data_dir=$(builtin cd ../data; pwd)
out_dir=${data_dir}/timing_files
mkdir $out_dir

# make subject list
subj_list=($(ls $proj_dir | grep "sub-*"))
for subj in ${subj_list[@]}; do

    # set output path
    subj_out=${out_dir}/${subj}/ses-S1

    # set subject paths, find tsv files for study, test
    subj_study=${proj_dir}/${subj}/ses-S1/func
    study_list=($(ls ${subj_study}/*_task-study_*events.tsv))

    subj_test=${proj_dir}/${subj}/ses-S2/func
    test_list=($(ls ${subj_test}/*_task-test_*events.tsv))

    # require 2 study and 3 test run files
    if [ ${#study_list[@]} != 2 ] || [ ${#test_list[@]} != 3 ]; then
        continue
    fi

    # submit R work
    echo "Making TFs for $subj ..."
    mkdir -p $subj_out
    Rscript funcX_studyTF.R \
        $proj_dir \
        $subj \
        $subj_out \
        "${study_list[@]}" \
        "${test_list[@]}"
done
