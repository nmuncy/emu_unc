#!/bin/bash

# check that all jobs are done
num_jobs=`squeue -u $(whoami) | wc -l`
if [ $num_jobs -gt 1 ]; then
    echo "Jobs still running, exiting ..."
    exit 0
fi

# update timing files
#source step3_submit.sh

# start decon script
~/miniconda3/bin/python step4_submit.py
