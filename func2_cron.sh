#!/bin/bash

# check that all jobs are done
num_jobs=`squeue -u $(whoami) | wc -l`
if [ $num_jobs -gt 1 ]; then
    echo "Jobs still running, exiting ..."
    exit 0
fi

# run wrapper
~/miniconda3/bin/python func2_submit.py