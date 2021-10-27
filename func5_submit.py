"""Wrap func5_test.py

This is merely a wrapper script, used for submitting
func5_test.py with required options.

Example
-------
python func5_submit.py
"""

import subprocess

sbatch_job = """
    sbatch \
        -J "pMVM" \
        -t 06:00:00 \
        --mem=4000 \
        --ntasks-per-node=1 \
        -p IB_44C_512G  \
        -o out.txt -e err.txt \
        --account iacc_madlab \
        --qos pq_madlab \
        --wrap="~/miniconda3/bin/python func5_test.py \
            -p /scratch/madlab/emu_UNC \
            -t /home/data/madlab/singularity-images/templateflow/tpl-MNIPediatricAsym/cohort-5 \
            -s tpl-MNIPediatricAsym_cohort-5_res-1 \
            -m mvm_plan.json"
"""
sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
job_id = sbatch_submit.communicate()[0]
print(job_id.decode("utf-8"))
