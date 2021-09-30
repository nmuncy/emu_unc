import os
import fnmatch
from datetime import datetime
import time
import subprocess


def main():

    # set necessary paths and variables
    code_dir = "/home/nmuncy/compute/emu_unc"
    deriv_dir = "/scratch/madlab/emu_UNC/derivatives"
    afni_dir = os.path.join(deriv_dir, "afni")

    task = "test"
    sess = "ses-S2"
    num_runs = 3

    # make slurm out dir, afni dir
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir,
        f"""Slurm_out/afniDcn_{current_time.strftime("%y-%m-%d_%H:%M")}""",
    )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # list of afni subjs
    subj_afni = [x for x in os.listdir(afni_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_afni.sort()

    # those who need deconvolution
    subj_list = []
    for subj in subj_afni:
        if not os.path.exists(
            os.path.join(afni_dir, subj, sess, f"{task}_decon_stats_REML+tlrc.HEAD")
        ):
            subj_list.append(subj)

    if len(subj_list) == 0:
        return

    # submit jobs
    for subj in subj_list[:10]:

        h_out = os.path.join(out_dir, f"out_{subj}.txt")
        h_err = os.path.join(out_dir, f"err_{subj}.txt")
        h_job = f"""dcn{subj.split("-")[1]}"""

        sbatch_job = f"""
            sbatch \
                -J "{h_job}" \
                -t 00:30:00 \
                --mem=4000 \
                --ntasks-per-node=1 \
                -p IB_44C_512G  \
                -o {h_out} -e {h_err} \
                --account iacc_madlab \
                --qos pq_madlab \
                --wrap="~/miniconda3/bin/python {code_dir}/func4_deconvolve.py \
                    {subj} \
                    {sess} \
                    {task} \
                    {num_runs} \
                    {afni_dir}"
        """
        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id.decode("utf-8"))
        time.sleep(1)


if __name__ == "__main__":
    main()
