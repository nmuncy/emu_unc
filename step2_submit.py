import os
import fnmatch
from datetime import datetime
import time
import subprocess


def main():

    # set necessary paths and variables
    code_dir = "/home/nmuncy/compute/emu_unc"
    deriv_dir = "/scratch/madlab/emu_UNC/derivatives"
    prep_dir = os.path.join(deriv_dir, "fmriprep")
    afni_dir = os.path.join(deriv_dir, "afni")

    task = "test"
    sess = "ses-S2"
    num_runs = 3
    space = "space-MNIPediatricAsym_cohort-5_res-2"

    # make slurm out dir, afni dir
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir,
        f"""Slurm_out/afniPP_{current_time.strftime("%y-%m-%d_%H:%M")}""",
    )
    for h_dir in [out_dir, afni_dir]:
        if not os.path.exists(h_dir):
            os.makedirs(h_dir)

    # list of fmriprep subjs
    subj_fmriprep = [
        x.split(".")[0] for x in os.listdir(prep_dir) if fnmatch.fnmatch(x, "*.html")
    ]
    subj_fmriprep.sort()

    # those who have fmriprep and need to finish pre-processing
    subj_list = []
    for subj in subj_fmriprep:
        prep_bool = os.path.exists(
            os.path.join(
                prep_dir,
                subj,
                sess,
                "func",
                f"{subj}_{sess}_task-{task}_run-1_{space}_desc-preproc_bold.nii.gz",
            )
        )
        afni_bool = os.path.exists(
            os.path.join(afni_dir, subj, sess, f"run-1_{task}_scale+tlrc.HEAD")
        )
        if prep_bool and not afni_bool:
            subj_list.append(subj)

    if len(subj_list) == 0:
        return

    # submit jobs
    for subj in subj_list:

        h_out = os.path.join(out_dir, f"out_{subj}.txt")
        h_err = os.path.join(out_dir, f"err_{subj}.txt")
        h_job = f"""afni{subj.split("-")[1]}"""

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
                --wrap="~/miniconda3/bin/python {code_dir}/step2_finish_preproc.py \
                    {subj} \
                    {sess} \
                    {task} \
                    {num_runs} \
                    {prep_dir} \
                    {afni_dir}"
        """
        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id.decode("utf-8"))
        time.sleep(1)


if __name__ == "__main__":
    main()
