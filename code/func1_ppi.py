"""Title

Desc.
"""

# %%
import os
import fnmatch
import subprocess


# %%
def submit_hpc_subprocess(bash_command):
    """Submit quick job as subprocess.

    Run a quick job as a subprocess and capture stdout/err. Use
    for AFNI and c3d commands.

    Parameters
    ----------
    bash_command : str
        Bash syntax, to be executed

    Returns
    -------
    (job_out, job_err) : tuple of str
        job_out = subprocess stdout
        job_out = subprocess stderr

    Example
    -------
    submit_hpc_subprocess("afni -ver")

    """
    h_cmd = f"""
        module load afni-20.2.06
        module load c3d-1.0.0-gcc-8.2.0
        {bash_command}
    """
    h_sp = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    job_out, job_err = h_sp.communicate()
    h_sp.wait()
    return (job_out, job_err)


def submit_hpc_sbatch(command, wall_hours, mem_gig, num_proc, job_name, out_dir):
    """Submit job to slurm scheduler (sbatch).

    Sbatch submit a larger job with scheduled resources. Waits for
    job_name to no longer be found in squeue. Stderr/out written to
    out_dir/sbatch_<job_name>.err/out. Supports AFNI and c3d commands.

    Parameters
    ----------
    command : str
        Bash code to be scheduled
    wall_hours : int
        number of desired walltime hours
    mem_gig : int
        amount of desired RAM
    num_proc : int
        number of desired processors
    job_name : str
        job name
    out_dir : str
        location for <job_name>.err/out

    Returns
    -------
    (job_name, job_id) : tuple of str
        job_name = scheduled job name
        job_id = scheduled job ID

    Example
    -------
    submit_hpc_sbatch("afni -ver")
    """
    sbatch_job = f"""
        sbatch \
        -J {job_name} \
        -t {wall_hours}:00:00 \
        --cpus-per-task={num_proc} \
        --mem-per-cpu={mem_gig}000 \
        -p IB_44C_512G \
        -o {out_dir}/{job_name}.out \
        -e {out_dir}/{job_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wait \
        --wrap="module load afni-20.2.06
            module load c3d-1.0.0-gcc-8.2.0
            {command}
        "
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_response.communicate()[0].decode("utf-8")
    return (job_name, job_id)


# %%
def clean_data(subj_dir, decon_str):

    # get TR
    h_cmd = f"3dinfo -tr {subj_dir}/*_run-1_*_desc-scaled_bold.nii.gz"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_tr = float(h_out.decode("utf-8").strip())

    # get proper brick length
    #   REML appends an extra brick because
    #   "reasons". Account for AFNIs random
    #   0-1 indexing
    h_cmd = f"3dinfo -nv {subj_dir}/decon_{decon_str}_stats_REML+tlrc"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_wrong = h_out.decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # list all scale files
    scale_list = [
        x
        for x in os.listdir(subj_dir)
        if fnmatch.fnmatch(x, "*desc-scaled_bold.nii.gz")
    ]
    scale_list.sort()

    # list undesirable sub-bricks (those starting with Run or mot)
    no_int = []
    with open(os.path.join(subj_dir, f"X.decon_{decon_str}.xmat.1D")) as f:
        h_file = f.readlines()
        for line in h_file:
            if line.__contains__("ColumnLabels"):
                col_list = (
                    line.replace("#", "").split('"')[1].replace(" ", "").split(";")
                )
                for i, j in enumerate(col_list):
                    if fnmatch.fnmatch(j, "Run*") or fnmatch.fnmatch(j, "mot*"):
                        no_int.append(f"{str(i)}")

    # strip extra sub-brick, make clean data by removing
    #   effects of no interest from concatenated runs
    h_cmd = f"""
        cd {subj_dir}
        3dTcat -prefix tmp_cbucket -tr {len_tr} "decon_{decon_str}_cbucket_REML+tlrc[0..{len_right}]"
        3dSynthesize -prefix tmp_effNoInt -matrix X.decon_{decon_str}.xmat.1D \
            -cbucket tmp_cbucket+tlrc -select {" ".join(no_int)} -cenfill nbhr
        3dTcat -prefix tmp_all_runs -tr {len_tr} {" ".join(scale_list)}
        3dcalc -a tmp_all_runs+tlrc -b tmp_effNoInt+tlrc -expr 'a-b' -prefix CleanData    """


# %%
def main():

    # For testing
    data_dir = "/home/data/madlab/McMakin_EMUR01/derivatives/afni"
    deriv_dir = "/scratch/madlab/emu_unc/derivatives/afni_ppi"
    subj = "sub-4001"
    sess = "ses-S2"
    decon_str = "task-test_UniqueBehs"

    subj_dir = os.path.join(data_dir, subj, sess, "func")
    # clean_data(subj_dir, decon_str)


if __name__ == "__main__":
    main()
