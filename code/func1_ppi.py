"""Title

Desc.
"""

# %%
import os
import fnmatch
import subprocess
import pandas as pd


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
def clean_data(subj, subj_data, subj_out, decon_str):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    # get TR
    h_cmd = f"3dinfo -tr {subj_data}/*_run-1_*_desc-scaled_bold.nii.gz"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_tr = float(h_out.decode("utf-8").strip())

    # get proper brick length
    #   REML appends an extra brick because
    #   "reasons". Account for AFNIs random
    #   0-1 indexing
    h_cmd = f"3dinfo -nv {subj_data}/decon_{decon_str}_cbucket_REML+tlrc"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_wrong = h_out.decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # list all scale files
    scale_list = [
        os.path.join(subj_data, x)
        for x in os.listdir(subj_data)
        if fnmatch.fnmatch(x, "*desc-scaled_bold.nii.gz")
    ]
    scale_list.sort()

    # make clean data
    clean_file = os.path.join(subj_out, "CleanData+tlrc.HEAD")
    if not os.path.exists(clean_file):

        #  list undesirable sub-bricks (those starting with Run or mot)
        no_int = []
        with open(os.path.join(subj_data, f"X.decon_{decon_str}.xmat.1D")) as f:
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
            cd {subj_out}

            3dTcat \
                -prefix tmp_cbucket \
                -tr {len_tr} \
                {subj_data}/"decon_{decon_str}_cbucket_REML+tlrc[0..{len_right}]"

            3dSynthesize \
                -prefix tmp_effNoInt \
                -matrix {subj_data}/X.decon_{decon_str}.xmat.1D \
                -cbucket tmp_cbucket+tlrc \
                -select {" ".join(no_int)} \
                -cenfill nbhr

            3dTcat \
                -prefix tmp_all_runs \
                -tr {len_tr} \
                {" ".join(scale_list)}

            3dcalc \
                -a tmp_all_runs+tlrc \
                -b tmp_effNoInt+tlrc \
                -expr 'a-b' \
                -prefix CleanData && rm tmp_*
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}cle", subj_out)

    assert os.path.exists(clean_file), "CleanData not found."
    file_dict = {
        "TR Length": len_tr,
        "Clean Data": clean_file,
        "Scale Data": scale_list,
    }
    return file_dict


# %%
def func_upsample(timeseries, resolution, tog_inter):
    """
    Receives a numeric list and returns and
        list upsampled by factor = 1000

    Will interpolate or repeat according
        to tog_inter (1 = interpolate)

    Paremeters
    ----------

    Returns
    -------

    """
    df = pd.DataFrame(
        data=timeseries, index=pd.RangeIndex(len(timeseries)) * resolution
    )
    df = df.rename(columns={0: "TS"})
    df.index = df.index.map(datetime.datetime.utcfromtimestamp)
    df = df.reindex(pd.date_range(df.index.min(), df.index.max(), freq="1L"))
    if tog_inter == 1:
        df.interpolate(method="polynomial", order=7, inplace=True)
    else:
        df = df.ffill()
    out_list = df["TS"].tolist()
    return out_list


def seed_timeseries(subj, subj_out, file_dict, seed_dict, dur=2):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    # find smallest multiplier that returns int for resampling
    res_multiplier = 2
    status = True
    len_tr = file_dict["TR Length"]

    while status:
        if ((len_tr * res_multiplier) % 2) == 0:
            status = False
        else:
            res_multiplier += 1

    # set status
    resample_decision = "easy" if res_multiplier <= 32 else "hard"

    # make ideal HRF, use same model as deconvolution
    hrf_file = os.path.join(subj_out, "HRF_model.1D")
    if not os.path.exists(hrf_file):
        h_cmd = f"""
            cd {subj_out}

            3dDeconvolve\
                -polort -1 \
                -nodata {round((1 / len_tr) * 19)} {len_tr} \
                -num_stimts 1 \
                -stim_times 1 1D:0 'TWOGAMpw(4,5,0.2,12,7,{dur})' \
                -x1D HRF_model.1D \
                -x1D_stop
        """
        h_out, h_err = submit_hpc_subprocess(h_cmd)
    assert os.path.exists(hrf_file), "HRF model failed."

    clean_data = file_dict["Clean Data"]
    for seed, coord in seed_dict.items():

        if not os.path.exists(os.path.join(subj_out, f"Seed_{seed}_neural_us.1D")):

            # make seed, get timeseries
            h_cmd = f"""
                cd {subj_out}

                echo {coord} | 3dUndump \
                    -xyz \
                    -srad 3 \
                    -master {clean_data} \
                    -prefix Seed_{seed} -

                3dmaskave \
                    -quiet \
                    -mask Seed_{seed}+tlrc \
                    {clean_data} > Seed_{seed}_orig.1D

                3dTfitter \
                    -RHS Seed_{seed}_orig.1D \
                    -FALTUNG {hrf_file} \
                    tmp.1D 012 0

                1dtranspose \
                    tmp.1D > Seed_{seed}_neural.1D
            """
            subj_num = subj.split("-")[1]
            h_out, h_err = submit_hpc_sbatch(
                h_cmd, 1, 4, 1, f"{subj_num}mksd", subj_out
            )

            # upsample seed timeseries
            if resample_decision == "easy":
                h_cmd = f"""
                    1dUpsample \
                        {res_multiplier} \
                        {subj_out}/Seed_{seed}_neural.1D \
                        > {subj_out}/Seed_{seed}_neural_us.1D
                """
                h_out, h_err = submit_hpc_subprocess(h_cmd)

            elif resample_decision == "hard":

                seed_file = open(f"{subj_out}/Seed_{seed}_neural.1D", "r")
                h_ts = seed_file.readlines()
                seed_ts = [float(x.strip()) for x in h_ts]
                seed_us = func_upsample(seed_ts, len_tr, 1)

                with open(f"{subj_out}/Seed_{seed}_neural_us.1D", "w") as seed_out:
                    for value in seed_us:
                        seed_out.write("%s\n" % value)


# %%
def main():

    # For testing
    data_dir = "/home/data/madlab/McMakin_EMUR01/derivatives/afni"
    deriv_dir = "/scratch/madlab/emu_unc/derivatives/afni_ppi"
    subj = "sub-4001"
    sess = "ses-S2"
    decon_str = "task-test_UniqueBehs"

    # make dict for seed construction from coordinates
    seed_dict = {"LHC": "-24 -12 -22"}

    subj_data = os.path.join(data_dir, subj, sess, "func")
    subj_out = os.path.join(deriv_dir, subj, sess, "func")
    if not os.path.exists(subj_out):
        os.makedirs(subj_out)

    file_dict = clean_data(subj, subj_data, subj_out, decon_str)
    seed_timeseries(subj, subj_out, file_dict, seed_dict)


if __name__ == "__main__":
    main()
