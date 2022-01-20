"""Title

Desc.

Examples
--------
sbatch --job-name=p1234 \\
    --output=p1234 \\
    --mem-per-cpu=4000 \\
    --partition=IB_44C_512G \\
    --account=iacc_madlab \\
    --qos=pq_madlab \\
    func1_ppi.py \\
    -s sub-1234 \\
    -d decon_task-test_UniqueBehs
"""

# %%
import os
import sys
import glob
import fnmatch
import shutil
import subprocess
import pandas as pd
import textwrap
from argparse import ArgumentParser, RawTextHelpFormatter


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


def clean_data(subj, subj_out, afni_data):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    # get TR
    h_cmd = f"3dinfo -tr {afni_data['scaled_files'][0]}"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_tr = float(h_out.decode("utf-8").strip())

    # get files from afni_data
    decon_file = afni_data["decon_file"].split(".")[0]
    decon_matrix = afni_data["decon_matrix"]
    scaled_list = afni_data["scaled_files"]

    # Get proper brick length - REML appends an extra brick because of
    # "reasons". Account for AFNIs random 0-1 indexing
    h_cmd = f"3dinfo -nv {decon_file}"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    len_wrong = h_out.decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # make clean_data
    clean_file = os.path.join(subj_out, "CleanData+tlrc.HEAD")
    if not os.path.exists(clean_file):
        print(f"Making clean_data for {subj}")

        #  list undesirable sub-bricks (those starting with Run or mot)
        no_int = []
        with open(decon_matrix) as f:
            h_file = f.readlines()
            for line in h_file:
                if line.__contains__("ColumnLabels"):
                    col_list = (
                        line.replace("#", "").split('"')[1].replace(" ", "").split(";")
                    )
                    for i, j in enumerate(col_list):
                        if fnmatch.fnmatch(j, "Run*") or fnmatch.fnmatch(j, "mot*"):
                            no_int.append(f"{str(i)}")

        # strip extra sub-brick, make clean_data by removing
        #   effects of no interest from concatenated runs
        h_cmd = f"""
            cd {subj_out}

            3dTcat \
                -prefix tmp_cbucket \
                -tr {len_tr} \
                {decon_file}"[0..{len_right}]"

            3dSynthesize \
                -prefix tmp_effNoInt \
                -matrix {decon_matrix} \
                -cbucket tmp_cbucket+tlrc \
                -select {" ".join(no_int)} \
                -cenfill nbhr

            3dTcat \
                -prefix tmp_all_runs \
                -tr {len_tr} \
                {" ".join(scaled_list)}

            3dcalc \
                -a tmp_all_runs+tlrc \
                -b tmp_effNoInt+tlrc \
                -expr 'a-b' \
                -prefix CleanData && rm tmp_*
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}cle", subj_out)

    assert os.path.exists(clean_file), "CleanData not found."
    afni_data["len_tr"] = len_tr
    afni_data["clean_data"] = clean_file
    return afni_data


def hrf_model(subj_out, afni_data, dur=2):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    # make ideal HRF, use same model as deconvolution
    len_tr = afni_data["len_tr"]
    hrf_file = os.path.join(subj_out, "HRF_model.1D")
    if not os.path.exists(hrf_file):
        print("\nMaking ideal HRF")
        h_cmd = f"""
            3dDeconvolve\
                -polort -1 \
                -nodata {round((1 / len_tr) * 19)} {len_tr} \
                -num_stimts 1 \
                -stim_times 1 1D:0 'TWOGAMpw(4,5,0.2,12,7,{dur})' \
                -x1D {hrf_file} \
                -x1D_stop
        """
        h_out, h_err = submit_hpc_subprocess(h_cmd)
    assert os.path.exists(hrf_file), "HRF model failed."

    afni_data["hrf_model"] = hrf_file
    return afni_data


def seed_timeseries(subj, subj_out, afni_data, seed_coord):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    hrf_file = afni_data["hrf_model"]
    clean_data = afni_data["clean_data"]

    seed_dict = {}
    for seed, coord in seed_coord.items():
        seed_dict[seed] = {"seed_ts": "", "seed_neural": ""}

        seed_ts = os.path.join(subj_out, f"Seed_{seed}_timeSeries.1D")
        seed_neural = os.path.join(subj_out, f"Seed_{seed}_neural.1D")
        if not os.path.exists(seed_ts) or not os.path.exists(seed_neural):

            # make seed, get timeseries
            print(f"\nMaking seed timeseries for {subj} - {seed}")
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
                    {clean_data} > {seed_ts}

                3dTfitter \
                    -RHS {seed_ts} \
                    -FALTUNG {hrf_file} \
                    tmp.1D 012 0

                1dtranspose \
                    tmp.1D > {seed_neural} && \
                    rm tmp.1D
            """
            subj_num = subj.split("-")[1]
            h_out, h_err = submit_hpc_sbatch(
                h_cmd, 1, 4, 1, f"{subj_num}mksd", subj_out
            )

        # update seed_dict with upsampled seed
        seed_dict[seed]["seed_ts"] = seed_ts
        seed_dict[seed]["seed_neural"] = seed_neural

    # update afni_data with all seed_files
    afni_data["seed_files"] = seed_dict
    return afni_data


def behavior_timeseries(subj, sess, task, subj_out, afni_data, stim_dur=2):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    # get afni_data values
    len_tr = afni_data["len_tr"]
    seed_dict = afni_data["seed_files"]
    tf_list = afni_data["timing_files"]

    # list of run length in seconds, total number of volumes (sum_vol)
    run_len = []
    num_vol = []
    for scale_file in afni_data["scaled_files"]:
        h_cmd = f"3dinfo -ntimes {scale_file}"
        h_out, h_err = submit_hpc_subprocess(h_cmd)
        h_vol = int(h_out.decode("utf-8").strip())
        num_vol.append(h_vol)
        run_len.append(str(h_vol * len_tr))
    sum_vol = sum(num_vol)

    # start dict of final ts files
    final_dict = {}
    for seed_name in seed_dict:
        final_dict[seed_name] = {}

    # get seed ts for e/behavior
    for timing_file in tf_list:

        # get binary behavior file, in TR time
        h_beh = timing_file.split("desc-")[1].split("_")[0]
        beh_file = os.path.join(subj_out, f"Beh_{h_beh}_bin.1D")
        if not os.path.exists(beh_file):
            print(f"\nMaking behavior files for {h_beh}")
            h_cmd = f"""
                timing_tool.py \
                    -timing {timing_file} \
                    -tr {len_tr} \
                    -stim_dur {stim_dur} \
                    -run_len {" ".join(run_len)} \
                    -min_frac 0.3 \
                    -timing_to_1D {beh_file}
            """
            h_out, h_err = submit_hpc_subprocess(h_cmd)
        assert os.path.exists(beh_file), f"Failed to properly make {beh_file}"

        # multiply beh_file by neural seed ts to get interaction term
        for seed_name in seed_dict:
            final_file = os.path.join(
                subj_out, f"Final_{seed_name}_{h_beh}_timeSeries.1D"
            )
            if not os.path.exists(final_file):
                print(f"\nMaking {final_file}")
                seed_neural = seed_dict[seed_name]["seed_neural"]
                h_cmd = f"""
                    cd {subj_out}

                    1deval \
                        -a {seed_neural} \
                        -b {beh_file} \
                        -expr 'a*b' \
                        > Seed_{seed_name}_{h_beh}_neural.1D

                    waver \
                        -FILE {len_tr} HRF_model.1D \
                        -peak 1 \
                        -TR {len_tr} \
                        -input Seed_{seed_name}_{h_beh}_neural.1D \
                        -numout {sum_vol} \
                        > {final_file}
                """
                subj_num = subj.split("-")[1]
                h_out, h_err = submit_hpc_sbatch(
                    h_cmd, 1, 4, 1, f"{subj_num}final", subj_out
                )
            assert os.path.exists(final_file), f"Failed to make {final_file}"
            final_dict[seed_name][f"{seed_name}_{h_beh}"] = final_file

    afni_data["final_files"] = final_dict
    return afni_data


def mot_files(subj_out, afni_data):
    """Constuct motion and censor files

    Mine <fMRIprep>_desc-confounds_timeseries.tsv for motion events, make
    motion files for mean (6df) and derivative (6df) motion events. Also,
    create motion censor file. Volume preceding a motion event is also
    censored. Finally, report the number of censored volumes.

    I'm not sure if motion is demeaned or not, given that
    it is output by fMRIprep (mined from confounds.tsv file).

    Parameters
    ----------
    subj_out : str
        /path/to/project_dir/derivatives/afni/sub-1234/ses-A
    afni_data : dict
        contains names for various files

    Returns
    -------
    afni_data : dict
        updated with names of motion files
        mot-mean = motion mean file
        mot-deriv = motion derivative file
        mot-censor = binary censory vector

    Notes
    -----
    Requires afni_data["conf_files"].

    As runs do not have an equal number of volumes, motion/censor files
    for each run are concatenated into a single file rather than managing
    zero padding.

    Writes 1D files - AFNI reads tsv as containing a header!
    """

    # determine relevant col labels
    mean_labels = [
        "trans_x",
        "trans_y",
        "trans_z",
        "rot_x",
        "rot_y",
        "rot_z",
    ]

    drv_labels = [
        "trans_x_derivative1",
        "trans_y_derivative1",
        "trans_z_derivative1",
        "rot_x_derivative1",
        "rot_y_derivative1",
        "rot_z_derivative1",
    ]

    # get list of confound timeseries files, set output strings
    mot_list = afni_data["conf_files"]
    mot_str = os.path.join(subj_out, mot_list[0].split("/")[-1].replace("run-1_", ""))
    mean_str = f"""{mot_str.replace("confounds", "mean").replace(".tsv", ".1D")}"""
    deriv_str = f"""{mot_str.replace("confounds", "deriv").replace(".tsv", ".1D")}"""
    cens_str = f"""{mot_str.replace("confounds", "censor").replace(".tsv", ".1D")}"""

    if not os.path.exists(cens_str):
        print("Making motion mean, derivative, censor files ...")

        # start empty lists to append
        mean_cat = []
        deriv_cat = []
        censor_cat = []

        # each confound/motion file
        for mot_file in mot_list:

            # read in data
            df_all = pd.read_csv(mot_file, sep="\t")

            # make motion mean file, round to 6 sig figs
            df_mean = df_all[mean_labels].copy()
            df_mean = df_mean.round(6)
            mean_cat.append(df_mean)

            # make motion deriv file
            df_drv = df_all[drv_labels].copy()
            df_drv = df_drv.fillna(0)
            df_drv = df_drv.round(6)
            deriv_cat.append(df_drv)

            # make motion censor file - sum columns,
            # invert binary, exclude preceding volume
            df_cen = df_all.filter(regex="motion_outlier")
            df_cen["sum"] = df_cen.iloc[:, :].sum(1)
            df_cen = df_cen.astype(int)
            df_cen = df_cen.replace({0: 1, 1: 0})
            zero_pos = df_cen.index[df_cen["sum"] == 0].tolist()
            zero_fill = [x - 1 for x in zero_pos]
            if -1 in zero_fill:
                zero_fill.remove(-1)
            df_cen.loc[zero_fill, "sum"] = 0
            censor_cat.append(df_cen)

        # cat files into singule file rather than pad zeros for e/run
        df_mean_cat = pd.concat(mean_cat, ignore_index=True)
        df_deriv_cat = pd.concat(deriv_cat, ignore_index=True)
        df_censor_cat = pd.concat(censor_cat, ignore_index=True)

        # determine BIDS string, write tsvs, make sure
        # output value is float (not scientific notation)
        df_mean_cat.to_csv(
            mean_str,
            sep="\t",
            index=False,
            header=False,
            float_format="%.6f",
        )
        df_deriv_cat.to_csv(
            deriv_str,
            sep="\t",
            index=False,
            header=False,
            float_format="%.6f",
        )
        df_censor_cat.to_csv(
            cens_str,
            sep="\t",
            index=False,
            header=False,
            columns=["sum"],
        )

    # update afni_data
    afni_data["mot-mean"] = mean_str
    afni_data["mot-deriv"] = deriv_str
    afni_data["mot-censor"] = cens_str

    return afni_data


def write_ppi_decon(subj, decon_ppi, subj_out, seed, afni_data, dur=2):
    """Generate deconvolution script.

    Write a deconvolution script using the pre-processed data, motion, and
    censored files passed by afni_data. Uses a 2GAM basis function
    (AFNI's TWOGAMpw). This script is used to generate X.files and the
    foo_stats.REML_cmd.

    Timing files should contain AFNI-formatted onset times (duration is hardcoded),
    using the asterisk for runs in which a certain behavior does not occur.

    Parameters
    ----------



    decon_name: str
        name of deconvolution, useful when conducting multiple
        deconvolutions on same session. Will be appended to
        BIDS task name (decon_<task-name>_<decon_name>).
    afni_data : dict
        contains names for various files
    work_dir : str
        /path/to/project_dir/derivatives/afni/sub-1234/ses-A
    dur : int/float/str
        duration of event to model

    Returns
    -------
    afni_data : dict
        updated with REML commands
        {"dcn-<decon_name>": foo_stats.REML_cmd}

    Notes
    -----
    Requires afni_data["epi-scale*"], afni_data["mot-mean"],
        afni_data["mot-deriv"], and afni_data["mot-censor"].
    Deconvolution files will be written in AFNI format, rather
        than BIDS. This includes the X.files (cue spooky theme), script,
        and deconvolved output. Files names will have the format:
            decon_<bids-task>_<decon_name>
    """

    print("\nBuilding decon script")

    # make timing file dict
    tf_dict = {}
    for tf in afni_data["timing_files"]:
        beh = tf.split("/")[-1].split("desc-")[1].split("_")[0]
        tf_dict[beh] = tf

    # set regressors - baseline
    reg_base = [
        f"""-ortvec {afni_data["mot-mean"]} mot_mean""",
        f"""-ortvec {afni_data["mot-deriv"]} mot_deriv""",
    ]

    # set regressors - behavior
    reg_beh = []
    for c_beh, beh in enumerate(tf_dict):

        # add stim_time info, order is
        #   -stim_times 1 tf_beh.txt basisFunction
        reg_beh.append("-stim_times")
        reg_beh.append(f"{c_beh + 1}")
        reg_beh.append(f"{tf_dict[beh]}")
        reg_beh.append(f"'TWOGAMpw(4,5,0.2,12,7,{dur})'")

        # add stim_label info, order is
        #   -stim_label 1 beh
        reg_beh.append("-stim_label")
        reg_beh.append(f"{c_beh + 1}")
        reg_beh.append(beh)

    # increase behavior counter for 0-indexing, use with ppi files
    c_beh += 2

    # set regressors - seed timeseries
    seed_ts = afni_data["seed_files"][seed]["seed_ts"]
    reg_beh.append(f"-stim_file {c_beh} {seed_ts}")
    reg_beh.append(f"-stim_label {c_beh} Seed_ts")

    # set regressors - seed behavior timeseries
    final_list = afni_data["final_files"][seed]
    for final_name, final_file in final_list.items():
        c_beh += 1
        reg_beh.append(f"-stim_file {c_beh} {final_file}")
        reg_beh.append(f"-stim_label {c_beh} S{final_name.split('_')[1]}")

    # input files
    epi_list = afni_data["scaled_files"]

    # build full decon command
    cmd_decon = f"""
        3dDeconvolve \\
            -x1D_stop \\
            -GOFORIT \\
            -input {" ".join(epi_list)} \\
            -censor {afni_data["mot-censor"]} \\
            {" ".join(reg_base)} \\
            -polort A \\
            -float \\
            -local_times \\
            -num_stimts {c_beh} \\
            {" ".join(reg_beh)} \\
            -jobs 1 \\
            -x1D {subj_out}/X.{decon_ppi}.xmat.1D \\
            -xjpeg {subj_out}/X.{decon_ppi}.jpg \\
            -x1D_uncensored {subj_out}/X.{decon_ppi}.nocensor.xmat.1D \\
            -bucket {subj_out}/{decon_ppi}_stats \\
            -errts {subj_out}/{decon_ppi}_errts
    """

    # write for review
    decon_script = os.path.join(subj_out, f"{decon_ppi}.sh")
    with open(decon_script, "w") as script:
        script.write(cmd_decon)

    # generate x-matrices, reml command
    out_file = os.path.join(subj_out, f"{decon_ppi}_stats.REML_cmd")
    if not os.path.exists(out_file):
        print(f"Running 3dDeconvolve for {decon_ppi}")
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(cmd_decon, 1, 4, 1, f"{subj_num}dcn", subj_out)

    # check, update afni_data
    assert os.path.exists(out_file), f"Failed to write {out_file}"
    afni_data[decon_ppi] = out_file
    return afni_data


def run_ppi_reml(subj, subj_out, decon_ppi, afni_data):
    """Deconvolve EPI timeseries.

    Generate an idea of nuissance signal from the white matter and
    include this in the generated 3dREMLfit command.

    Parameters
    ----------
    work_dir : str
        /path/to/project_dir/derivatives/afni/sub-1234/ses-A
    afni_data : dict
        contains names for various files

    Returns
    -------
    afni_data : dict
        updated for nuissance, deconvolved files
        epi-nuiss = nuissance signal file
        rml-<decon_name> = deconvolved file (<decon_name>_stats_REML+tlrc)
    """

    subj_num = subj.split("-")[-1]
    epi_list = afni_data["scaled_files"]
    eroded_mask = afni_data["wme_mask"]
    h_nuiss = epi_list[0].split("/")[-1].replace("_run-1", "")
    nuiss_file = os.path.join(
        subj_out, h_nuiss.replace("desc-scaled", "desc-nuissance")
    )

    # generate WM timeseries (nuissance file) for task
    if not os.path.exists(nuiss_file):
        print(f"Making nuissance file {nuiss_file} ...")
        tcat_file = "tmp_tcat.sub".join(nuiss_file.rsplit("sub", 1))
        epi_eroded = "tmp_epi.sub".join(eroded_mask.rsplit("sub", 1))
        h_cmd = f"""
            3dTcat -prefix {tcat_file} {" ".join(epi_list)}

            3dcalc \
                -a {tcat_file} \
                -b {eroded_mask} \
                -expr 'a*bool(b)' \
                -datum float \
                -prefix {epi_eroded}

            3dmerge \
                -1blur_fwhm 20 \
                -doall \
                -prefix {nuiss_file} \
                {epi_eroded}
        """
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}wts", subj_out)
    assert os.path.exists(nuiss_file), f"Failed to write {nuiss_file}"
    afni_data["epi-nuiss"] = nuiss_file

    # run each planned deconvolution
    reml_script = afni_data[decon_ppi]
    reml_out = reml_script.replace(".REML_cmd", "_REML+tlrc.HEAD")
    if not os.path.exists(reml_out):
        print(f"Starting 3dREMLfit for {decon_ppi} ...")
        h_cmd = f"""
            tcsh \
                -x {reml_script} \
                -dsort {afni_data["epi-nuiss"]} \
                -GOFORIT
        """
        h_out, h_err = submit_hpc_sbatch(h_cmd, 25, 4, 6, f"{subj_num}rml", subj_out)
    assert os.path.exists(reml_out), f"Failed to write {reml_out}"
    afni_data[f"rml-{decon_ppi}"] = reml_out.split(".")[0]
    return afni_data


def copy_data(subj_out, subj_data, decon_ppi):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """
    h_cmd = f"cp {subj_out}/{{,X.}}{decon_ppi}.* {subj_data}"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    check_file = os.path.join(subj_data, f"{decon_ppi}_stats_REML+tlrc.HEAD")
    assert os.path.exists(check_file), f"Missing PPI decon file in {subj_data}"


def get_args():
    """Get and parse arguments"""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--data-dir",
        type=str,
        default="/home/data/madlab/McMakin_EMUR01/derivatives/afni",
        help=textwrap.dedent(
            """\
            Path to BIDS-formatted derivatives directory containing output
            of github.com/emu-project/func_processing/cli/run_afni.py
            (default : %(default)s)
            """
        ),
    )
    parser.add_argument(
        "--deriv-dir",
        type=str,
        default="/scratch/madlab/emu_unc/derivatives/afni_ppi",
        help=textwrap.dedent(
            """\
            Path to desired scratch location, for intermediates.
            (default : %(default)s)
            """
        ),
    )
    parser.add_argument(
        "--sess",
        type=str,
        default="ses-S2",
        help=textwrap.dedent(
            """\
            BIDS-formatted session string
            (default : %(default)s)
            """
        ),
    )
    parser.add_argument(
        "--task",
        type=str,
        default="task-test",
        help=textwrap.dedent(
            """\
            BIDS-formatted task string
            (default : %(default)s)
            """
        ),
    )

    required_args = parser.add_argument_group("Required Arguments")
    required_args.add_argument(
        "-s",
        "--subj",
        help="BIDS subject str (sub-1234)",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-d",
        "--decon-name",
        help="Prefix of decon file (decon_task-test_UniqueBehs)",
        type=str,
        required=True,
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser


# %%
def main():

    # # For testing
    # data_dir = "/home/data/madlab/McMakin_EMUR01/derivatives/afni"
    # deriv_dir = "/scratch/madlab/emu_unc/derivatives/afni_ppi"
    # subj = "sub-4001"
    # sess = "ses-S2"
    # task = "task-test"
    # decon_str = f"decon_{task}_UniqueBehs"

    # check for correct conda env
    assert (
        "emuR01_unc_env" in sys.executable
    ), "Please activate emuR01_unc conda environment."

    # get passed args
    args = get_args().parse_args()
    data_dir = args.data_dir
    deriv_dir = args.deriv_dir
    subj = args.subj
    sess = args.sess
    task = args.task
    decon_str = args.decon_name

    # make dict for seed construction from coordinates
    seed_coord = {"LHC": "-24 -12 -22"}

    # setup paths
    subj_data = os.path.join(data_dir, subj, sess, "func")
    subj_out = os.path.join(deriv_dir, subj, sess, "func")
    if not os.path.exists(subj_out):
        os.makedirs(subj_out)

    # get required files produced by
    # github.com/emu-project/func_processing/cli/run_afni.py
    afni_data = {}
    afni_data["scaled_files"] = sorted(
        glob.glob(f"{subj_data}/*desc-scaled_bold.nii.gz")
    )
    afni_data["decon_file"] = f"{subj_data}/{decon_str}_cbucket_REML+tlrc.HEAD"
    afni_data["decon_matrix"] = f"{subj_data}/X.{decon_str}.xmat.1D"
    afni_data["timing_files"] = glob.glob(
        f"{subj_data}/{subj}_{sess}_{task}_desc-*_events.1D"
    )
    afni_data["conf_files"] = sorted(
        glob.glob(f"{subj_data}/*_run-*_desc-confounds_timeseries.tsv")
    )
    subj_anat = os.path.join(data_dir, subj, sess, "anat")
    afni_data["wme_mask"] = glob.glob(f"{subj_anat}/*desc-WMe_mask.nii.gz")

    # check that required files exist
    for h_type, h_file in afni_data.items():
        if type(h_file) == list:
            assert os.path.exists(h_file[0]), f"Missing files for {h_type}"
        else:
            assert os.path.exists(h_file), f"Missing file for {h_type}: {h_file}"

    # generate timeseries files
    afni_data = clean_data(subj, subj_out, afni_data)
    afni_data = hrf_model(subj_out, afni_data)
    afni_data = seed_timeseries(subj, subj_out, afni_data, seed_coord)
    afni_data = behavior_timeseries(subj, sess, task, subj_out, afni_data)

    # do decons for each seed
    afni_data = mot_files(subj_out, afni_data)
    for seed in seed_coord:
        decon_ppi = f"{decon_str}_PPI-{seed}"
        write_ppi_decon(subj, decon_ppi, subj_out, seed, afni_data)
        run_ppi_reml(subj, subj_out, decon_ppi, afni_data)
        copy_data(subj_out, subj_data, decon_ppi)

    # clean up
    shutil.rmtree(os.path.join(deriv_dir, subj))


# %%
if __name__ == "__main__":
    main()
