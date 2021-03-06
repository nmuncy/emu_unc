#!/usr/bin/env python

r"""Conduct PPI Analysis.

Find averaged timeseries of seed, solve for RHS to get
neural timeseries, make behavior vectors, then convolve
to get behavior timeseries of seed.

Supports coordinates for seed construction, or passed
ROI masks.

Examples
--------
sbatch --job-name=p1234 \
    --output=p1234 \
    --mem-per-cpu=4000 \
    --partition=IB_44C_512G \
    --account=iacc_madlab \
    --qos=pq_madlab \
    func4_ppi.py \
    -s sub-1234 \
    -d decon_task-study_precTest \
    -r LHC \
    -i "-24 -12 -22"

sbatch --job-name=p1234 \
    --output=p1234 \
    --mem-per-cpu=4000 \
    --partition=IB_44C_512G \
    --account=iacc_madlab \
    --qos=pq_madlab \
    func4_ppi.py \
    -s sub-1234 \
    -d decon_task-test_UniqueBehs \
    -r amgL \
    -i "/path/to/amgL_mask.ni.gz"
"""

# %%
import os
import sys
import glob
import fnmatch
import shutil
import subprocess
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


def clean_data(subj, subj_work, afni_data):
    """Construct Clean Data.

    Concatenate pre-processed EPI files and remove
    "effects of no interest" i.e. data from deconvolved
    sub-bricks beginning with "Run" or "mot" (these are
    the polynomial and motion baseline regressors).

    Writes <subj_work>/CleanData+tlrc

    Parameters
    ----------
    subj : str
        BIDS-formatted subject string
    subj_work : str
        Path to subject scratch directory
    afni_data : dict
        contains the fields decon_file, decon_matrix, and
        scaled_files

    Returns
    -------
    afni_data : dict
        added fields of len_tr, clean_data

    Raises
    ------
    AssertionError
        missing <subj_work>/CleanData+tlrc.HEAD
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
    clean_file = os.path.join(subj_work, "CleanData+tlrc.HEAD")
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
            cd {subj_work}

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
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}cle", subj_work)

    assert os.path.exists(clean_file), "CleanData not found."
    afni_data["len_tr"] = len_tr
    afni_data["clean_data"] = clean_file
    return afni_data


def hrf_model(subj_work, afni_data, dur=2):
    """Create HRF model.

    Generate an idealized HRF function using the same
    basis function used in the deconvolution (2GAM).

    Writes <subj_work>/HRF_model.1D

    Parameters
    ----------
    subj_work : str
        Path to subject scratch directory
    afni_data : dict
        contains the fields len_tr
    dur : int/float
        duration of event to model
        (default : 2)

    Returns
    -------
    afni_data : dict
        added field of hrf_model

    Raises
    ------
    AssertionError
        missing <subj_work>/HRF_model.1D
    """
    # make ideal HRF, use same model as deconvolution
    len_tr = afni_data["len_tr"]
    hrf_file = os.path.join(subj_work, "HRF_model.1D")
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


def seed_timeseries(subj, subj_work, afni_data, seed_tuple):
    """Get the seed timeseries.

    Extract an averaged timeseries for the voxels labeled
    by a seed. Will construct a seed with srad=3 if coordinates
    are supplied. The modeled HRF is then deconvolved from the
    timeseries.

    Writes <subj_work>/Seed_<seed>_<timeSeries|neural>.1D

    Parameters
    ----------
    subj : str
        BIDS-formatted subject string
    subj_work : str
        Path to subject scratch directory
    afni_data : dict
        contains the fields hrf_model, clean_data
    seed_tuple : tuple
        [0] = seed name (str)
        [1] = seed info (str), coordinates or path to mask

    Returns
    -------
    afni_data : dict
        added field seed_files

    Raises
    ------
    AssertionError
        missing <subj_work>/Seed_<seed>+tlrc.HEAD
    """
    hrf_file = afni_data["hrf_model"]
    clean_data = afni_data["clean_data"]
    seed = seed_tuple[0]
    seed_info = seed_tuple[1]

    seed_out = {}
    seed_out[seed] = {"seed_ts": "", "seed_neural": ""}
    seed_ts = os.path.join(subj_work, f"Seed_{seed}_timeSeries.1D")
    seed_neural = os.path.join(subj_work, f"Seed_{seed}_neural.1D")
    if not os.path.exists(seed_ts) or not os.path.exists(seed_neural):

        # make seed if coordinates supplied as seed_info
        if len(seed_info.split(" ")) != 1:
            seed_file = os.path.join(subj_work, f"Seed_{seed}+tlrc.HEAD")
            if not os.path.exists(seed_file):
                h_cmd = f"""
                    echo {seed_info} | 3dUndump \
                        -xyz \
                        -srad 3 \
                        -master {clean_data} \
                        -prefix {subj_work}/Seed_{seed} -
                """
                subj_num = subj.split("-")[1]
                h_out, h_err = submit_hpc_sbatch(
                    h_cmd, 1, 1, 1, f"{subj_num}mksd", subj_work
                )
        else:
            seed_file = seed_info
        assert os.path.exists(seed_file), "Seed file not found."

        # get timeseries
        print(f"\nMaking seed timeseries for {subj} - {seed}")
        h_cmd = f"""
            cd {subj_work}

            3dmaskave \
                -quiet \
                -mask {seed_file} \
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
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}mksd", subj_work)

    # update seed_out with upsampled seed
    seed_out[seed]["seed_ts"] = seed_ts
    seed_out[seed]["seed_neural"] = seed_neural

    # update afni_data with all seed_files
    afni_data["seed_files"] = seed_out
    return afni_data


def behavior_timeseries(subj, sess, subj_work, afni_data, stim_dur=2):
    """Extract behavior timeseries.

    Convert timing files to binary vectors, with one binary value for each
    EPI volume. Multiply by the seed neural timeseries to get seed-behavior
    neural timeseries, and then convolve with HRF model to get seed-behavior
    HRF timeseries.

    Writes <subj_work>/Final_<seed>_<behavior>_timeSeries.1D

    Parameters
    ----------
    subj : str
        BIDS-formatted subject string
    sess : str
        BIDS-formatted session string
    subj_work : str
        Path to subject scratch directory
    afni_data : dict
        contains the fields len_tr, seed_files, timing_files,
        and scaled_files
    stim_dur : int/float
        duration of stimulus to model
        (default : 2)

    Returns
    -------
    afni_data : dict
        added field final_files

    Raises
    ------
    AssertionError
        missing <subj_work>/Beh_<behavior>_bin.1D
        missing <subj_work>/Final_<seed>_<behavior>_timeSeries.1D
    """
    # get afni_data values
    len_tr = afni_data["len_tr"]
    seed_out = afni_data["seed_files"]
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
    for seed_name in seed_out:
        final_dict[seed_name] = {}

    # get seed ts for e/behavior
    for timing_file in tf_list:

        # get binary behavior file, in TR time
        h_beh = timing_file.split("desc-")[1].split("_")[0]
        beh_file = os.path.join(subj_work, f"Beh_{h_beh}_bin.1D")
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
        for seed_name in seed_out:
            final_file = os.path.join(
                subj_work, f"Final_{seed_name}_{h_beh}_timeSeries.1D"
            )
            if not os.path.exists(final_file):
                print(f"\nMaking {final_file}")
                seed_neural = seed_out[seed_name]["seed_neural"]
                h_cmd = f"""
                    cd {subj_work}

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
                    h_cmd, 1, 4, 1, f"{subj_num}final", subj_work
                )
            assert os.path.exists(final_file), f"Failed to make {final_file}"
            final_dict[seed_name][f"{seed_name}_{h_beh}"] = final_file

    afni_data["final_files"] = final_dict
    return afni_data


def write_ppi_decon(subj, decon_ppi, subj_work, seed, afni_data, dur=2):
    """Generate deconvolution script.

    Write a deconvolution script using the pre-processed data, motion, PPI, and
    censored files passed by afni_data. Uses a 2GAM basis function
    (AFNI's TWOGAMpw).

    Timing files should contain AFNI-formatted onset times (duration is hardcoded),
    using the asterisk for runs in which a certain behavior does not occur.

    Writes X.<decon_ppi> files, <decon_ppi>.sh for review, and <decon_ppi>.REML_cmd.

    Parameters
    ----------
    subj : str
        BIDS-formatted subject string
    decon_ppi : str
        prefix of output PPI deconvolved file,
        e.g. decon_<task-name>_<deconName>_PPI-<seed>
    subj_work : str
        Path to subject scratch directory
    seed : str
        name of seed
    afni_data : dict
        contains the fields timing_files, mot-mean, mot-deriv, mot-censor,
        seed_files, final_files, scaled_files
    dur : int/float
        duration of behavior to model
        (default : 2)

    Returns
    -------
    afni_data : dict
        updated with REML commands
        {decon_ppi: foo_stats.REML_cmd}

    Raises
    ------
    AssertionError
        missing <subj_work>/<decon_ppi>_stats.REML_cmd
    """
    print("\nBuilding decon script")

    # make timing file dict
    tf_dict = {}
    for tf in afni_data["timing_waves"]:
        beh = tf.split("/")[-1].split("desc-")[1].split("Cens_")[0]
        tf_dict[beh] = tf

    # set regressors - baseline
    reg_base = [
        f"""-ortvec {afni_data["mot-mean"]} mot_mean""",
        f"""-ortvec {afni_data["mot-deriv"]} mot_deriv""",
    ]

    # set inverted censor as baseline regressor, start regressor count
    c_beh = 1
    reg_cens = [f"-stim_file {c_beh} {afni_data['mot-censorInv']}"]
    reg_cens.append(f"-stim_base {c_beh}")
    reg_cens.append(f"-stim_label {c_beh} mot_cens")

    # set behavior regressors
    reg_beh = []
    for h_beh, h_tf in tf_dict.items():
        c_beh += 1
        reg_beh.append(f"-stim_file {c_beh} {h_tf}")
        reg_beh.append(f"-stim_label {c_beh} {h_beh}")

    # set regressors - seed timeseries
    c_beh += 1
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
            {" ".join(reg_base)} \\
            -polort A \\
            -float \\
            -local_times \\
            -num_stimts {c_beh} \\
            {" ".join(reg_cens)} \\
            {" ".join(reg_beh)} \\
            -jobs 1 \\
            -x1D {subj_work}/X.{decon_ppi}.xmat.1D \\
            -xjpeg {subj_work}/X.{decon_ppi}.jpg \\
            -x1D_uncensored {subj_work}/X.{decon_ppi}.nocensor.xmat.1D \\
            -bucket {subj_work}/{decon_ppi}_stats \\
            -errts {subj_work}/{decon_ppi}_errts
    """

    # write for review
    decon_script = os.path.join(subj_work, f"{decon_ppi}.sh")
    with open(decon_script, "w") as script:
        script.write(cmd_decon)

    # generate x-matrices, reml command
    out_file = os.path.join(subj_work, f"{decon_ppi}_stats.REML_cmd")
    if not os.path.exists(out_file):
        print(f"Running 3dDeconvolve for {decon_ppi}")
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(
            cmd_decon, 1, 4, 1, f"{subj_num}dcn", subj_work
        )

    # check, update afni_data
    assert os.path.exists(out_file), f"Failed to write {out_file}"
    afni_data[decon_ppi] = out_file
    return afni_data


def run_ppi_reml(subj, subj_work, decon_ppi, afni_data):
    """Deconvolve EPI timeseries.

    Generate an idea of nuissance signal from the white matter and
    include this in the generated 3dREMLfit command.

    Writes <subj_work>/<decon_ppi>_stats_REML+tlrc

    Parameters
    ----------
    subj : str
        BIDS-formatted subject string
    subj_work : str
        Path to subject scratch directory
    decon_ppi : str
        prefix of decon PPI file
    afni_data : dict
        contains fields scaled_files, wme_mask

    Returns
    -------
    afni_data : dict
        updated for nuissance, deconvolved files
        epi-nuiss = nuissance signal file
        rml-<decon_ppi> = deconvolved file (<decon_ppi>_stats_REML+tlrc)

    Raises
    ------
    AssertionError
        missing <subj_work>/nuissance file
        missing <subj_work>/<decon_ppi>_stats_REML+tlrc.HEAD
    """
    subj_num = subj.split("-")[-1]
    epi_list = afni_data["scaled_files"]
    eroded_mask = afni_data["wme_mask"][0]
    h_nuiss = epi_list[0].split("/")[-1].replace("_run-1", "")
    nuiss_file = os.path.join(
        subj_work, h_nuiss.replace("desc-scaled", "desc-nuissance")
    )

    # generate WM timeseries (nuissance file) for task
    if not os.path.exists(nuiss_file):
        print(f"Making nuissance file {nuiss_file} ...")
        tcat_file = "tmp_tcat.sub".join(nuiss_file.rsplit("sub", 1))
        h_eroded = eroded_mask.split("/")[-1]
        epi_eroded = os.path.join(subj_work, f"tmp_epi.{h_eroded}")
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
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}wts", subj_work)
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
        h_out, h_err = submit_hpc_sbatch(h_cmd, 25, 4, 8, f"{subj_num}rml", subj_work)
    assert os.path.exists(reml_out), f"Failed to write {reml_out}"
    afni_data[f"rml-{decon_ppi}"] = reml_out.split(".")[0]
    return afni_data


def copy_data(subj_work, subj_func, decon_ppi):
    """Copy final files to storage.

    Move only the final files to storage, those with
    <decon_ppi> in the file name.

    Parameters
    ----------
    subj_work : str
        Path to subject scratch directory
    subj_func : str
        Path to subect's project derivative directory
    decon_ppi : str
        prefix of deconvolved file

    Raises
    ------
    AssertionError
        missing <subj_func>/<decon_ppi>_stats_REML+tlrc.HEAD
    """
    h_cmd = f"cp {subj_work}/{{,X.}}{decon_ppi}* {subj_func}"
    h_out, h_err = submit_hpc_subprocess(h_cmd)
    check_file = os.path.join(subj_func, f"{decon_ppi}_stats_REML+tlrc.HEAD")
    assert os.path.exists(check_file), f"Missing PPI decon file in {subj_func}"


def get_args():
    """Get and parse arguments."""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--deriv-dir",
        type=str,
        default="/home/data/madlab/McMakin_EMUR01/derivatives/emu_unc",
        help=textwrap.dedent(
            """\
            Path to BIDS-formatted derivatives directory containing output
            of github.com/emu-project/func_processing/cli/run_afni.py
            (default : %(default)s)
            """
        ),
    )
    parser.add_argument(
        "--work-dir",
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
        "--timing-dir",
        type=str,
        default=None,
        help=textwrap.dedent(
            """\
            Path to subject timing files, if they are not
            located in derivatives/foo/subj-1234/ses-1/func.
            (default : %(default)s)
            """
        ),
    )

    required_args = parser.add_argument_group("Required Arguments")
    required_args.add_argument(
        "-p",
        "--subj",
        help="BIDS subject str (sub-1234)",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-s", "--sess", type=str, help="BIDS-formatted session string", required=True
    )
    required_args.add_argument(
        "-t", "--task", type=str, help="BIDS-formatted task string", required=True
    )
    required_args.add_argument(
        "-d",
        "--decon-name",
        help="Prefix of decon file (decon_task-test_UniqueBehs)",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-r",
        "--seed-name",
        help="Name of seed supplied/to be made.",
        type=str,
        required=True,
    )
    required_args.add_argument(
        "-i",
        "--seed-info",
        help="Coordinates of or path to seed.",
        type=str,
        required=True,
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser


# %%
def main():
    """Coordinate work."""
    # get passed args
    args = get_args().parse_args()
    deriv_dir = args.deriv_dir
    work_dir = args.work_dir
    timing_dir = args.timing_dir
    subj = args.subj
    sess = args.sess
    task = args.task
    decon_str = args.decon_name
    seed_name = args.seed_name
    seed_info = args.seed_info

    # setup paths/dicts
    subj_func = os.path.join(deriv_dir, subj, sess, "func")
    subj_anat = os.path.join(deriv_dir, subj, sess, "anat")
    subj_work = os.path.join(work_dir, subj, sess, "func")
    subj_time = timing_dir if timing_dir else subj_func
    print(subj_time)
    if not os.path.exists(subj_work):
        os.makedirs(subj_work)

    # Get required files produced by
    # github.com/emu-project/func_processing/cli/run_afni_subj.py
    # using deconvolve.write_new_decon function.
    afni_data = {}
    h_dcn = decon_str.split("_")[-1]
    timing_all = sorted(
        glob.glob(f"{subj_time}/*_{task}_decon-{h_dcn}_desc-*_events.*")
    )
    afni_data["timing_files"] = [
        x for x in timing_all if not fnmatch.fnmatch(x, "*Cens_events.txt")
    ]
    afni_data["timing_waves"] = sorted(
        glob.glob(f"{subj_time}/*_{task}_decon-{h_dcn}_desc-*Cens_events.*")
    )
    afni_data["scaled_files"] = sorted(
        glob.glob(f"{subj_func}/*{sess}_{task}*desc-scaled_bold.nii.gz")
    )
    afni_data["decon_file"] = f"{subj_func}/{decon_str}_cbucket_REML+tlrc.HEAD"
    afni_data["decon_matrix"] = f"{subj_func}/X.{decon_str}.xmat.1D"
    afni_data["mot-mean"] = f"{subj_func}/{subj}_{sess}_{task}_desc-mean_timeseries.1D"
    afni_data[
        "mot-deriv"
    ] = f"{subj_func}/{subj}_{sess}_{task}_desc-deriv_timeseries.1D"
    afni_data[
        "mot-censorInv"
    ] = f"{subj_func}/{subj}_{sess}_{task}_desc-censorInv_timeseries.1D"
    afni_data["wme_mask"] = glob.glob(f"{subj_anat}/*desc-WMe_mask.nii.gz")
    print(afni_data)

    # check that required files exist
    for h_type, h_file in afni_data.items():
        if type(h_file) == list:
            assert os.path.exists(h_file[0]), f"Missing files for {h_type}"
        else:
            assert os.path.exists(h_file), f"Missing file for {h_type}: {h_file}"

    # generate timeseries files
    seed_tuple = (seed_name, seed_info)
    afni_data = clean_data(subj, subj_work, afni_data)
    afni_data = hrf_model(subj_work, afni_data)
    afni_data = seed_timeseries(subj, subj_work, afni_data, seed_tuple)
    afni_data = behavior_timeseries(subj, sess, subj_work, afni_data)

    # do decons for each seed
    decon_ppi = f"{decon_str}_PPI-{seed_name}"
    afni_data = write_ppi_decon(subj, decon_ppi, subj_work, seed_name, afni_data)
    afni_data = run_ppi_reml(subj, subj_work, decon_ppi, afni_data)
    copy_data(subj_work, subj_func, decon_ppi)

    # clean up
    shutil.rmtree(os.path.join(work_dir, subj))


# %%
if __name__ == "__main__":
    env_found = [x for x in sys.path if "emuR01_unc" in x]
    if not env_found:
        print("\nERROR: madlab conda emuR01_unc_env required.")
        print("\tHint: $madlab_env emuR01_unc\n")
        sys.exit()
    main()
