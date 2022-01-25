#!/usr/bin/env python

"""Make kmeans bla masks.

Use bla amygdala masks from kmeans to make
resampled masks in template space.

Examples
--------
sbatch --job-name=p1234 \\
    --output=p1234 \\
    --mem-per-cpu=4000 \\
    --partition=IB_44C_512G \\
    --account=iacc_madlab \\
    --qos=pq_madlab \\
    func0_masks.py \\
    -s sub-1234
"""

# %%
import os
import sys
import json
import glob
import textwrap
import subprocess
from argparse import ArgumentParser, RawTextHelpFormatter


# %%
def submit_hpc_sbatch(command, wall_hours, mem_gig, num_proc, job_name, work_dir):
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
        -o {work_dir}/{job_name}.out \
        -e {work_dir}/{job_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wait \
        --wrap="module load afni-20.2.06
            module load c3d-1.0.0-gcc-8.2.0
            module load ants-2.3.5
            {command}
        "
    """
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_response.communicate()[0].decode("utf-8")
    return (job_name, job_id)


def ants_warp(subj, work_dir, file_dict):
    """Register DWI to T2w template.

    Calculate symmetric normalization warp vectors
    to move from native diffusion to template space.

    Parameters
    ----------
    subj : str
        BIDS subject string (sub-1234)
    work_dir : str
        path to scratch working directory for subj
    file_dict : dict
        key, path to needed files
        mni-atlas = path to MNI PediatricAsym cohort-5 T2w res-2
        ntv-b0 = path to subject b0 file

    Returns
    -------
        file_dict : dict
            updated with paths to new files
            mni-b0 = path to warped b0 file
            wrp-calc = path to normalizations calculations
    """
    ants_warp = os.path.join(work_dir, "ants_Warped.nii.gz")
    if not os.path.exists(ants_warp):
        h_cmd = f"""
            antsRegistrationSyN.sh \
                -d 3 \
                -f {file_dict["mni-atlas"]} \
                -m {file_dict["ntv-b0"]} \
                -o {work_dir}/ants_ \
                -n 4
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 2, 4, 4, f"{subj_num}ants", work_dir)

    assert os.path.exists(ants_warp), "ANTs Registration failed."
    file_dict["mni-b0"] = ants_warp
    file_dict["wrp-calc"] = os.path.join(work_dir, "ants_0GenericAffine.mat")
    return file_dict


# %%
def make_mask(subj, sess, work_dir, out_dir, file_dict):
    """Generate BLA masks in template space.

    Extract BLA mask, warp it into template space using
    calculations from ants_warp, resample into
    functional dimensions, and then binarize.

    Parameters
    ----------
    subj : str
        BIDS subject string (sub-1234)
    sess : str
        BIDS session string (ses-S2)
    work_dir : str
        path to scratch working directory for subj
    out_dir : str
        path to final output location for mask
    file_dict : dict
        key, path to needed files
        ntv-mask = kmeans mask in native space
        bla-num = BLA label value in kmeans mask
        mni-atlas = path to MNI atlas
        mni-b0 = path to warped b0 file
        wrp-calc = path to normalizations calculations
        epi-ref = path to reference EPI file

    Returns
    -------
    file_dict : dict
        updated with path to new file
        mni-mask = path to bla mask in template space
    """

    bla_num = file_dict["bla-num"]
    space = file_dict["mni-atlas"].split("/")[-1].split("tpl-")[1].split("_")[0]
    cohort = file_dict["mni-atlas"].split("/")[-1].split("cohort-")[1].split("_")[0]
    bla_mask = os.path.join(
        out_dir,
        f"{subj}_{sess}_space-{space}_cohort-{cohort}_res-2_desc-bla_mask.nii.gz",
    )
    if not os.path.exists(bla_mask):
        h_cmd = f"""
            c3d \
                {file_dict["ntv-mask"]} \
                -thresh {bla_num} {bla_num} 1 0 \
                -o {work_dir}/tmp_bla_thresh.nii.gz

            WarpImageMultiTransform \
                3 \
                {work_dir}/tmp_bla_thresh.nii.gz \
                {work_dir}/tmp_bla_warp.nii.gz \
                {file_dict["wrp-calc"]} \
                -R {file_dict["mni-b0"]}

            3dresample \
                -master {file_dict["epi-ref"]} \
                -rmode NN \
                -input {work_dir}/tmp_bla_warp.nii.gz \
                -prefix {work_dir}/tmp_bla_res.nii.gz

            c3d \
                {work_dir}/tmp_bla_res.nii.gz \
                -thresh 0.3 1 1 0 \
                -o {bla_mask}

            3drefit \
                -space MNI \
                {bla_mask}
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}wrp", work_dir)

    assert os.path.exists(bla_mask), f"Failed to make {bla_mask}"
    file_dict["mni-mask"] = bla_mask
    return file_dict


# %%
def get_args():
    """Get and parse arguments"""
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--proj-dir",
        type=str,
        default="/home/data/madlab/McMakin_EMUR01",
        help=textwrap.dedent(
            """\
            Path to BIDS-formatted derivatives directory containing output
            of github.com/emu-project/func_processing/cli/run_afni.py
            (default : %(default)s)
            """
        ),
    )
    parser.add_argument(
        "--scratch-dir",
        type=str,
        default="/scratch/madlab/emu_unc",
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

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser


# %%
def main():

    # # For testing
    # proj_dir = "/home/data/madlab/McMakin_EMUR01"
    # scratch_dir = "/scratch/madlab/emu_unc"
    # subj = "sub-4001"
    # sess = "ses-S2"
    # task = "task-test"

    # check for correct conda env
    assert (
        "emuR01_unc_env" in sys.executable
    ), "Please activate emuR01_unc conda environment."

    # get passed args
    args = get_args().parse_args()
    proj_dir = args.proj_dir
    scratch_dir = args.scratch_dir
    subj = args.subj
    sess = args.sess
    task = args.task

    # set up
    atlas_dir = "/home/data/madlab/atlases/templateflow/tpl-MNIPediatricAsym/cohort-5"
    kmean_dir = os.path.join(proj_dir, "derivatives/kmeans", subj, sess, "dwi")
    b0_dir = os.path.join(proj_dir, "derivatives/dwi_bedpostx", subj, "emubpx/b0_file")
    epi_dir = os.path.join(proj_dir, "derivatives/afni", subj, sess, "func")
    work_dir = os.path.join(scratch_dir, "derivatives/kmeans_warp", subj, sess, "dwi")
    out_dir = os.path.join(proj_dir, "derivatives/afni", subj, sess, "dwi")
    for h_dir in [work_dir, out_dir]:
        if not os.path.exists(h_dir):
            os.makedirs(h_dir)

    # get label definitions
    with open(
        os.path.join(proj_dir, "code/amygdala_kmeans/initial_clustering.json")
    ) as jf:
        label_defs = json.load(jf)

    # set up file dict
    file_dict = {}
    file_dict["bla-num"] = label_defs[subj.split("-")[1]]["bla"][0]
    file_dict["ntv-b0"] = glob.glob(f"{b0_dir}/{subj}_{sess}_*desc-eddyCorrected*")[0]
    file_dict["ntv-mask"] = glob.glob(
        f"{kmean_dir}/{subj}_{sess}_*desc-kmeansCluster*"
    )[0]
    file_dict["mni-atlas"] = os.path.join(
        atlas_dir, "tpl-MNIPediatricAsym_cohort-5_res-2_T2w.nii.gz"
    )
    file_dict["epi-ref"] = os.path.join(
        epi_dir,
        f"{subj}_{sess}_{task}_run-1_space-MNIPediatricAsym_cohort-5_res-2_desc-scaled_bold.nii.gz",
    )

    file_dict = ants_warp(subj, work_dir, file_dict)
    file_dict = make_mask(subj, sess, work_dir, out_dir, file_dict)


if __name__ == "__main__":
    main()
