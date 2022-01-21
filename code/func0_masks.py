"""Title

Desc.

Examples
--------
"""

# %%
import os
import sys
import json
import glob
from func1_ppi import submit_hpc_sbatch


# %%
def ants_warp(subj, work_dir, file_dict):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """
    ants_warp = os.path.join(work_dir, "ants_Warped.nii.gz")
    if not os.path.exists(ants_warp):
        h_cmd = f"""
            module load ants-2.3.5
            antsRegistrationSyN.sh \
                -d 3 \
                -f {file_dict["mni-atlas"]} \
                -m {file_dict["ntv-b0"]} \
                -o ${work_dir}/ants_ \
                -n 4
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 2, 4, 4, f"{subj_num}ants", work_dir)

    assert os.path.exists(ants_warp), "ANTs Registration failed."
    file_dict["mni-b0"] = ants_warp
    file_dict["wrp-calc"] = os.path.join(work_dir, "ants_0GenericAffine.mat")
    return file_dict


# %%
def make_mask(subj, sess, work_dir, file_dict):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """

    bla_num = file_dict["bla-num"]
    space = file_dict["mni-atlas"].split("/")[-1].split("tpl-")[1].split("_")[0]
    cohort = file_dict["mni-atlas"].split("/")[-1].split("cohort-")[1].split("_")[0]
    bla_mask = os.path.join(
        work_dir,
        f"{subj}_{sess}_space-{space}_cohort-{cohort}_res-2_desc-bla_mask.nii.gz",
    )
    if not os.path.exists(bla_mask):
        h_cmd = f"""
            module load ants-2.3.5

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
        """
        subj_num = subj.split("-")[1]
        h_out, h_err = submit_hpc_sbatch(h_cmd, 1, 4, 1, f"{subj_num}wrp", work_dir)

    assert os.path.exists(bla_mask), f"Failed to make {bla_mask}"
    file_dict["mni-mask"] = bla_mask
    return file_dict


# %%
def main():

    # For testing
    proj_dir = "/home/data/madlab/McMakin_EMUR01"
    scratch_dir = "/scratch/madlab/emu_unc"
    subj = "sub-4001"
    sess = "ses-S2"
    task = "task-test"

    # check for correct conda env
    assert (
        "emuR01_unc_env" in sys.executable
    ), "Please activate emuR01_unc conda environment."

    # set up
    atlas_dir = "/home/data/madlab/atlases/templateflow/tpl-MNIPediatricAsym/cohort-5"
    kmean_dir = os.path.join(proj_dir, "derivatives/kmeans", subj, sess, "dwi")
    b0_dir = os.path.join(proj_dir, "derivatives/dwi_bedpostx", subj, "emubpx/b0_file")
    epi_dir = os.path.join(proj_dir, "dset", subj, sess, "func")
    work_dir = os.path.join(scratch_dir, "derivatives/kmeans_warp", subj, sess, "dwi")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

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
        epi_dir, f"{subj}_{sess}_{task}_run-1_bold.nii.gz"
    )

    file_dict = ants_warp(subj, work_dir, file_dict)
    file_dict = make_mask(subj, sess, work_dir, file_dict)


if __name__ == "__main__":
    main()
