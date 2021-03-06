"""Title.

Desc.

"""

# %%
import os
import fnmatch
import subprocess


# %%
def subj_roi(roi, roi_pieces, subj_dir):
    """Title.

    Desc.

    Parameters
    ----------

    """
    # extract rois
    for h_roi in roi_pieces:
        h_cmd = f"""
            module load c3d-1.0.0-gcc-8.2.0

            c3d {subj_dir}/aparc.nii.gz \
                -thresh \
                    {roi_pieces[h_roi]} \
                    {roi_pieces[h_roi]} \
                    {roi_pieces[h_roi]} \
                    0 \
                -o {subj_dir}/tmp_{h_roi}_single.nii.gz
        """
        h_thr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_thr.wait()

    # add rois
    h_cmd = f"""
        module load c3d-1.0.0-gcc-8.2.0

        c3d {subj_dir}/tmp_*_single.nii.gz \
            -accum -add -endaccum \
            -o {subj_dir}/prior_{roi}.nii.gz
        rm {subj_dir}/tmp*
    """
    h_add = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_add.wait()


def control_subj_roi(subj_list, jlf_dir, fs_dir, proj_dir, sess, roi_dict):
    """Title.

    Desc.

    Parameters
    ----------

    Returns
    -------
    """
    jlf_dict = {}
    for subj in subj_list:

        # set up subj_dir
        subj_dir = os.path.join(jlf_dir, subj)
        if not os.path.exists(subj_dir):
            os.makedirs(subj_dir)

        # get t1, seg files
        subj_t1 = os.path.join(
            proj_dir, "dset", subj, sess, "anat", f"{subj}_{sess}_T1w.nii.gz"
        )
        subj_fs = os.path.join(fs_dir, subj, "mri", "aparc.a2009s+aseg.mgz")
        if os.path.exists(subj_t1) and os.path.exists(subj_fs):
            jlf_dict[subj] = [subj_t1]
            h_cmd = f"""
                module load freesurfer-7.1
                cp {subj_t1} {subj_dir}
                mri_convert {subj_fs} {subj_dir}/aparc.nii.gz
            """
            h_cp = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_cp.wait()
        else:
            continue

        # extract rois
        for roi in roi_dict:
            roi_file = f"{subj_dir}/prior_{roi}.nii.gz"
            if not os.path.exists(roi_file):
                roi_pieces = roi_dict[roi]
                subj_roi(roi, roi_pieces, subj_dir)

            if os.path.exists(roi_file):
                jlf_dict[subj].append(roi_file)

    return jlf_dict


# %%
def main():
    # For testing
    sess = "ses-S2"
    roi_dict = {"LRAmg": {"LAmg": 18, "RAmg": 54}}
    proj_dir = "/scratch/madlab/emu_UNC"

    # set up
    fs_dir = os.path.join(proj_dir, "derivatives/freesurfer")
    jlf_dir = os.path.join(proj_dir, "derivatives/jlf")
    if not os.path.exists(jlf_dir):
        os.makedirs(jlf_dir)
    subj_list = [x for x in os.listdir(fs_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    # %%
    # make subject priors
    jlf_dict = control_subj_roi(
        subj_list[0:5], jlf_dir, fs_dir, proj_dir, sess, roi_dict
    )


if __name__ == "__main__":
    main()
