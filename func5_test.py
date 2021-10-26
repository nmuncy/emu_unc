"""Conduct Group-level Analyses.

Desc.

Examples
--------

Notes
-----

"""

# %%
import os
import json
import fnmatch
import subprocess
import pathlib
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from func2_finish_preproc import func_sbatch


# %%
def make_gmInt_mask(
    subj_list,
    deriv_dir,
    sess,
    phase,
    tplflow_dir,
    tplflow_str,
    group_dir,
    frac_value=0.8,
):
    """Make a gray matter * group intersection mask.

    Create a group intersection mask by determining where meaningful data
    exists in both functional and structural data across the entire project.
    Then restrict this mask by multipying it by a gray matter mask. This
    helps reduce the number of voxels to run statistics and simulations on.

    Parameters
    ----------
    subj_list : list
        list of subjects to include in group intersection mask
    deriv_dir : str
        /path/to/BIDS/derivatives/afni
    sess : str
        BIDS session string (ses-S1)
    phase : str
        aspect of experiment being analysed (study, test)
    tplflow_dir : str
        /path/to/templateflow/atlas_dir
    tplflow_str : str
        atlas file string (tpl-MNIPediatricAsym_cohort-5_res-1)
        used to find T1 and GM files
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    frac_value = float [default=0.8]
        proportion of participants that must have data in a voxel for
        the voxel to be tested at the group-level.

    Notes
    -----
    tplflow_dir needs to contain tissue class segmentation
        masks (label-GM_probseg)

    Output written to group_dir

    MRI output : structural
        Group_intersect_mean.nii.gz
        Group_intersect_mask.nii.gz
        Group_GM_intersect_mask.nii.gz
        Group_GM_intersect_mask+tlrc
    """

    # set ref file for resampling
    ref_file = os.path.join(deriv_dir, subj_list[1], sess, f"run-1_{phase}_scale+tlrc")

    # make group intersection mask
    if not os.path.exists(os.path.join(group_dir, "Group_intersect_mask.nii.gz")):

        mask_list = []
        for subj in subj_list:
            mask_file = os.path.join(deriv_dir, subj, sess, "mask_epi_anat+tlrc")
            if os.path.exists(f"{mask_file}.HEAD"):
                mask_list.append(mask_file)

        h_cmd = f"""
            module load afni-20.2.06
            cd {group_dir}

            cp {tplflow_dir}/{tplflow_str}_T1w.nii.gz .
            3dMean -prefix Group_intersect_mean.nii.gz {" ".join(mask_list)}
            3dmask_tool \
                -input {" ".join(mask_list)} \
                -frac {frac_value} \
                -prefix Group_intersect_mask.nii.gz
        """
        h_mask = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_mask.wait()

    # make GM intersection mask
    if not os.path.exists(os.path.join(group_dir, "Group_GM_intersect_mask+tlrc.HEAD")):
        h_cmd = f"""
            module load afni-20.2.06
            module load c3d-1.0.0-gcc-8.2.0
            cd {group_dir}

            c3d \
                {tplflow_dir}/{tplflow_str}_label-GM_probseg.nii.gz \
                -thresh 0.3 1 1 0 \
                -o tmp_GM.nii.gz

            3dresample \
                -master {ref_file} \
                -rmode NN \
                -input tmp_GM.nii.gz \
                -prefix tmp_GM_mask.nii.gz

            c3d \
                tmp_GM_mask.nii.gz Group_intersect_mask.nii.gz \
                -multiply \
                -o tmp_GM_intersect_prob_mask.nii.gz

            c3d \
                tmp_GM_intersect_prob_mask.nii.gz \
                -thresh 0.1 1 1 0 \
                -o Group_GM_intersect_mask.nii.gz

            3dcopy Group_GM_intersect_mask.nii.gz Group_GM_intersect_mask+tlrc
            3drefit -space MNI Group_GM_intersect_mask+tlrc

            if [ -f Group_GM_intersect_mask+tlrc.HEAD ]; then
                rm tmp_*
            fi
        """
        h_GMmask = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_GMmask.wait()


# %%
def get_subbrick(subj_file, beh_list):
    """Find sub-brick of a behavior.

    As every participant may not have every behavior (targMS),
    deconvolvution occured with available behaviors and so
    do not have a consistant sub-brick assignment. This will find
    the sub-brick associated with the desired behavior (beh).

    Parameters
    ----------
    subj_file : str
        deconvolved subject file (test_decon_stats_REML+tlrc)
    # beh : str
    #     behavior string (lureCR)

    Returns
    -------
    # h_out : str
    #     string of sub-brick integer
    """

    subbrick_dict = {}
    for beh in beh_list:
        h_cmd = f"""
            module load afni-20.2.06
            3dinfo \
                -subbrick_info {subj_file} |
                grep "{beh}#0_Coef" |
                awk '{{print $4}}' |
                sed 's/\#//'
        """
        h_brick = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_brick.wait()
        h_out = h_brick.communicate()[0]
        subbrick = h_out.decode("utf-8").replace("\n", "")
        if subbrick:
            subbrick_dict[beh] = subbrick

    return subbrick_dict


# %%
def get_pars(df_group, subj):
    """Determine subject PARS6 group tertile.

    Parameters
    ----------
    df_group : pandas.DataFrame
        summary emuR01 dataframe with PARS6 scores
    subj : str
        BIDS subject

    Returns
    -------
    subj_pars_group : str
        subject Low, Med, High grouping
    """

    subj_num = int(subj.split("-")[-1])
    idx_subj = df_group.index[df_group["emu_study_id"] == subj_num]
    subj_pars = df_group.iloc[idx_subj]["pars_6"].item()
    if subj_pars <= 3.0:
        subj_pars_group = "Low"
    elif subj_pars > 3.0 and subj_pars <= 12.0:
        subj_pars_group = "Med"
    elif subj_pars > 12.0:
        subj_pars_group = "High"
    else:
        subj_pars_group = np.nan
    return subj_pars_group


# %%
def make_seed(seed_name, coord_list, group_dir, ref_file):
    """Make seed from supplied coordinates.

    Make a ROI seed of size=1voxel from a set of coordinates.
    Used to test if a subject has EPI data at the desired
    coordinate (account for subject-specific fallout).

    Parameters
    ----------
    seed_name : str
        name of seed ROI (LHC)
    coord_list : str
        coordinates (-25 7 48)
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    ref_file : str
        /path/to/BIDS/derivatives/afni/sub-1234/ses-A/run-1_scale+tlrc
        used for writing header

    Notes
    -----
    MRI output : structural
        group_dir/Check_<seed_name>.nii.gz
    """

    h_cmd = f"""
        module load afni-20.2.06
        echo "{coord_list}" > tmp_coord
        3dUndump \
            -prefix {group_dir}/Check_{seed_name}.nii.gz \
            -master {ref_file} \
            -srad 0 \
            -xyz tmp_coord
        rm tmp_coord
    """
    h_seed = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_seed.wait()


# %%
def check_coord(coord_dict, group_dir, subj_int_mask, subj):
    """Return True if meaningful data exists at all checked coordinates.

    Parameters
    ----------
    coord_dict : dict
        dictionary of {roi: coordinates}
        e.g. {"LAmg": "-24 -5 -29"}
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    subj_int_mask : str
        subject intersection mask
        /path/to/BIDS/derivatives/afni/sub-1234/ses-A/mask_epi_anat+tlrc
    subj : str
        BIDS subject string

    Returns
    -------
    out_status : bool
        True if mask value == 1.0 at seed location
    """

    for seed in coord_dict:
        mask = os.path.join(group_dir, f"Check_{seed}.nii.gz")
        print(f"Checking {subj} for data in {seed} ...")
        if not os.path.exists(mask):
            make_seed(seed, coord_dict[seed], group_dir, subj_int_mask)
        h_cmd = f"""
            module load afni-20.2.06
            3dROIstats \
                -quiet \
                -mask {mask} \
                {subj_int_mask} |
                tr -d '[:space:]'
        """
        h_check = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        h_check.wait()
        vox_value = h_check.communicate()[0]
        out_status = True if float(vox_value.decode("utf-*")) == 1.0 else False
        if not out_status:
            return out_status
    return out_status


# %%
def group_analysis(
    mvm_dict,
    subj_list,
    sess,
    phase,
    group_dir,
    deriv_dir,
    df_group,
):
    """Conduct group-level analyses.

    Write a 3dMVM command by building the dataTable of subjects who
    have data at pre-specified coordinates.

    Parameters
    ----------
    # beh_list = list
    #     list of behaviors to test
    # glt_dict = dict
    #     dictionary of post-hoc t-tests
    subj_list = list
        list of subjects to potentially include
    sess = str
        BIDS session string
    phase : str
        aspect of experiment being analysed (study, test)
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    deriv_dir : str
        /path/to/BIDS/derivatives/afni
    df_group : pandas.DataFrame
        summary emuR01 dataframe with PARS6 scores
    # coord_dict : dict
    #     dictionary of {roi: coordinates}
    #     e.g. {"LAmg": "-24 -5 -29"}

    Notes
    -----
    MRI output : functional
        group_dir/MVM_<phase>+tlrc
    """

    # create data_table of subjects who have data at coordinates
    # and have the behaviors of interest
    for mvm_title in mvm_dict:

        # get planned info
        beh_list = mvm_dict[mvm_title]["behaviors"]
        glt_dict = mvm_dict[mvm_title]["post-hoc"]
        coord_dict = mvm_dict[mvm_title]["coord_check"]
        decon_file = mvm_dict[mvm_title]["decon_file"]

        data_table = []
        mvm_subj_dict = {}
        for subj in subj_list:

            # check that data exists at desired coords
            subj_int_mask = os.path.join(deriv_dir, subj, sess, "mask_epi_anat+tlrc")
            data_exists = check_coord(coord_dict, group_dir, subj_int_mask, subj)
            if not data_exists:
                continue

            # determine group
            subj_pars_group = get_pars(df_group, subj)

            # find behavior sub-brick, write subject row if subject
            # has all behaviors
            subj_file = os.path.join(deriv_dir, subj, sess, decon_file)
            beh_dict = get_subbrick(subj_file, beh_list)
            if len(beh_dict.keys()) != len(beh_list):
                continue
            else:
                mvm_subj_dict[subj] = beh_dict
                for beh in beh_dict:
                    data_table.append(subj)
                    data_table.append(subj_pars_group)
                    data_table.append(beh)
                    data_table.append(f"""'{subj_file}[{beh_dict[beh]}]'""")

        # write out subjects included in mvm
        with open(os.path.join(group_dir, f"subj_{mvm_title}.json"), "w") as jf:
            json.dump(mvm_subj_dict, jf)

        # set up post-hoc tests
        glt_list = []
        for count, test in enumerate(glt_dict):
            glt_list.append(f"-gltLabel {count + 1} {test}")
            glt_list.append(f"-gltCode {count + 1} {glt_dict[test]}")

        # write command for review
        h_cmd = f"""
            3dMVM \\
                -prefix MVM_{mvm_title} \\
                -jobs 10 \\
                -mask Group_GM_intersect_mask+tlrc \\
                -bsVars PARS6 \\
                -wsVars 'WSVARS' \\
                -num_glt {len(glt_dict.keys())} \\
                {" ".join(glt_list)} \\
                -dataTable \\
                Subj PARS6 WSVARS InputFile \\
                {" ".join(data_table)}
        """
        mvm_script = os.path.join(group_dir, f"MVM_{mvm_title}.sh")
        with open(mvm_script, "w") as ms:
            ms.write(h_cmd)

    # execute scripts
    mvm_list = [x for x in os.listdir(group_dir) if fnmatch.fnmatch(x, "MVM_*.sh")]
    for mvm_script in mvm_list:
        h_cmd = f"""
            cd {group_dir}
            source {mvm_script}
        """
        func_sbatch(h_cmd, 2, 6, 10, "cMVM", group_dir)


# %%
def model_noise(subj, subj_file, group_dir, acf_file):
    """Model noise.

    Parameters
    ----------
    subj : str
        BIDS subject string
    subj_file : str
        deconvolution residual of a subject
        /path/to/BIDS/derivatives/afni/sub-1234/ses-A/<phase>_decon_errts_REML+tlrc
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    acf_file : str
        output acf txt to append (group_dir/ACF_subj_all.txt)

    Notes
    -----
    Parameter estimations are captured in acf_file.
    """

    h_cmd = f"""
        cd {group_dir}

        3dFWHMx \
            -mask Group_GM_intersect_mask+tlrc \
            -input {subj_file} \
            -acf >> {acf_file}
    """
    func_sbatch(h_cmd, 2, 4, 1, f"a{subj.split('-')[-1]}", group_dir)


# %%
def run_montecarlo(group_dir, acf_file, mc_file):
    """Conduct monte carlo simulations.

    Incorporate noise estimations in MC simulations. Restrict
    simulations to group gray matter * intersection mask.

    Parameters
    ----------
    group_dir : str
        /path/to/BIDS/derivatives/afni/analyses
    acf_file : str
        output of model_noise (group_dir/ACF_subj_all.txt)
        used to sum across 3 parameters
    mc_file : str
        captures simulation output/recommendations
        group_dir/MC_thresholds.txt

    Notes
    -----
    3dClustSim output are captured in group_dir/mc_file.
    """

    # average relevant parameter estimations of 3dFWHMx
    df_acf = pd.read_csv(acf_file, sep=" ", header=None)
    df_acf = df_acf.dropna(axis=1)
    df_acf = df_acf.loc[(df_acf != 0).any(axis=1)]
    mean_list = list(df_acf.mean())

    # conduct MC simulations
    h_cmd = f"""
        cd {group_dir}

        3dClustSim \
            -mask Group_GM_intersect_mask+tlrc \
            -LOTS \
            -iter 10000 \
            -acf {mean_list[0]} {mean_list[1]} {mean_list[2]} \
            > {mc_file}
    """
    func_sbatch(h_cmd, 6, 4, 10, "mc", group_dir)


# %%
def get_args():
    parser = ArgumentParser("Receive Bash args from wrapper")
    parser.add_argument("sess", help="Session")
    parser.add_argument("phase", help="Phase")
    parser.add_argument(
        "tplflow_dir", help="Location of fMRIprep output-spaces template"
    )
    parser.add_argument("parent_dir", help="Location of Project Directory")

    return parser


# %%
def main():

    # TODO update for multiple MVMs

    # for testing
    parent_dir = "/scratch/madlab/emu_UNC"
    tplflow_dir = "/home/data/madlab/singularity-images/templateflow/tpl-MNIPediatricAsym/cohort-5"
    tplflow_str = "tpl-MNIPediatricAsym_cohort-5_res-1"
    phase = "test"
    sess = "ses-S2"

    # # get args
    # args = get_args().parse_args()
    # sess = args.sess
    # phase = args.phase
    # tplflow_dir = args.tplflow_dir
    # parent_dir = args.parent_dir

    code_dir = pathlib.Path().resolve()

    #
    deriv_dir = os.path.join(parent_dir, "derivatives/afni")
    group_dir = os.path.join(deriv_dir, "analyses")
    if not os.path.exists(group_dir):
        os.makedirs(group_dir)
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    with open(os.path.join(code_dir, "mvm_plan.json")) as json_file:
        mvm_dict = json.load(json_file)

    """ make group gray matter intersect mask """
    if not os.path.exists(os.path.join(group_dir, "Group_GM_intersect_mask+tlrc.HEAD")):
        frac_value = 0.8
        make_gmInt_mask(
            subj_list,
            deriv_dir,
            sess,
            phase,
            tplflow_dir,
            tplflow_str,
            frac_value,
            group_dir,
        )

    """ run MVM """
    if not os.path.exists(os.path.join(group_dir, "MVM+tlrc.HEAD")):

        # TODO pull this csv from SharePoint
        df_sum = pd.read_csv("emuR01_summary_latest.csv")
        df_group = df_sum[["emu_study_id", "pinf_random", "pars_6"]]

        # coord_dict = {"LAmg": "-24 -5 -29", "RAmg": "22 -3 -29"}
        # group_analysis(beh_dict, glt_dict, subj_list, sess, phase, group_dir, deriv_dir)

    """ get subj acf estimates """
    # define, start file
    acf_file = os.path.join(group_dir, "ACF_subj_all.txt")
    if not os.path.exists(acf_file):
        open(acf_file, "a").close()

    # if file is empty, run model_noise for e/subj
    acf_size = os.path.getsize(acf_file)
    if acf_size == 0:
        for subj in subj_list:
            subj_file = os.path.join(
                deriv_dir, subj, sess, f"{phase}_single_errts_REML+tlrc"
            )
            model_noise(subj, subj_file, group_dir, acf_file)

    """ do clust simulations """
    mc_file = os.path.join(group_dir, "MC_thresholds.txt")
    if not os.path.exists(mc_file):
        open(mc_file, "a").close()

    mc_size = os.path.getsize(mc_file)
    if mc_size == 0:
        run_montecarlo(group_dir, acf_file, mc_file)


if __name__ == "__main__":
    main()
