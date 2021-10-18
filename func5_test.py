# %%
import os
import json
import fnmatch
import subprocess
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from func2_finish_preproc import func_sbatch


# %%
def make_gmInt_mask(
    subj_list, deriv_dir, sess, phase, tplflow_dir, tplflow_str, frac_value, group_dir
):
    """ Make a gray matter * group intersection mask """

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
def get_subbrick(subj_file, beh):
    """ Find sub-brick of beh """

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
    return h_out.decode("utf-8").replace("\n", "")


# %%
def get_pars(df_group, subj):
    """ Determine subj pars group tertile"""

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
    """ Make seed from supplied coordinate """

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
def check_coord(coord_dict, group_dir, subj_int_mask):
    """ Return True if meaninful data exists at all checked coordinates """

    for seed in coord_dict:
        mask = os.path.join(group_dir, f"Check_{seed}.nii.gz")
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
    beh_list,
    glt_dict,
    subj_list,
    sess,
    phase,
    group_dir,
    deriv_dir,
    df_group,
    coord_dict,
):

    data_table = []
    for subj in subj_list:
        subj_file = os.path.join(
            deriv_dir, subj, sess, f"{phase}_decon_stats_REML+tlrc"
        )

        # check that data exists at desired coords
        subj_int_mask = os.path.join(deriv_dir, subj, sess, "mask_epi_anat+tlrc")
        data_exists = check_coord(coord_dict, group_dir, subj_int_mask)
        if not data_exists:
            continue

        subj_pars_group = get_pars(df_group, subj)
        for beh in beh_list:
            subj_brick = get_subbrick(subj_file, beh)
            if subj_brick:
                data_table.append(subj)
                data_table.append(subj_pars_group)
                data_table.append(beh)
                data_table.append(f"""'{subj_file}[{int(subj_brick)}]'""")

    glt_list = []
    for count, test in enumerate(glt_dict):
        glt_list.append(f"-gltLabel {count + 1} {test}")
        glt_list.append(f"-gltCode {count + 1}")
        glt_list.append(f"'WSVARS: 1*{glt_dict[test][0]} -1*{glt_dict[test][1]}'")

    h_cmd = f"""
        cd {group_dir}

        3dMVM \
            -prefix MVM \
            -jobs 10 \
            -mask Group_GM_intersect_mask+tlrc \
            -bsVars PARS6 \
            -wsVars 'WSVARS' \
            -num_glt {len(list(glt_dict.keys()))} \
            {" ".join(glt_list)} \
            -dataTable \
            Subj PARS6 WSVARS InputFile \
            {" ".join(data_table)}
    """
    func_sbatch(h_cmd, 2, 6, 10, "cMVM", group_dir)


# %%
def model_noise(subj, subj_file, group_dir, acf_file):
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

    df_acf = pd.read_csv(acf_file, sep=" ", header=None)
    df_acf = df_acf.dropna(axis=1)
    df_acf = df_acf.loc[(df_acf != 0).any(axis=1)]
    mean_list = list(df_acf.mean())

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

    #
    deriv_dir = os.path.join(parent_dir, "derivatives/afni")
    group_dir = os.path.join(deriv_dir, "analyses")
    if not os.path.exists(group_dir):
        os.makedirs(group_dir)
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    # with open(os.path.join(group_dir, "beh_dict.json")) as json_file:
    #     beh_dict = json.load(json_file)

    # with open(os.path.join(group_dir, "glt_dict.json")) as json_file:
    #     glt_dict = json.load(json_file)

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
        beh_list = ["negLF", "negLC", "neuLF", "neuLC"]
        glt_dict = {
            "negLC-negLF": ["negLC", "negLF"],
            "negLF-neuLF": ["negLF", "neuLF"],
        }

        # TODO pull this csv from SharePoint
        df_sum = pd.read_csv("emuR01_summary_latest.csv")
        df_group = df_sum[["emu_study_id", "pinf_random", "pars_6"]]

        coord_dict = {"LAmg": "-24 -5 -29", "RAmg": "22 -3 -29"}
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
