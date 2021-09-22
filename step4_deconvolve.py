# %%
import os
import subprocess
import pandas as pd
import fnmatch


# %%
def mot_files(work_dir, num_run):

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

    for run in range(0, num_run):

        h_run = f"run-{run + 1}"
        df_all = pd.read_csv(
            os.path.join(work_dir, f"{h_run}_motion_all.tsv"), sep="\t"
        )

        df_mean = df_all[mean_labels].copy()
        df_mean = df_mean.round(6)
        df_mean.to_csv(
            os.path.join(work_dir, f"mot_demean_{task}_{h_run}.1D"),
            sep=" ",
            index=False,
            header=False,
        )

        df_drv = df_all[drv_labels].copy()
        df_drv = df_drv.fillna(0)
        df_drv = df_drv.round(6)
        df_drv.to_csv(
            os.path.join(work_dir, f"mot_deriv_{task}_{h_run}.1D"),
            sep=" ",
            index=False,
            header=False,
        )

        df_out = df_all.filter(regex="motion_outlier")
        df_out["sum"] = df_out.sum(axis=1)
        df_out = df_out.astype(int)
        df_out.to_csv(
            os.path.join(work_dir, f"mot_censor_{task}_{h_run}.1D"),
            sep=" ",
            index=False,
            header=False,
            columns=["sum"],
        )

    # cat files into singular rather than pad zeros
    #   bash for simplicity
    h_cmd = f"""
        cd {work_dir}

        > mot_demean_{task}_concat.1D
        > mot_deriv_{task}_concat.1D
        > mot_censor_{task}_concat.1D

        cat mot_demean_{task}_run*.1D >> mot_demean_{task}_concat.1D
        cat mot_deriv_{task}_run*.1D >> mot_deriv_{task}_concat.1D
        cat mot_censor_{task}_run*.1D >> mot_censor_{task}_concat.1D
    """
    h_cat = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_cat.wait()


def func_write_decon(task, work_dir, epi_list, dur, tf_dict, decon_type):
    """
    Notes:
        This function generates a 3dDeconvolve command.
        It supports GAM, 2GAM, TENT, and dmBLOCK basis functions.
        TENT does not currently include duration.
    """

    # # build censor arguments
    # reg_base = []
    # for cmot, mot in enumerate(dmn_list):
    #     reg_base.append(f"-ortvec {mot} mot_dmn_{cmot + 1}")

    # for cmot, mot in enumerate(drv_list):
    #     reg_base.append(f"-ortvec {mot} mot_drv_{cmot + 1}")

    reg_base = [
        f"-ortvec mot_demean_{task}_concat.1D mot_dmn_1",
        f"-ortvec mot_deriv_{task}_concat.1D mot_drv_1",
    ]

    # determine tr
    h_cmd = f"""
        module load afni-20.2.06
        3dinfo -tr {work_dir}/{epi_list[0]}
    """
    h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_len_tr = h_tr.communicate()[0]
    len_tr = float(h_len_tr.decode("utf-8").strip())

    # determine, build behavior regressors
    switch_dict = {
        "dmBLOCK": ["'dmBLOCK(1)'", "-stim_times_AM1"],
        "GAM": ["'GAM'", "-stim_times"],
        "2GAM": [f"'TWOGAMpw(4,5,0.2,12,7,{dur})'", "-stim_times"],
    }

    reg_beh = []
    for c_beh, beh in enumerate(tf_dict):
        if decon_type == "dmBLOCK" or decon_type == "GAM" or decon_type == "2GAM":

            # add stim_time info, order is
            #   -stim_times 1 tf_beh.txt basisFunction
            reg_beh.append(switch_dict[decon_type][1])
            reg_beh.append(f"{c_beh + 1}")
            reg_beh.append(f"timing_files/{tf_dict[beh]}")
            reg_beh.append(switch_dict[decon_type][0])

            # add stim_label info, order is
            #   -stim_label 1 beh
            reg_beh.append("-stim_label")
            reg_beh.append(f"{c_beh + 1}")
            reg_beh.append(beh)

        elif decon_type == "TENT":

            # extract duration, account for no behavior in 1st run
            tmp_str = tf_dict[beh].replace("tf", "dur")
            dur_file = open(os.path.join(work_dir, "timing_files", tmp_str)).readlines()

            if "*" not in dur_file[0]:
                tent_len = round(12 + float(dur_file[0]))
            else:
                with open(os.path.join(work_dir, "timing_files", tmp_str)) as f:
                    for line in f:
                        s = re.search(r"\d+", line)
                        if s:
                            tmp_num = s.string.split("\t")[0]
                tent_len = round(12 + float(tmp_num))
            tent_args = ["0", str(tent_len), str(round(tent_len / len_tr))]

            # stim_time
            reg_beh.append("-stim_times")
            reg_beh.append(f"{c_beh + 1}")
            reg_beh.append(f"timing_files/{tf_dict[beh]}")
            reg_beh.append(f"""'TENT({",".join(tent_args)})'""")

            # stim_label
            reg_beh.append("-stim_label")
            reg_beh.append(f"{c_beh + 1}")
            reg_beh.append(beh)

    # set output str
    h_out = f"{phase}_{desc}"

    # build full decon command
    cmd_decon = f"""
        3dDeconvolve \\
            -x1D_stop \\
            -GOFORIT \\
            -input {" ".join(run_list)} \\
            -censor {cen_file} \\
            {" ".join(reg_base)} \\
            -polort A \\
            -float \\
            -local_times \\
            -num_stimts {len(tf_dict.keys())} \\
            {" ".join(reg_beh)} \\
            -jobs 1 \\
            -x1D X.{h_out}.xmat.1D \\
            -xjpeg X.{h_out}.jpg \\
            -x1D_uncensored X.{h_out}.nocensor.xmat.1D \\
            -bucket {h_out}_stats \\
            -cbucket {h_out}_cbucket \\
            -errts {h_out}_errts
    """
    print(cmd_decon)
    return cmd_decon


# %%
# For testing
afni_dir = "/scratch/madlab/emu_UNC/derivatives/afni"
subj = "sub-4002"
sess = "ses-S2"
task = "test"
num_run = 3

work_dir = os.path.join(afni_dir, subj, sess)

# motion and censor
if not os.path.exists(os.path.join(work_dir, f"mot_censor_{task}_concat.1D")):
    mot_files(work_dir, num_run)

# write decon script
tf_list = os.listdir(os.path.join(work_dir, "timing_files"))
tf_list.sort()
tf_dict = {}
for tf in tf_list:
    beh = tf.split("_")[-1].split(".")[0]
    tf_dict[beh] = tf

epi_list = [
    x.split(".")[0]
    for x in os.listdir(work_dir)
    if fnmatch.fnmatch(x, "*scale+tlrc.HEAD")
]
epi_list.sort()


# %%
