# %%
import os
import glob
import subprocess
import time
import fnmatch
import math
import shutil
from argparse import ArgumentParser


# %%
def _copyfile_patched(src, dst, length=16 * 1024 * 1024):
    """Patches shutil method to improve copyfile speed"""
    while 1:
        buf = src.read(length)
        if not buf:
            break
        dst.write(buf)
    shutil.copyfile = _copyfile_patched


# %%
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_str, work_dir):
    """ Submit jobs to slurm, wait for job to finish """
    full_name = f"{work_dir}/sbatch_writeOut_{h_str}"
    sbatch_job = f"""
        sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p IB_44C_512G -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load afni-20.2.06 \n {command}"
    """
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_response.communicate()[0]
    print(job_id, h_str, sbatch_job)

    while_count = 0
    break_status = False
    while not break_status:

        check_cmd = "squeue -u $(whoami)"
        sq_check = subprocess.Popen(check_cmd, shell=True, stdout=subprocess.PIPE)
        out_lines = sq_check.communicate()[0]
        b_decode = out_lines.decode("utf-8")

        if h_str not in b_decode:
            break_status = True
        else:
            while_count += 1
            print(f"Wait count for sbatch job {h_str}: {while_count}")
            time.sleep(3)

    print(f"Sbatch job {h_str} finished")


# %%
def copy_data(prep_dir, work_dir, subj, sess, task, num_runs):
    """ get relevant fmriprep files, rename """
    # set vars, dict
    tpl_ref = "space-MNIPediatricAsym_cohort-5_res-2"
    anat_str = f"{subj}_{sess}_{tpl_ref}"
    func_str = f"{subj}_{sess}_task-{task}"
    copy_dict = {
        "anat": {
            "tmp_struct_head.nii.gz": f"{anat_str}_desc-preproc_T1w.nii.gz",
            "tmp_struct_mask.nii.gz": f"{anat_str}_desc-brain_mask.nii.gz",
            "final_mask_GM_prob.nii.gz": f"{anat_str}_label-GM_probseg.nii.gz",
            "final_mask_WM_prob.nii.gz": f"{anat_str}_label-WM_probseg.nii.gz",
        },
        "func": {},
    }

    # add preproc bold, confound TS to copy_dicts
    for run in range(0, num_runs):
        h_run = f"run-{run+1}"

        copy_dict["func"][
            f"tmp_{h_run}_{task}_preproc.nii.gz"
        ] = f"{func_str}_{h_run}_{tpl_ref}_desc-preproc_bold.nii.gz"

        copy_dict["func"][
            f"{h_run}_motion_all.tsv"
        ] = f"{func_str}_{h_run}_desc-confounds_timeseries.tsv"

    # copy data
    for scan_type in copy_dict:
        source_dir = os.path.join(prep_dir, subj, sess, scan_type)
        for h_file in copy_dict[scan_type]:
            in_file = os.path.join(source_dir, copy_dict[scan_type][h_file])
            out_file = os.path.join(work_dir, h_file)
            if not os.path.exists(os.path.join(work_dir, h_file)):
                shutil.copyfile(in_file, out_file)

    # 3dcopy data
    tmp_list = [x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, "tmp_*")]
    for tmp_file in tmp_list:
        in_file = os.path.join(work_dir, tmp_file)
        h_str = tmp_file.split(".")[0].split("_", 1)[1]
        out_file = os.path.join(work_dir, f"{h_str}+tlrc")
        if not os.path.exists(f"{out_file}.HEAD"):
            h_cmd = f"""
                module load afni-20.2.06
                3dcopy {in_file} {out_file} && rm {in_file}
            """
            h_cp = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_cp.wait()


# %%
def mk_epi_list(work_dir):
    epi_list = [
        x.split("_pre")[0]
        for x in os.listdir(work_dir)
        if fnmatch.fnmatch(x, "*preproc+tlrc.HEAD")
    ]
    epi_list.sort()
    return epi_list


# %%
def blur_epi(work_dir, subj_num, blur_mult):

    # Blur
    epi_list = mk_epi_list(work_dir)

    for run in epi_list:
        if not os.path.exists(os.path.join(work_dir, f"{run}_blur+tlrc.HEAD")):

            # calc voxel dim i
            h_cmd = f"""
                module load afni-20.2.06
                3dinfo -dk {work_dir}/{run}_preproc+tlrc
            """
            h_gs = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_gs_out = h_gs.communicate()[0]
            grid_size = h_gs_out.decode("utf-8").strip()
            blur_size = math.ceil(blur_mult * float(grid_size))

            # do blur
            h_cmd = f"""
                cd {work_dir}
                3dmerge \
                    -1blur_fwhm {blur_size} \
                    -doall \
                    -prefix {run}_blur \
                    {run}_preproc+tlrc
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}blur", work_dir)


# %%
def make_masks(work_dir, subj_num, tiss_list):

    epi_list = mk_epi_list(work_dir)

    # Make EPI-T1 intersection mask (mask_epi_anat)
    if not os.path.exists(os.path.join(work_dir, "mask_epi_anat+tlrc.HEAD")):

        for run in epi_list:
            if not os.path.exists(
                os.path.join(work_dir, f"tmp_mask.{run}_blur+tlrc.HEAD")
            ):
                h_cmd = f"""
                    module load afni-20.2.06
                    cd {work_dir}
                    3dAutomask -prefix tmp_mask.{run} {run}_blur+tlrc
                """
                h_mask = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
                h_mask.wait()

        h_cmd = f"""
            cd {work_dir}

            3dmask_tool \
                -inputs tmp_mask.*+tlrc.HEAD \
                -union \
                -prefix tmp_mask_allRuns

            3dmask_tool \
                -input tmp_mask_allRuns+tlrc struct_mask+tlrc \
                -inter \
                -prefix mask_epi_anat
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}uni", work_dir)

    # Make tissue-class masks
    for tiss in tiss_list:
        if not os.path.exists(
            os.path.join(work_dir, f"final_mask_{tiss}_eroded+tlrc.HEAD")
        ):
            h_cmd = f"""
                module load c3d-1.0.0-gcc-8.2.0
                cd {work_dir}

                c3d \
                    final_mask_{tiss}_prob.nii.gz \
                    -thresh 0.5 1 1 0 \
                    -o tmp_{tiss}_bin.nii.gz

                3dmask_tool \
                    -input tmp_{tiss}_bin.nii.gz \
                    -dilate_input -1 \
                    -prefix final_mask_{tiss}_eroded

                3drefit -space MNI final_mask_{tiss}_eroded+orig
                3drefit -view tlrc final_mask_{tiss}_eroded+orig
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}tiss", work_dir)


# %%
def scale_epi(work_dir, subj_num, task):

    epi_list = mk_epi_list(work_dir)

    # TODO update here for multiple phases/tasks
    for run in epi_list:
        if not os.path.exists(os.path.join(work_dir, f"tmp_{run}_mask+tlrc.HEAD")):
            h_cmd = f"""
                module load afni-20.2.06
                cd {work_dir}

                3dcalc \
                    -overwrite \
                    -a {run}_preproc+tlrc \
                    -expr 1 \
                    -prefix tmp_{run}_mask

                3dTstat \
                    -min \
                    -prefix tmp_{run}_min \
                    tmp_{run}_mask+tlrc
            """
            h_min = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            h_min.wait()

    if not os.path.exists(os.path.join(work_dir, f"{task}_minVal_mask+tlrc.HEAD")):
        h_cmd = f"""
            cd {work_dir}

            3dMean \
                -datum short \
                -prefix tmp_mean_{task} \
                tmp_run-{{1..{len(epi_list)}}}_{task}_min+tlrc

            3dcalc \
                -a tmp_mean_{task}+tlrc \
                -expr 'step(a-0.999)' \
                -prefix {task}_minVal_mask
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}min", work_dir)

    # TODO update here for multiple phases/tasks
    for run in epi_list:
        if not os.path.exists(os.path.join(work_dir, f"{run}_scale+tlrc.HEAD")):
            h_cmd = f"""
                cd {work_dir}

                3dTstat -prefix tmp_tstat_{run} {run}_blur+tlrc

                3dcalc -a {run}_blur+tlrc \
                    -b tmp_tstat_{run}+tlrc \
                    -c {task}_minVal_mask+tlrc \
                    -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
                    -prefix {run}_scale
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}scale", work_dir)


def get_args():
    parser = ArgumentParser("Receive Bash args from wrapper")
    parser.add_argument("h_subj", help="Subject ID")
    parser.add_argument("h_sess", help="Session")
    parser.add_argument("h_task", help="Task String")
    parser.add_argument("h_num", help="Number of Runs")
    parser.add_argument("h_prep", help="/path/to/derivatives/fmriprep")
    parser.add_argument("h_afni", help="/path/to/derivatives/afni")
    return parser


# %%
def main():

    # # For testing
    # proj_dir = "/scratch/madlab/emu_UNC"
    # prep_dir = os.path.join(proj_dir, "derivatives/fmriprep")
    # afni_dir = os.path.join(proj_dir, "derivatives/afni")

    # subj = "sub-4020"
    # sess = "ses-S2"
    # task = "test"
    # num_runs = 3

    """ get passed arguments """
    args = get_args().parse_args()
    subj = args.h_subj
    sess = args.h_sess
    task = args.h_task
    num_runs = int(args.h_num)
    prep_dir = args.h_prep
    afni_dir = args.h_afni

    """ setup afni directory """
    work_dir = os.path.join(afni_dir, subj, sess)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # get fMRIprep data
    if not os.path.exists(os.path.join(work_dir, "struct_head+tlrc.HEAD")):
        copy_data(prep_dir, work_dir, subj, sess, task, num_runs)

    # blur data
    subj_num = subj.split("-")[-1]
    blur_mult = 1.5
    if not os.path.exists(os.path.join(work_dir, f"run-1_{task}_blur+tlrc.HEAD")):
        blur_epi(work_dir, subj_num, blur_mult)

    # make subject intersection mask
    tiss_list = ["GM", "WM"]
    if not os.path.exists(
        os.path.join(work_dir, f"final_mask_{tiss_list[0]}_eroded+tlrc.HEAD")
    ):
        make_masks(work_dir, subj_num, tiss_list)

    # scale data
    if not os.path.exists(os.path.join(work_dir, f"run-1_{task}_scale+tlrc.HEAD")):
        scale_epi(work_dir, subj_num, task)

    # clean
    if os.path.exists(os.path.join(work_dir, f"run-1_{task}_scale+tlrc.HEAD")):
        for tmp_file in glob.glob(f"{work_dir}/tmp*"):
            os.remove(tmp_file)


if __name__ == "__main__":
    main()
