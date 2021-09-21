# %%
import os
import subprocess
import time
import fnmatch
import math
import shutil


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
# Submit jobs to slurm, wait for job to finish
def func_sbatch(command, wall_hours, mem_gig, num_proc, h_str, work_dir):

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
def func_scale(work_dir, phase_list, subj_num):

    # """
    # Step 10: Scale data

    # 1) Data is scaled by mean signal
    # """

    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        if not os.path.exists(os.path.join(work_dir, f"{phase}_minVal_mask+tlrc.HEAD")):
            h_cmd = f"""
                cd {work_dir}

                3dMean \
                    -datum short \
                    -prefix tmp_mean_{phase} \
                    tmp_run-{{1..{len(epi_list)}}}_{phase}_min+tlrc

                3dcalc \
                    -a tmp_mean_{phase}+tlrc \
                    -expr 'step(a-0.999)' \
                    -prefix {phase}_minVal_mask
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}min", work_dir)

    for phase in phase_list:
        epi_list = func_epi_list(phase, work_dir)
        for run in epi_list:
            if not os.path.exists(os.path.join(work_dir, f"{run}_scale+tlrc.HEAD")):
                h_cmd = f"""
                    cd {work_dir}

                    3dTstat -prefix tmp_tstat_{run} {run}_blur+tlrc

                    3dcalc -a {run}_blur+tlrc \
                        -b tmp_tstat_{run}+tlrc \
                        -c {phase}_minVal_mask+tlrc \
                        -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
                        -prefix {run}_scale
                """
                func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}scale", work_dir)


# %%
def copy_data(prep_dir, work_dir, subj, sess, task, num_runs):

    # set vars, dict
    tpl_ref = "space-MNIPediatricAsym_cohort-5_res-2"
    anat_str = f"{subj}_{sess}_{tpl_ref}"
    func_str = f"{subj}_{sess}_task-{task}"
    copy_dict = {
        "anat": {
            "tmp_struct_ns.nii.gz": f"{anat_str}_desc-preproc_T1w.nii.gz",
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
        if not os.path.exists(out_file):
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
def make_masks(work_dir, subj_num, atropos_dict, atropos_dir):

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

            # 3dresample \
            #     -master tmp_mask_allRuns+tlrc \
            #     -input struct_ns+tlrc \
            #     -prefix tmp_anat_resamp

            # 3dmask_tool \
            #     -dilate_input 5 -5 \
            #     -fill_holes \
            #     -input tmp_anat_resamp+tlrc \
            #     -prefix tmp_mask_struct

            3dmask_tool \
                -input tmp_mask_allRuns+tlrc struct_mask+tlrc \
                -inter \
                -prefix mask_epi_anat

            # 3dABoverlap \
            #     -no_automask tmp_mask_allRuns+tlrc tmp_mask_struct+tlrc | \
            #     tee out.mask_ae_overlap.txt
        """
        func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}uni", work_dir)

    # Make tissue-class masks
    #   I like Atropos better than AFNI's way, so use those priors
    h_ref = f"{epi_list[0]}_blur+tlrc"

    for key in atropos_dict:
        h_tiss = atropos_dict[key]
        if not os.path.exists(
            os.path.join(work_dir, f"final_mask_{h_tiss}_eroded+tlrc.HEAD")
        ):
            h_cmd = f"""
                module load c3d-1.0.0-gcc-8.2.0
                cd {work_dir}

                c3d \
                    {atropos_dir}/Prior{key}.nii.gz \
                    -thresh 0.3 1 1 0 \
                    -o tmp_{h_tiss}_bin.nii.gz

                3dresample \
                    -master {h_ref} \
                    -rmode NN \
                    -input tmp_{h_tiss}_bin.nii.gz \
                    -prefix final_mask_{h_tiss}+tlrc

                3dmask_tool \
                    -input tmp_{h_tiss}_bin.nii.gz \
                    -dilate_input -1 \
                    -prefix tmp_mask_{h_tiss}_eroded

                3dresample \
                    -master {h_ref} \
                    -rmode NN \
                    -input tmp_mask_{h_tiss}_eroded+orig \
                    -prefix final_mask_{h_tiss}_eroded
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}atr", work_dir)


# %%
def main():

    # For testing
    proj_dir = "/scratch/madlab/emu_UNC"
    prep_dir = os.path.join(proj_dir, "derivatives/old_fmriprep")
    afni_dir = os.path.join(proj_dir, "derivatives/afni")

    subj = "sub-4020"
    sess = "ses-S2"
    task = "test"
    num_runs = 3

    work_dir = os.path.join(afni_dir, subj, sess)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # get fMRIprep data
    if not os.path.exists(os.path.join(work_dir, "struct_ns+tlrc.HEAD")):
        copy_data(prep_dir, work_dir, subj, sess, task, num_runs)

    # blur data
    subj_num = subj.split("-")[-1]
    blur_mult = 1.5
    if not os.path.exists(os.path.join(work_dir, f"run-1_{task}_blur+tlrc.HEAD")):
        blur_epi(work_dir, subj_num, blur_mult)

    # make subject intersection mask


if __name__ == "__main__":
    main()
