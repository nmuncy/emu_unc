import os
import fnmatch

proj_dir = "/scratch/madlab/emu_UNC"
prep_dir = os.path.join(proj_dir, "derivatives/fmriprep")
afni_dir = os.path.join(proj_dir, "derivatives/afni")

if not os.path.exists(work_dir):
    os.makedirs(work_dir)

subj_list = [
    x.split(".")[0] for x in os.listdir(prep_dir) if fnmatch.fnmatch(x, "*.html")
]
