#!/bin/bash

# pyAFQ requires ipywidgets

proj_dir=/scratch/madlab/emu_UNC
diff_dir=${proj_dir}/derivatives/dwi_preproc

# # get json AFQ needs
# ref_file=${proj_dir}/dset/dataset_description.json
# cp $ref_file $proj_dir
# cp $ref_file $diff_dir

# submit job
python diff1_prob_setup.py $proj_dir $diff_dir
