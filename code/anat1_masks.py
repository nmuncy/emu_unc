"""Title

Desc.

Examples
--------
"""

# %%
import os
import json
import fnmatch


# %%
def main():

    # For testing
    proj_dir = "/home/data/madlab/McMakin_EMUR01"
    scratch_dir = "/scratch/madlab/emu_unc"

    # set up
    kmean_dir = os.path.join(proj_dir, "derivatives/kmeans")
    work_dir = os.path.join(scratch_dir, "derivatives/kmeans")
    with open(
        os.path.join(proj_dir, "code/amygdala_kmeans/initial_clustering.json")
    ) as jf:
        label_defs = json.load(jf)

    subj_list = [x for x in label_defs.keys()]


if __name__ == "__main__":
    main()
