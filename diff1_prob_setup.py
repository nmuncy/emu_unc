import os
import toml
import json
import sys


"""
Notes:

This script will edit the BIDs json file dataset_description
    with values needed by AFQ. It will then edit the toml
    file config with project-specific values.
"""


def main():

    # receive arguments
    bids_dir = str(sys.argv[1])
    deriv_dir = str(sys.argv[2])

    # append data_description with new field
    json_file = os.path.join(deriv_dir, "dataset_description.json")
    new_field = {"PipelineDescription": {"Name": "dwi_preproc"}}

    with open(json_file) as jf:
        json_content = json.load(jf)

    if "PipelineDescription" not in json_content.keys():
        json_content.update(new_field)

    with open(json_file, "w") as jf:
        json.dump(json_content, jf)

    # edit config.toml
    toml_file = "config_prob.toml"
    toml_dict = toml.load(toml_file)
    toml_dict["TRACTOGRAPHY"]["directions"] = "prob"
    toml_dict["BIDS"]["bids_path"] = bids_dir
    toml_dict["files"]["dmriprep_folder"] = deriv_dir
    toml_dict["SEGMENTATION"]["seg_algo"] = "reco80"

    tf = open(toml_file, "w")
    toml.dump(toml_dict, tf)
    tf.close()


if __name__ == "__main__":
    main()
