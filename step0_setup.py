import os
import shutil
import fnmatch
import json


def main():

    # set paths for source data
    code_dir = "/home/nmuncy/compute/emu_unc"
    source_dir = "/home/data/madlab/McMakin_EMUR01"
    source_dset = os.path.join(source_dir, "dset")
    source_deriv_dwi = os.path.join(source_dir, "derivatives/dwi_preproc")

    # set paths for data destination, make proj dir
    proj_dir = "/scratch/madlab/emu_UNC"
    proj_dset = os.path.join(proj_dir, "dset")
    proj_deriv = os.path.join(proj_dir, "derivatives")
    for h_dir in [proj_dset, proj_deriv]:
        if not os.path.exists(h_dir):
            os.makedirs(h_dir)

    # get bids json files
    for h_file in [
        "dataset_description.json",
        "participants.tsv",
        "task-test_bold.json",
    ]:
        shutil.copyfile(
            os.path.join(source_dset, h_file), os.path.join(proj_dset, h_file)
        )

    # determine relevant subjects
    subj_list = [x for x in os.listdir(source_deriv_dwi) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    # specify session
    sess = "ses-S2"

    # work for e/subject
    for subj in subj_list:

        print(f"\n Copying data for {subj} ...")

        # set up subj directories
        subj_dwi = os.path.join(proj_deriv, "dwi_preproc", subj, sess, "dwi")
        subj_anat = os.path.join(proj_dset, subj, sess, "anat")
        subj_func = os.path.join(proj_dset, subj, sess, "func")
        subj_fmap = os.path.join(proj_dset, subj, sess, "fmap")
        for h_dir in [subj_dwi, subj_anat, subj_func, subj_fmap]:
            if not os.path.exists(h_dir):
                os.makedirs(h_dir)

        # start copy dict
        copy_dict = {}
        for hold in ["dwi", "anat", "func", "events", "fmap"]:
            copy_dict[hold] = {"input": [], "output": []}

        """ get dwi - bval exists in different location"""
        source_dwi = os.path.join(source_deriv_dwi, subj, sess, "dwi")

        for h_suff in ["nii.gz", "bvec"]:
            copy_dict["dwi"]["input"].append(
                os.path.join(
                    source_dwi, f"{subj}_{sess}_run-1_desc-eddyCorrected_dwi.{h_suff}"
                )
            )
        copy_dict["dwi"]["input"].append(
            os.path.join(
                source_dset, subj, sess, "dwi", f"{subj}_{sess}_run-1_dwi.bval"
            )
        )

        # set output
        for h_suff in ["nii.gz", "bvec", "bval"]:
            copy_dict["dwi"]["output"].append(
                os.path.join(subj_dwi, f"{subj}_{sess}_dwi.{h_suff}")
            )

        # avoid repeating work
        if os.path.exists(copy_dict["dwi"]["output"][0]):
            continue

        """ get t1 - gathered in S1 so rename"""
        source_anat = os.path.join(source_dset, subj, "ses-S1", "anat")
        for h_suff in ["nii.gz", "json"]:
            copy_dict["anat"]["input"].append(
                os.path.join(source_anat, f"{subj}_ses-S1_run-2_T1w.{h_suff}")
            )
            copy_dict["anat"]["output"].append(
                os.path.join(subj_anat, f"{subj}_{sess}_T1w.{h_suff}")
            )

        # get func - epi and json only
        source_func = os.path.join(source_dset, subj, sess, "func")
        func_list = [
            x
            for x in os.listdir(source_func)
            if fnmatch.fnmatch(x, "*task-test*") and "events" not in x
        ]
        func_list.sort()
        for h_file in func_list:
            copy_dict["func"]["input"].append(os.path.join(source_func, h_file))
            copy_dict["func"]["output"].append(os.path.join(subj_func, h_file))

        # get new events, rename
        source_event = os.path.join(
            code_dir, "beh_data/pre_covid/task_files", subj, sess
        )
        event_list = [x for x in os.listdir(source_event)]
        for h_file in event_list:
            copy_dict["events"]["input"].append(os.path.join(source_event, h_file))
            copy_dict["events"]["output"].append(os.path.join(subj_func, h_file))

        # get relevant fmaps - only task, not rest or dwi
        source_fmap = os.path.join(source_dset, subj, sess, "fmap")
        json_list = [
            os.path.join(source_fmap, x)
            for x in os.listdir(source_fmap)
            if fnmatch.fnmatch(x, "*acq-func*json")
        ]
        # fmap_list = []
        for h_json in json_list:
            with open(h_json, "r") as jf:
                json_dict = json.loads(jf.read())
                res = [val for key, val in json_dict.items() if "task-test" in str(val)]
                if len(res) != 0:
                    copy_dict["fmap"]["input"].append(h_json)
                    copy_dict["fmap"]["input"].append(
                        f"""{h_json.split(".")[0]}.nii.gz"""
                    )

                    copy_dict["fmap"]["output"].append(
                        os.path.join(subj_fmap, h_json.split("/")[-1])
                    )
                    copy_dict["fmap"]["output"].append(
                        os.path.join(
                            subj_fmap,
                            f"""{h_json.split("/")[-1].split(".")[0]}.nii.gz""",
                        )
                    )

        # verify data
        break_status = False
        for file_type in copy_dict:
            num_items = len(copy_dict[file_type]["input"])
            for h_ind in range(0, num_items):
                h_input = copy_dict[file_type]["input"][h_ind]
                if not os.path.exists(h_input):
                    break_status = True
                    break
            if break_status:
                break

        # remove subject if break status
        if break_status:
            print(f"\t Missing data for {subj}, removing ...")
            shutil.rmtree(os.path.join(proj_dset, subj))
            shutil.rmtree(os.path.join(proj_deriv, "dwi_preproc", subj))
            continue

        # copy data
        for file_type in copy_dict:
            print(f"\t Copying {file_type} data ...")
            num_items = len(copy_dict[file_type]["input"])
            for h_ind in range(0, num_items):
                h_input = copy_dict[file_type]["input"][h_ind]
                h_output = copy_dict[file_type]["output"][h_ind]
                shutil.copyfile(h_input, h_output)


if __name__ == "__main__":
    main()
