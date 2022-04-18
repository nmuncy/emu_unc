#!/bin/bash

function Usage {
    cat <<USAGE

    Wrapper for func4_ppi.py. Checks for which subjects do not have
    PPI output, and submits N sbatch jobs for detected subjects.

    Conda environment emuR01_unc required.

    Required Arguments:
        -d <project_derivatives> = path to project derivatives location
        -f <decon_string> = prefix of decon file
        -i <mask_info> = path to a mask, or coordinates for mask to be
            constructed. Only one input argument accepted for this option.
        -n <number> = number of subjects to submit jobs
        -r <seed_name> = seed\'s name, for writing files, sub-bricks
        -s <session> = BIDS session string
        -t <task> = BIDS task string
        -w <scratch_directory> = path to scratch/working directory

    Example Usage:
        deriv_dir=/home/data/madlab/McMakin_EMUR01/derivatives/emu_unc
        $0 \\
            -w /scratch/madlab/emu_unc/derivatives/afni_ppi \\
            -d \$deriv_dir \\
            -f decon_task-study_rVal \\
            -s ses-S1 \\
            -t task-study \\
            -r amgL \\
            -i \${deriv_dir}/template/tpl-MNIPediatricAsym_cohort-5_res-2_desc-amgLClean_mask.nii.gz \\
            -n 8

        deriv_dir=/home/data/madlab/McMakin_EMUR01/derivatives/emu_unc
        $0 \\
            -w /scratch/madlab/emu_unc/derivatives/afni_ppi \\
            -d \$deriv_dir \\
            -f decon_task-test_UniqueBehs \\
            -s ses-S2 \\
            -t task-study \\
            -r LHC \\
            -i "-24 -12 -22" \\
            -n 8

USAGE
}

# capture arguments
while getopts ":d:f:i:n:r:s:t:w:h" OPT; do
    case $OPT in
    d)
        deriv_dir=${OPTARG}
        if [ ! -d $deriv_dir ]; then
            echo -e "\n\t ERROR: -d directory not found.\n" >&2
            Usage
            exit 1
        fi
        ;;
    f)
        decon_str=${OPTARG}
        ;;
    i)
        seed_info=${OPTARG}
        if [ ${#seed_info} -lt 5 ]; then
            echo -e "\n\t ERROR: argument -i <seed-info> not configured correctly." >&2
            Usage
            exit 1
        fi
        ;;
    n)
        num_subj=${OPTARG}
        if [ $num_subj -lt 1 ]; then
            echo -e "\n\t ERROR: -n value must be greater than 0.\n" >&2
            Usage
            exit 1
        fi
        ;;
    r)
        seed_name=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    t)
        task=${OPTARG}
        ;;
    w)
        scratch_dir=${OPTARG}
        ;;
    h)
        Usage
        exit 0
        ;;
    \?)
        echo "\n\t Error: invalid option -${OPTARG}."
        Usage
        exit 1
        ;;
    :)
        echo "\n\t Error: -${OPTARG} requires an argument."
        Usage
        exit 1
        ;;
    esac
done

# print help if no arg
if [ $OPTIND == 1 ]; then
    Usage
    exit 0
fi

# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
    deriv_dir)
        h_ret="-d"
        ;;
    decon_str)
        h_ret="-f"
        ;;
    seed_info)
        h_ret="-i"
        ;;
    num_subj)
        h_ret="-n"
        ;;
    seed_name)
        h_ret="-r"
        ;;
    sess)
        h_ret="-s"
        ;;
    task)
        h_ret="-t"
        ;;
    scratch_dir)
        h_ret="-w"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in deriv_dir decon_str seed_info num_subj seed_name sess task scratch_dir; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for conda env
which python | grep "emuR01_unc" >/dev/null 2>&1
if [ $? != 0 ]; then
    echo -e "\n\t ERROR: Please conda activate emuR01_unc_env and try again.\n" >&2
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, options captured:
        -w : $scratch_dir
        -d : $deriv_dir
        -f : $decon_str
        -s : $sess
        -t : $task
        -n : $num_subj
        -r : $seed_name
        -i : $seed_info

EOF

# find subjects missing ppi output
echo -e "Building subject lists ...\n"
subj_list=()
subj_all=($(ls $deriv_dir | grep "sub-*"))
for subj in ${subj_all[@]}; do
    ppi_file=${deriv_dir}/${subj}/${sess}/func/${decon_str}_PPI-${seed_name}_stats_REML+tlrc.HEAD
    if [ ! -f $ppi_file ]; then
        subj_list+=($subj)
    fi
done

# patch - remove subjs w/missing data, redeclare array
# problem_list=(sub-4{011,020,021,055,056,063,090,172,197})  # for test
problem_list=(sub-4001) # for study
ind_problem=()
for ind in ${!subj_list[@]}; do
    for prob in ${problem_list[@]}; do
        if [[ "${subj_list[$ind]}" == "${prob}" ]]; then
            ind_problem+=($ind)
        fi
    done
done

for ind in "${ind_problem[@]}"; do
    unset "subj_list[$ind]"
done
declare -a subj_list=(${subj_list[@]})

# submit N jobs
time=$(date '+%Y-%m-%d_%H:%M')
out_dir=${scratch_dir}/slurm_out/ppi_${time}
mkdir -p $out_dir

# reference subj timing dir
parent_dir="$(dirname $(pwd))"

cat <<-EOF
    Submitting sbatch job

        func3_ppi.py \\
            -p <subj> \\
            -s $sess \\
            -t $task \\
            -d $decon_str \\
            -r $seed_name \\
            -i "$seed_info" \\
            --timing-dir ${parent_dir}/data/timing_files/<subj>/$sess

    with the following subjects:

        ${subj_list[@]:0:$num_subj}

EOF

c=0
while [ $c -lt $num_subj ]; do

    subj=${subj_list[$c]}
    subj_time=${parent_dir}/data/timing_files/${subj}/$sess
    sbatch \
        --job-name=p${subj#*-} \
        --output=${out_dir}/${subj}.out \
        --mem-per-cpu=4000 \
        --partition=IB_44C_512G \
        --account=iacc_madlab \
        --qos=pq_madlab \
        func3_ppi.py \
        -p $subj \
        -s $sess \
        -t $task \
        -d $decon_str \
        -r $seed_name \
        -i "$seed_info" \
        --timing-dir $subj_time

    sleep 1
    let c+=1
done
