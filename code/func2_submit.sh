#!/bin/bash

function Usage {
    cat <<USAGE
    Foo

    Example Usage:
        $0 \\
            -w /scratch/madlab/emu_unc/derivatives/afni_ppi \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives/afni \\
            -f decon_task-test_UniqueBehs_PPI-LHC \\
            -s ses-S2 \\
            -r blaL \\
            -i /home/data/madlab/McMakin_EMUR01/derivatives/emu_unc/tpl-MNIPediatricAsym_cohort-5_res-2_desc-blaL_mask.nii.gz \\
            -n 8

        $0 \\
            -w /scratch/madlab/emu_unc/derivatives/afni_ppi \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives/afni \\
            -f decon_task-test_UniqueBehs_PPI-LHC \\
            -s ses-S2 \\
            -r LHC \\
            -i "-24 -12 -22" \\
            -n 8

USAGE
}

# capture arguments
while getopts ":d:f:i:n:r:s:w:" OPT; do
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
        ppi_str=${OPTARG}
        ;;
    i)
        seed_info=${OPTARG}
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
    w)
        scratch_dir=${OPTARG}
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

# check args
if [ -z $scratch_dir ]; then
    echo -e "\n\t ERROR: please specify -w directory.\n" >&2
    Usage
    exit 1
fi

if [ -z $ppi_str ]; then
    echo -e "\n\t ERROR: please specify -f PPI file name.\n" >&2
    Usage
    exit 1
fi

if [ -z $sess ]; then
    echo -e "\n\t ERROR: please specify -s session.\n" >&2
    Usage
    exit 1
fi

if [ -z $seed_name ]; then
    echo -e "\n\t ERROR: please specify -r <seed name>.\n" >&2
    Usage
    exit 1
fi

if [ ${#seed_info} -lt 8 ]; then
    echo -e "\n\t ERROR: argument -i <seed-info> not configured correctly.\n" >&2
    Usage
    exit 1
fi

# check for conda env
which python | grep "emuR01_unc" >/dev/null 2>&1
if [ $? != 0 ]; then
    echo -e "\n\t ERROR: Please conda activate emuR01_unc_env and try again.\n" >&2
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, captured the options:
        -w : $scratch_dir
        -d : $deriv_dir
        -f : $ppi_str
        -s : $sess
        -n : $num_subj
        -r : $seed_name
        -i : $seed_info

EOF

# find subjects missing ppi output
subj_list=()
subj_all=($(ls $deriv_dir | grep "sub-*"))
for subj in ${subj_all[@]}; do
    ppi_file=${deriv_dir}/${subj}/${sess}/func/${ppi_str}_stats_REML+tlrc.HEAD
    if [ ! -f $ppi_file ]; then
        subj_list+=($subj)
    fi
done

# submit N jobs
time=$(date '+%Y-%m-%d_%H:%M')
out_dir=${scratch_dir}/slurm_out/ppi_${time}
mkdir -p $out_dir

echo "Submitting jobs for:"
echo -e "\t${subj_list[@]:0:$num_subj}\n"

d_arg=${ppi_str%_*}
c=0
while [ $c -lt $num_subj ]; do

    subj=${subj_list[$c]}
    sbatch \
        --job-name=p${subj#*-} \
        --output=${out_dir}/${subj}.out \
        --mem-per-cpu=4000 \
        --partition=IB_44C_512G \
        --account=iacc_madlab \
        --qos=pq_madlab \
        func2_ppi.py \
        -s $subj \
        -d $d_arg \
        -r $seed_name \
        -i "$seed_info"

    sleep 1
    let c+=1
done
