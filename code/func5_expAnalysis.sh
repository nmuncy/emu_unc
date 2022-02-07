#!/bin/bash

#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -p IB_44C_512G
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 16000
#SBATCH --job-name ppiETAC

# load relevant modules
module load afni-20.2.06
export TMPDIR=/scratch/madlab/emu_unc/derivatives/etac

# help
function Usage {
    cat <<USAGE
    Generate and run paired t-test via 3dttest++ ETAC.

    Find subjects who have PPI output (ref func2_ppi.py), then build ETAC
    command for subjects who have both input behaviors. Assumes example
    strings/locations described in required args, and requires group
    intersection mask produced by func3_roiAnalysis.sh.

    Generates an ETAC script in project derivatives/emu_unc called etac_<ppi_seed>.sh,
    and output written to same directory.

    Required Arguments:
        -d </path/to/dir> = location of project derivatives directory, should contain both
            afni and emu_unc sub-directories.
        -p <ppi_seed> = identifying PPI seed name
            (e.g. amgL to find <data_dir>/<subj>/<sess>/func/decon_task-test_UniqueBehs_PPI-amgL_stats_REML+tlrc.HEAD)
        -s <session> = BIDS session string
        <behaviors> = remaining args are sub-brick behaviors to extract (ref 3dinfo -verb)
            note - exactly 2 must be given

    Example Usage:
        sbatch func5_expAnalysis.sh \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives \\
            -p amgL \\
            -s ses-S2 \\
            SnegLF SneuLF

USAGE
}

# receive args
while getopts ":d:p:s:h" OPT; do
    case $OPT in
    d)
        proj_dir=${OPTARG}
        if [ ! -d ${proj_dir}/afni ] || [ ! -d ${proj_dir}/emu_unc ]; then
            echo -e "\n\t ERROR: did not detect $proj_dir or required sub-directories." >&2
            Usage
            exit 1
        fi
        ;;
    p)
        ppi_seed=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    h)
        Usage
        exit 0
        ;;
    :)
        echo -e "\n\t ERROR: option '$OPTARG' missing argument." >&2
        Usage
        exit 1
        ;;
    \?)
        echo -e "\n\t ERROR: invalid option '$OPTARG'." >&2
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

# get remaning args as behaviors
shift "$((OPTIND - 1))"
beh_list=("$@")

# make sure required args have values - determine which (first) arg is empty
function emptyArg {
    case $1 in
    proj_dir)
        h_ret="-d"
        ;;
    sess)
        h_ret="-s"
        ;;
    ppi_seed)
        h_ret="-p"
        ;;
    *)
        echo -n "Unknown option."
        ;;
    esac
    echo -e "\n\n \t ERROR: Missing input parameter for \"${h_ret}\"." >&2
    Usage
    exit 1
}

for opt in proj_dir sess ppi_seed; do
    h_opt=$(eval echo \${$opt})
    if [ -z $h_opt ]; then
        emptyArg $opt
    fi
done

# check for passed behaviors
if [ ${#beh_list[@]} != 2 ]; then
    echo -e "\n\n\t ERROR: Incorrect number of behavior inputs." >&2
    Usage
    exit 1
fi

# print report
cat <<-EOF

    Checks passed, options captured:
        -d : $proj_dir
        -p : $ppi_seed
        -s : $sess
        beh : ${beh_list[@]}

EOF

# # for testing
# sess=ses-S2
# ppi_seed=amgL
# proj_dir=/home/data/madlab/McMakin_EMUR01/derivatives
# beh_list=(SnegLF SneuLF)

# set up
out_dir=${proj_dir}/emu_unc
data_dir=${proj_dir}/afni
group_mask=${out_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_desc-grpIntx_mask.nii.gz

# find subjs with PPI output
subj_list_all=($(ls $data_dir | grep "sub-*"))
subj_list=()
for subj in ${subj_list_all[@]}; do
    check_file=${data_dir}/${subj}/${sess}/func/decon_task-test_UniqueBehs_PPI-${ppi_seed}_stats_REML+tlrc.HEAD
    if [ -f $check_file ]; then
        subj_list+=($subj)
    fi
done
echo -e "\nSubject list:\n\t${subj_list[@]}\n"

# build etac structure
echo -e "Building ETAC command ...\n"
out_str=FINAL_PPI-${ppi_seed}_${beh_list[0]}-${beh_list[1]}
etacCmd=(3dttest++
    -paired
    -mask $group_mask
    -prefix $out_str
    -prefix_clustsim ${out_str}_clustsim
    -ETAC
    -ETAC_opt NN=2:sid=2:hpow=0:pthr=0.01,0.005,0.002,0.001:name=etac
)

# build set A, B
setA=()
setB=()
for subj in ${subj_list[@]}; do
    ppi_file=${data_dir}/${subj}/${sess}/func/decon_task-test_UniqueBehs_PPI-${ppi_seed}_stats_REML+tlrc
    beh_arr=()
    for beh in ${beh_list[@]}; do
        h_brick=$(3dinfo -label2index "${beh}#0_Coef" $ppi_file)
        beh_arr+=($h_brick)
    done
    if [ ${#beh_arr[@]} != 2 ]; then
        echo -e "\n\t Missing behavior for $subj, skipping ..."
        continue
    fi
    setA+=($subj "${ppi_file}[${beh_arr[0]}]")
    setB+=($subj "${ppi_file}[${beh_arr[1]}]")
done

# update etac command with sets A/B
etacCmd+=(-setA ${beh_list[0]} ${setA[@]})
etacCmd+=(-setB ${beh_list[1]} ${setB[@]})

# print etac command for review, run
echo "${etacCmd[@]}" >${out_dir}/etac_${ppi_seed}.sh
echo -e "Starting ETAC ..."
cd $out_dir
"${etacCmd[@]}"
