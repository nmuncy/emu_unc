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

    Find subjects who have PPI output (ref func3_ppi.py), then build ETAC
    command for subjects who have both input behaviors. Assumes example
    strings/locations described in required args, and requires group
    intersection mask produced by func4_roiAnalysis.sh.

    Generates an ETAC script in project derivatives/emu_unc called etac_<ppi_seed>.sh,
    and output written to same directory.

    Required Arguments:
        -c <path/to/code/dir> = path to clone of https://github.com/nmuncy/emu_unc,
            should contain "code" and "data" sub-directories.
        -d </path/to/dir> = location of project derivatives directory
        -n <decon_name> = identifying deconvolution name
        -p <ppi_seed> = identifying PPI seed name
            (e.g. amgL to find <deriv_dir>/<subj>/<sess>/func/decon_task-test_UniqueBehs_PPI-amgL_stats_REML+tlrc.HEAD)
        -s <session> = BIDS session string
        -t <task> = BIDS task string
        <behaviors> = remaining args are sub-brick behaviors to extract (ref 3dinfo -verb)
            note - exactly 2 must be given

    Example Usage:
        code_dir="\$(dirname "\$(pwd)")"
        sbatch func7_expAnalysis.sh \\
            -c \$code_dir \\
            -d /home/data/madlab/McMakin_EMUR01/derivatives/emu_unc \\
            -n precTest \\
            -p amgL \\
            -s ses-S1 \\
            -t task-study \\
            SPnegLF SPneuLF

USAGE
}

# receive args
while getopts ":c:d:n:p:s:t:h" OPT; do
    case $OPT in
    c)
        code_dir=${OPTARG}
        if [ ! -d ${code_dir} ] || [ ! -d ${code_dir}/data ]; then
            echo -e "\n\t ERROR: did not detect -c $code_dir or required sub-directory \"data\"." >&2
            Usage
            exit 1
        fi
        ;;
    d)
        deriv_dir=${OPTARG}
        if [ ! -d ${deriv_dir} ]; then
            echo -e "\n\t ERROR: did not detect -d $deriv_dir." >&2
            Usage
            exit 1
        fi
        ;;
    n)
        decon_name=${OPTARG}
        ;;
    p)
        ppi_seed=${OPTARG}
        ;;
    s)
        sess=${OPTARG}
        ;;
    t)
        task=${OPTARG}
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
    code_dir)
        h_ret="-c"
        ;;
    deriv_dir)
        h_ret="-d"
        ;;
    decon_name)
        h_ret="-n"
        ;;
    sess)
        h_ret="-s"
        ;;
    task)
        h_ret="-t"
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

for opt in code_dir deriv_dir decon_name sess task ppi_seed; do
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
        -c : $code_dir
        -d : $deriv_dir
        -n : $decon_name
        -p : $ppi_seed
        -s : $sess
        -t : $task
        beh : ${beh_list[@]}

EOF

# set up
data_dir=${code_dir}/data
template_dir=${deriv_dir}/template
analysis_dir=${deriv_dir}/analyses
group_mask=${template_dir}/tpl-MNIPediatricAsym_cohort-5_res-2_${sess}_${task}_desc-grpIntx_mask.nii.gz

# find subjs with PPI output - start subject list with NSlacc sufficient trial list, from func5_mkdf.R
# subj_list_all=($(ls $deriv_dir | grep "sub-*"))
ref_file=${data_dir}/df_${sess}_${task}_${ppi_seed}-NSlacc_subjs-sufficient-trials.txt
subj_list_all=($(cat $ref_file | awk '{print $1}'))
subj_list=()
for subj in ${subj_list_all[@]}; do
    check_file=${deriv_dir}/${subj}/${sess}/func/decon_${task}_${decon_name}_PPI-${ppi_seed}_stats_REML+tlrc.HEAD
    if [ -f $check_file ]; then
        subj_list+=($subj)
    fi
done
echo -e "\nSubject list:\n\t${subj_list[@]}\n"

# build etac structure
echo -e "Building ETAC command ...\n"
out_str=FINAL_${sess}_${task}_PPI-${ppi_seed}_${beh_list[0]}-${beh_list[1]}
etacCmd=(3dttest++
    -paired
    -mask $group_mask
    -prefix $out_str
    -prefix_clustsim ${out_str}_clustsim
    -ETAC
    -ETAC_blur 4 8
    -ETAC_opt NN=2:sid=2:hpow=0:pthr=0.01,0.005,0.002,0.001:name=etac
)

# build set A, B
setA=()
setB=()
for subj in ${subj_list[@]}; do
    ppi_file=${deriv_dir}/${subj}/${sess}/func/decon_${task}_${decon_name}_PPI-${ppi_seed}_stats_REML+tlrc
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
echo "${etacCmd[@]}" >${analysis_dir}/etac_${sess}_${task}_${ppi_seed}.sh
echo -e "Starting ETAC ..."
cd $analysis_dir
"${etacCmd[@]}"
