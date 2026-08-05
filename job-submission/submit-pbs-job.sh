#!/bin/bash

# ============================================================================ #
# DEFAULT VALUES                                                               #
# ============================================================================ #

PRJ_DIR="/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project"
PBS_TEMPLATE="$PRJ_DIR/job-submission/job_submission_pbs.sh"
DEFAULT_ENV="scanpy_env"

DEFAULT_MEM="250g"
DEFAULT_CORES="10"
DEFAULT_TIME="02:00:00"

show_help() {
    echo "Usage: $(basename "$0") -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]"
    echo ""
    echo "Required Parameters:"
    echo "  -e  Experiment folder name (EXP)"
    echo "  -s  Python script filename to run (SCRIPT)"
    echo ""
    echo "Optional Parameters:"
    echo "  -j  Custom PBS job name"
    echo "  -n  Conda environment name (ENV) [Default: $DEFAULT_ENV]"
    echo "  -p  Directory to project (PRJ_DIR) [Default: $PRJ_DIR]"
    echo ""
    echo "Optional Cluster Resources:"
    echo "  -m  Memory allocation [Default: $DEFAULT_MEM]"
    echo "  -c  Number of CPU cores (ppn) [Default: $DEFAULT_CORES]"
    echo "  -t  Walltime limit (HH:MM:SS) [Default: $DEFAULT_TIME]"
    exit 0
}

# ============================================================================ #
# GET VALUES                                                                   #
# ============================================================================ #

exp_val=""
script_val=""
job_name=""
env_val="$DEFAULT_ENV"
prj_dir="$PRJ_DIR"

# Resource variable initialization
mem_val="$DEFAULT_MEM"
cores_val="$DEFAULT_CORES"
time_val="$DEFAULT_TIME"

# Parse options

while getopts "he:s:n:j:p:m:c:t:" opt; do
    case "$opt" in
        h) show_help ;;
        e) exp_val=$OPTARG ;;
        s) script_val=$OPTARG ;;
        n) env_val=$OPTARG ;;
        j) job_name=$OPTARG ;;
        p) prj_dir=$OPTARG ;;
        m) mem_val=$OPTARG ;;  
        c) cores_val=$OPTARG ;;
        t) time_val=$OPTARG ;; 
        *) show_help ;;
    esac
done

# ============================================================================ #
# VALIDATION                                                                   #
# ============================================================================ #

if [[ -z "$exp_val" || -z "$script_val" ]]; then
    echo "Error: Both -e (experiment) and -s (script) are required."
    show_help
fi

if [[ "$script_val" != *.py ]]; then
    echo "Warning: The script '$script_val' does not end in '.py'. Is this correct?"
    read -p "Continue anyway? (y/N): " confirm
    [[ "$confirm" != [yY] ]] && exit 1
fi

EXP_DIR="${prj_dir}/${exp_val}"
if [ ! -d "$EXP_DIR" ]; then
    echo "Warning: Experiment directory does not exist:"
    echo "$EXP_DIR"
    read -p "Continue anyway? (y/N): " confirm
    [[ "$confirm" != [yY] ]] && exit 1
fi

# ============================================================================ #
# JOB NAME                                                                     #
# ============================================================================ #
# If no job name is specified, strip '.py' from the script filename
if [[ -z "$job_name" ]]; then
    # Removes the path if a user types folder/script.py, then strips .py
    base_script=$(basename "$script_val")
    job_name="${base_script%.py}"
fi

# ============================================================================ #
# SUBMIT JOB                                                                   #
# ============================================================================ #

echo "🚀 Submitting job '$job_name'..."
echo "   Resources: Cores=$cores_val | Mem=$mem_val | Time=$time_val"

job_id=$(qsub -N "$job_name" \
     -l mem="$mem_val" \
     -l nodes=1:ppn="$cores_val" \
     -l walltime="$time_val" \
     -v EXP="$exp_val",ENV="$env_val",SCRIPT="$script_val",PRJ_DIR="$prj_dir" \
     "$PBS_TEMPLATE")

clean_id=$(echo "$job_id" | grep -oE '^[0-9]+')
echo "✅ Success! Job ID: $clean_id is now in the queue."

qstat -a -u "$USER"