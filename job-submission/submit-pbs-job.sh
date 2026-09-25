#!/bin/bash

# ============================================================================ #
# DEFAULT VALUES                                                               #
# ============================================================================ #

PRJ_DIR="/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project"
PBS_TEMPLATE="$PRJ_DIR/job-submission/pbs_submit_job.sh"
DEFAULT_ENV="scanpy_env"

DEFAULT_MEM="300g"
DEFAULT_CORES="10"
DEFAULT_TIME="03:00:00"

show_help() {
    echo "Usage: $(basename "$0") -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]"
    echo ""
    echo "Required Parameters:"
    echo "  -e  Experiment folder name (EXP) [Accepts partial names like '01' or 'practice']"
    echo "  -s  Python script filename to run (SCR) [should end with .py]"
    echo ""
    echo "Optional Parameters:"
    echo "  -i  Launch the interactive setup wizard"
    echo "  -j  Custom job name"
    echo "  -n  Conda environment name [Default: $DEFAULT_ENV]"
    echo "  -p  Directory to project [Default: $PRJ_DIR]"
    echo ""
    echo "Optional Cluster Resources:"
    echo "  -m  Memory allocation [Default: $DEFAULT_MEM]"
    echo "  -c  Number of CPU cores [Default: $DEFAULT_CORES]"
    echo "  -t  Walltime limit (HH:MM:SS) [Default: $DEFAULT_TIME]"
    exit 0
}

# ============================================================================ #
# GET VALUES                                                                   #
# ============================================================================ #

exp_val=""
script_val=""
job_name=""
env_val=""
prj_dir=""
mem_val=""
cores_val=""
time_val=""
interactive_mode=false

# Parse options
while getopts "hie:s:n:j:p:m:c:t:" opt; do
    case "$opt" in
        h) show_help ;;
        i) interactive_mode=true ;;
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
# INTERACTIVE MODE                                                             #
# ============================================================================ #

if [ "$interactive_mode" = true ]; then
    echo "⚙️  Entering Interactive Job Submission Setup..."
    echo "   (Pre-filled flags are skipped. Press [Enter] to accept defaults)"
    echo "=========================================================="
    
    while [[ -z "$exp_val" ]]; do
        read -p "📁 Enter Experiment name (or partial match): " exp_val
        [[ -z "$exp_val" ]] && echo "❌ Error: Experiment name cannot be blank."
    done

    while [[ -z "$script_val" ]]; do
        read -p "📄 Enter Python script filename: " script_val
        [[ -z "$script_val" ]] && echo "❌ Error: Script name cannot be blank."
    done

    if [[ -z "$job_name" ]]; then
        read -p "⚙️  Enter Custom Job Name [defaults to script name]: " input_job
        job_name="$input_job"
    fi

    if [[ -z "$env_val" ]]; then
        read -p "⚙️  Enter Conda Env [$DEFAULT_ENV]: " input_env
        env_val="${input_env:-$DEFAULT_ENV}"
    fi

    if [[ -z "$prj_dir" ]]; then
        read -p "📁 Enter Project Directory [$PRJ_DIR]: " input_prj
        prj_dir="${input_prj:-$PRJ_DIR}"
    fi

    if [[ -z "$cores_val" ]]; then
        read -p "📊 Enter CPU Cores [$DEFAULT_CORES]: " input_cores
        cores_val="${input_cores:-$DEFAULT_CORES}"
    fi

    if [[ -z "$mem_val" ]]; then
        read -p "📊 Enter Memory [$DEFAULT_MEM]: " input_mem
        mem_val="${input_mem:-$DEFAULT_MEM}"
    fi

    if [[ -z "$time_val" ]]; then
        read -p "⏳ Enter Walltime [$DEFAULT_TIME]: " input_time
        time_val="${input_time:-$DEFAULT_TIME}"
    fi

    echo ""
    echo "submit-job -m $mem_val -c $cores_val -t $time_val -n $env_val -e $exp_val -j $job_name -s $script_val -p $prj_dir"
    echo "=========================================================="
fi

# ============================================================================ #
# FILL LATE DEFAULTS FOR PURE FLAG MODE                                        #
# ============================================================================ #

[[ -z "$env_val" ]]   && env_val="$DEFAULT_ENV"
[[ -z "$prj_dir" ]]   && prj_dir="$PRJ_DIR"

# Resources
[[ -z "$mem_val" ]]   && mem_val="$DEFAULT_MEM"
[[ -z "$time_val" ]]  && time_val="$DEFAULT_TIME"
[[ -z "$cores_val" ]] && cores_val="$DEFAULT_CORES"

# ============================================================================ #
# VALIDATION                                                                   #
# ============================================================================ #

if [[ -z "$exp_val" || -z "$script_val" ]]; then
    echo "❌ Error: Both -e (experiment) and -s (script) are required."
    show_help
fi

if [ ! -f "$PBS_TEMPLATE" ]; then
    echo "❌ Error: PBS Template file not found at:"
    echo "   $PBS_TEMPLATE"
    echo "Aborting."
    exit 1
fi

# DYNAMIC MATCHING ----------------------------------------------------------- #

if [ -d "$prj_dir" ]; then
    # Scan for directories matching the pattern, extracting just the folder names
    matches=($(find "$prj_dir" -maxdepth 1 -type d -name "*${exp_val}*" -exec basename {} \;))
    if [ ${#matches[@]} -eq 1 ]; then
        # Exactly one match found! Automatically replace input with full name
        echo "🔍 Found unique match for '$exp_val': ${matches[0]}"
        exp_val="${matches[0]}"
    elif [ ${#matches[@]} -gt 1 ]; then
        # Multiple matching folders found
        echo "⚠️  Ambiguous input '$exp_val'. Multiple folders match:"
        for match in "${matches[@]}"; do
            echo "   - $match"
        done
        echo "Aborting. Please be more specific."
        exit 1
    fi
fi
# ---------------------------------------------------------------------------- #

EXP_DIR="${prj_dir}/${exp_val}"
if [ ! -d "$EXP_DIR" ]; then
    echo "❌ Error: Experiment directory does not exist:"
    echo "$EXP_DIR"
    echo "Aborting."
    exit 1
fi

TARGET_SCRIPT="${EXP_DIR}/scripts/${script_val}"
if [[ ! -f "$TARGET_SCRIPT" ]]; then
    echo "❌ Error: Python script not found at the expected path:"
    echo "   $TARGET_SCRIPT"
    echo "Aborting."
    exit 1
fi

if [[ "$script_val" != *.py ]]; then
    echo "⚠️ Warning: The script '$script_val' does not end in '.py'. Is this correct?"
    read -p "   Continue anyway? (y/N): " confirm
    [[ "$confirm" != [yY] ]] && exit 1
fi

# ============================================================================ #
# JOB NAME                                                                     #
# ============================================================================ #

if [[ -z "$job_name" ]]; then
    base_script=$(basename "$script_val")
    job_name="${base_script%.py}"
fi

# ============================================================================ #
# SUBMIT JOB                                                                   #
# ============================================================================ #

echo "🚀 Submitting job '$job_name'"
echo "   Resources: Cores=$cores_val | Mem=$mem_val | Time=$time_val"

job_id=$(qsub -N "$job_name" \
     -l mem="$mem_val" \
     -l nodes=1:ppn="$cores_val" \
     -l walltime="$time_val" \
     -v EXP="$exp_val",ENV="$env_val",SCR="$script_val",DIR="$prj_dir",CORES="$cores_val",MEM="$mem_val",TIME="$time_val" \
     "$PBS_TEMPLATE")

clean_id=$(echo "$job_id" | grep -oE '^[0-9]+')
echo "✅ Success! Job ID: $clean_id is now in the queue."
echo
echo "📁 Log Directory: $EXP_DIR/logs"
echo "📄 Log File: $job_name.out"

qstat -a -u "$USER"