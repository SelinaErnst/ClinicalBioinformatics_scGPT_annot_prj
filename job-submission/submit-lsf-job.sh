#!/bin/bash

# ============================================================================ #
# DEFAULT VALUES                                                               #
# ============================================================================ #

PRJ_DIR="/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project"
DEFAULT_ENV="scanpy_env"

DEFAULT_MEM="150GB"
DEFAULT_CORES="10"
DEFAULT_TIME="02:00"

show_help() {
    echo "Usage: $(basename "$0") -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]"
    echo ""
    echo "Required Parameters:"
    echo "  -e  Experiment folder name (EXP) [Accepts partial names like '01' or 'practice']"
    echo "  -s  Python script filename to run (SCRIPT) [should end with .py]"
    echo ""
    echo "Optional Parameters:"
    echo "  -i  Launch the interactive setup wizard"
    echo "  -j  Custom job name"
    echo "  -n  Conda environment name (ENV) [Default: $DEFAULT_ENV]"
    echo "  -p  Directory to project (PRJ_DIR) [Default: $PRJ_DIR]"
    echo ""
    echo "Optional Cluster Resources:"
    echo "  -m  Memory allocation [Default: $DEFAULT_MEM]"
    echo "  -c  Number of CPU cores [Default: $DEFAULT_CORES]"
    echo "  -t  Walltime limit (HH:MM) [Default: $DEFAULT_TIME]"
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
    echo "submit-job -e $exp_val -s $script_val -n $env_val -j $job_name -m $mem_val -c $cores_val -t $time_val -p $prj_dir"
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
JOB_DIR="${EXP_DIR}/jobs"
LOG_DIR="${JOB_DIR}/logs"
SCR_DIR="${EXP_DIR}/scripts"

if [ ! -d "$EXP_DIR" ]; then
    echo "❌ Error: Experiment directory does not exist:"
    echo "   $EXP_DIR"
    exit 1
fi

TARGET_SCRIPT="${SCR_DIR}/${script_val}"
if [[ ! -f "$TARGET_SCRIPT" ]]; then
    echo "❌ Error: Python script not found at the expected path:"
    echo "   $TARGET_SCRIPT"
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

mkdir -p "$LOG_DIR"
SUBMIT_SCRIPT="$JOB_DIR/lsf_submit_$job_name.sh"

cat << EOF > "$SUBMIT_SCRIPT"
#!/bin/bash
#BSUB -J $job_name
#BSUB -q long
#BSUB -n $cores_val
#BSUB -R "rusage[mem=$mem_val] span[hosts=1]"
#BSUB -W $time_val
#BSUB -o $LOG_DIR/${job_name}.out
#BSUB -e $LOG_DIR/${job_name}.err
START_TIME=\$(date +%s)
source ~/.bashrc
cd $SCR_DIR
conda activate $env_val
python $script_val
END_TIME=\$(date +%s)
DURATION=\$((\$END_TIME - \$START_TIME))
printf "\n"; printf "%.0s#" {1..58}; printf "\n"
printf "⏱️ Total Runtime: %02d:%02d:%02d\n" \$((\$DURATION/3600)) \$((\$DURATION%3600/60)) \$((\$DURATION%60))
printf "%.0s#" {1..58}; printf "\n\n"
EOF

# ============================================================================ #
# FINAL PRE-FLIGHT VERIFICATION PROMPT                                         #
# ============================================================================ #

echo ""
echo "📋 Ready to submit job '$job_name' with configurations:"
echo "   - Cores:    $cores_val"
echo "   - Memory:   $mem_val"
echo "   - Time:     $time_val"
echo "   - Script:   $SUBMIT_SCRIPT"
echo ""

read -p "🚀 Proceed with submission to the cluster queue? (y/N): " run_confirm

if [[ "$run_confirm" != [yY] && "$run_confirm" != [yY][eE][sS] ]]; then
    echo "🛑 Submission aborted by user. The cluster queue was not contacted."
    echo "💾 The generated script has been saved for manual modification at:"
    echo "   $SUBMIT_SCRIPT"
    exit 0
fi

# ============================================================================ #
# EXECUTE SUBMISSION                                                           #
# ============================================================================ #

chmod +x $SUBMIT_SCRIPT
submit_output=$(bsub < $SUBMIT_SCRIPT)

clean_id=$(echo "$submit_output" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+')

echo "✅ Success! Job ID: $clean_id is now in the queue."
echo 
echo "📁 Log Directory: $LOG_DIR"
echo "📄 Log File: $job_name.out"


# ============================================================================ #
# KEEP SCRIPT PROMPT                                                           #
# ============================================================================ #

echo ""
# Ask user if they want to retain the file before we construct the pathways
read -p "❓ Delete the generated submission script file? (y/N): " del_confirm


# Evaluate retention decision
if [[ "$del_confirm" == [yY] || "$del_confirm" == [yY][eE][sS] ]]; then
    rm -f "$SUBMIT_SCRIPT"
else
    echo "📄 Saved Submission File: ${SUBMIT_SCRIPT}"
fi