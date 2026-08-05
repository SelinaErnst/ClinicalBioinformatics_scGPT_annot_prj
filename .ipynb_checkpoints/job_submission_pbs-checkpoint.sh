#!/bin/bash
#PBS -l mem=250g
#PBS -l nodes=1:ppn=10
#PBS -l walltime=02:00:00
#PBS -o /dev/null
#PBS -e /dev/null

START_TIME=$(date "+%Y-%m-%d %H:%M:%S")
START_SECONDS=$SECONDS

# ============================================================================ #
# PARAMETERS                                                                   #
# ============================================================================ #

EXP="${EXP:?Error: You must provide the experiment using -v EXP=experimentXY}"
ENV="${ENV:?Error: You must provide the conda env using -v EXP=env_name}"
PY="${SCRIPT:?Error: You must provide a script name using -v SCRIPT=your_script.py}"

# ============================================================================ #
# DIRECTORIES                                                                  #
# ============================================================================ #

PRJ_DIR="/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project"
EXP_DIR="${PRJ_DIR}/${EXP}"
LOG_DIR="${EXP_DIR}/logs"
SCR_DIR="${EXP_DIR}/scripts"

FILE_ERR="${LOG_DIR}/${PBS_JOBNAME}.err"
FILE_LOG="${LOG_DIR}/${PBS_JOBNAME}.out"

# ============================================================================ #
# HEADER                                                                       #
# ============================================================================ #

{
    echo "##########################################################"
    echo "                    PBS JOB METADATA                      "
    echo "##########################################################"
    echo "Job Name:        ${PBS_JOBNAME}"
    echo "Job ID:          ${PBS_JOBID}"
    echo "Experiment:      ${EXP}"
    echo "Conda Env:       ${ENV}"
    echo "Python Script:   ${PY}"
    echo "Project Dir:     ${PRJ_DIR}"
    echo "Started At:      ${START_TIME}"
    echo "Allocated Cores: ${PBS_NUM_PPN:-10}"
    echo "=========================================================="
    echo ""
} >> "$FILE_LOG"

# ============================================================================ #
# ACTION                                                                       #
# ============================================================================ #

cd $SCR_DIR
source ~/.bashrc

module load Miniconda3/4.5.11
conda activate "$ENV"

mkdir "$LOG_DIR"
python -u "$PY" >> "$FILE_LOG" 2> "$FILE_ERR"

# ============================================================================ #
# SUMMARY                                                                      #
# ============================================================================ #

EXIT_CODE=$?
END_TIME=$(date "+%Y-%m-%d %H:%M:%S")
DURATION=$((SECONDS - START_SECONDS))

# Convert duration seconds into hours, minutes, seconds
H=$((DURATION / 3600))
M=$(((DURATION % 3600) / 60))
S=$((DURATION % 60))

{
    echo ""
    echo "=========================================================="
    echo "                   EXECUTION SUMMARY                      "
    echo "=========================================================="
    echo "Finished At:     ${END_TIME}"
    echo "Total Runtime:   ${H}h ${M}m ${S}s"
    echo "Exit Status:     ${EXIT_CODE} ($( [ $EXIT_CODE -eq 0 ] && echo "SUCCESS" || echo "FAILED" ))"
    echo "##########################################################"
} >> "$FILE_LOG"

# ============================================================================ #
# CLEANUP                                                                      #
# ============================================================================ #

if [ ! -s "$FILE_ERR" ]; then
    rm "$FILE_ERR"
fi