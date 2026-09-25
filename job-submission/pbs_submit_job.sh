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

EXP="${EXP:?Error: You must provide the experiment using -v EXP=experiment}"
SCR="${SCR:?Error: You must provide a script name using -v SCR=your_script.py}"

[ -z "$ENV" ] && ENV="scanpy_env"
[ -z "$DIR" ] && DIR="/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project"

[ -z "$CORES" ] && CORES=$(grep -E '^#PBS -l nodes=' "$0" | grep -oE 'ppn=[0-9]+' | cut -d= -f2)
[ -z "$MEM" ]   && MEM=$(grep -E '^#PBS -l mem=' "$0" | cut -d= -f2)
[ -z "$TIME" ]  && TIME=$(grep -E '^#PBS -l walltime=' "$0" | cut -d= -f2)

# ============================================================================ #
# DIRECTORIES                                                                  #
# ============================================================================ #

EXP_DIR="${DIR}/${EXP}"
JOB_DIR="${EXP_DIR}/jobs"
LOG_DIR="${JOB_DIR}/logs"
SCR_DIR="${EXP_DIR}/scripts"

FILE_ERR="${LOG_DIR}/${PBS_JOBNAME}.err"
FILE_LOG="${LOG_DIR}/${PBS_JOBNAME}.out"

# ============================================================================ #
# HEADER                                                                       #
# ============================================================================ #

mkdir -p "$LOG_DIR"
{
    echo "##########################################################"
    echo "                    PBS JOB METADATA                      "
    echo "##########################################################"
    echo "Job Name:        ${PBS_JOBNAME}"
    echo "Job ID:          ${PBS_JOBID}"
    echo "Python Script:   ${SCR}"
    echo "Conda Env:       ${ENV}"
    echo "Project Dir:     ${DIR}"
    echo "Experiment:      ${EXP}"
    echo "Started At:      ${START_TIME}"
    echo "Allocated Cores: ${CORES}"
    echo "Requested Mem:   ${MEM}"
    echo "Requested Time:  ${TIME}"
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

python -u "$SCR" >> "$FILE_LOG" 2> "$FILE_ERR"

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