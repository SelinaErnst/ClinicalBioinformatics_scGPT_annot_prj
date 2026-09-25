
# job submission

I created bash scrits to help me write job-submission files. I use the alias `submit-job` to access these scripts.  
> ⚠️ **WARNING:** Since my project folder is organized into experiments with specific sub folders, these scripts would not work in a different setting.

To make it even easier to use, it is helpful to write a line to `~/.bashrc`.  
This way, you can use the command `submit-job` to automatically run scripts from specified folders.   

```bash
# PBS:
alias submit-job='/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/job-submission/job_submission_pbs.sh'
# LSF:
alias submit-job='/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/job-submission/job_submission_lsf.sh'
```

### LSF

The LSF version of `submit-job` is more recent and had some updates. It creates an actual job submission script which can be saved if required (recommended if process needs to be rerun).  

```
❌ Error: Both -e (experiment) and -s (script) are required.
Usage: submit-lsf-job.sh -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]

Required Parameters:
  -e  Experiment folder name (EXP) [Accepts partial names like '01' or 'practice']
  -s  Python script filename to run (SCRIPT) [should end with .py]

Optional Parameters:
  -i  Launch the interactive setup wizard
  -j  Custom job name
  -n  Conda environment name (ENV) [Default: scanpy_env]
  -p  Directory to project (PRJ_DIR) [Default: /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project]

Optional Cluster Resources:
  -m  Memory allocation [Default: 150GB]
  -c  Number of CPU cores [Default: 10]
  -t  Walltime limit (HH:MM) [Default: 02:00]
```
Here is an example:

```
submit-job -m 150GB -t 01:00 -c 10 -n scanpy_env -e practice -j plot_umaps -s 04_plot_umaps.py 
```

### PBS 

The script `submit-pbs-job.sh` relies on another script `pbs_submit_job.sh` that serves as a universal job submission script. This script outputs its own log files that contain run time and other information about the job. The error file also is created but it will be deleted if after running the job it is still empty. 

> ⚠️ **WARNING:** If the script requires any module to be loaded it is recommended to use LSF instead because the universal job submission would have to include these additional lines.  
> ⚠️ **WARNING:** Do not rerun jobs unless the former one has run through. The log files would get mixed up otherwise.


```
❌ Error: Both -e (experiment) and -s (script) are required.
Usage: submit-pbs-job.sh -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]

Required Parameters:
  -e  Experiment folder name (EXP) [Accepts partial names like '01' or 'practice']
  -s  Python script filename to run (SCR) [should end with .py]

Optional Parameters:
  -i  Launch the interactive setup wizard
  -j  Custom job name
  -n  Conda environment name [Default: scanpy_env]
  -p  Directory to project [Default: /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project]

Optional Cluster Resources:
  -m  Memory allocation [Default: 300g]
  -c  Number of CPU cores [Default: 10]
  -t  Walltime limit (HH:MM:SS) [Default: 03:00:00]
```

Here is an example:

```bash
submit-job -m 250g -t 02:00:00 -c 10 -n scanpy_env -e experiment01 -j downsample_20per -s 04_downsample.py 
```
