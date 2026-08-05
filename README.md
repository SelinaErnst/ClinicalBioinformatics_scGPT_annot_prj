|Folder | Comments|
|-----|---------|
|experiment01| Practicing with relevant python packages|

# job submission

### PBS 

For the kitz cluster (PBS), I created a bash script that can outomatically run jobs and create log files.  
To make it even easier to use, it is helpful to write a line to `~/.bashrc`

```
alias submit-job='/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/job-submission/submit-pbs-job.sh'
```

This way, you can use the command `submit-job` to automatically run scripts from specified folders.

```
Error: Both -e (experiment) and -s (script) are required.
Usage: submit-pbs-job.sh -e <experiment> -s <script.py> [-n <env_name>] [-j <job_name>]

Required Parameters:
  -e  Experiment folder name (EXP)
  -s  Python script filename to run (SCRIPT)

Optional Parameters:
  -j  Custom PBS job name
  -n  Conda environment name (ENV) [Default: scanpy_env]
  -p  Directory to project (PRJ_DIR) [Default: /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project]

Optional Cluster Resources:
  -m  Memory allocation [Default: 250g]
  -c  Number of CPU cores (ppn) [Default: 10]
  -t  Walltime limit (HH:MM:SS) [Default: 02:00:00]
```

Here is an example:

```bash
submit-job -m 400g -t 02:00:00 -c 10 -n scanpy_env -e experiment01 -j downsample_20per -s 04_downsample.py 
```