# Getting Work Done on Tufts HPC Cluster with SLURM

<img src="https://ormstonhouse.com/wp-content/uploads/2019/09/tetris.jpg" alt="teris" width=50%>



**ALL** work **MUST** to be performed on **compute nodes**!

If you see prompt like this 
`[your_utln@login-prod-01]`
`[your_utln@login-prod-02]`
`[your_utln@login-prod-03]`
**DON'T** run any programs, YET! 
**Get resource allocation first**!

<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/Cluster_20230516.png" alt="Memory-Storage" width=70%>

> **Things to think about before requesting resources:**
>
> - What program?
> - What kind of resources? CPUs only or with GPUs?
> - How much memory?
> - How many CPU cores?
> - How long does the job need to run?
> - Command line application? GUI?
> - Does the program need interactions?
> - Prototyping, debugging, or production?

## OnDemand

- **[OnDemand (https://ondemand.pax.tufts.edu)]( https://ondemand.pax.tufts.edu/)**

  - **`Interactive Apps`** --> RStudio, Matlab, JupyterLab, Jupyter Notebook, .etc

  - **`Clusters`** --> Tufts HPC Cluster  Shell Access 

  - **`Files`** 

  - **`Jobs`**

  - [OnDemand]( https://ondemand.pax.tufts.edu/)  **`Clusters`** --> Tufts HPC Cluster FastX11 Shell Access

    

## Slurm Information

- View information about Slurm nodes and partitions.

`$ sinfo`

With more specifc information and formated output:

`$ sinfo -o "%20N %10P %10c %10m %85f %10G "` - NODELIST, PARTITION, CPUS, MEMORY,AVAIL_FEATURES, GRES  

More  [sinfo](https://slurm.schedmd.com/sinfo.html) options

**You can only see the partitions you have access to.**

**For most users**, you will see `batch`,` mpi`,` gpu`,` largemem`, and `preempt` partitions.          

How to check **GPU, Memory, CPU** availability on the cluster?



## Utilize New `hpctools` Module !!!

Users can use `hpctools` module to check: **Free CPU resources, Free GPU resources, User Past and Active jobs, and Project space quota and usage. **

```[tutln01@login-prod-01 ~]$ module load hpctools
[tutln01@login-prod-01 ~]$ module load hpctools
	 command: hpctools
[tutln01@login-prod-01 ~]$ hpctools
 Please select from the following options:

  1. Checking Free Resources On Each Node in Given Partition(s)

  2. Checking Free GPU Resources On Each Node in Given Partition(s)

  3. Checking tutln01 Past Completed Jobs in Given Time Period

  4. Checking tutln01 Active Job informantion

  5. Checking Project Space Storage Quota Informantion

  6. Checking Any Directory Storage Usage Informantion

  Press q to quit

Your Selection:

```

Then follow the onscreen instructions to get desired information.



## Interactive Session

- Particularly good for debugging and working with software GUI. 

  `$ srun [options] --pty [command]`

- Command 

  - command to run an application, given the module is already loaded.
  - `bash` for a bash shell

- Options

  - Pseudo terminal `--pty`
  - Partition `-p` or `--partition=`
    - Default batch if not specified
  - Time `-t` or `--time=`
    - Default 15 minutes if not specified on non-interactive partition
    - Format: DD-HH:MM:SS
  - Number of CPU cores `-n` 
    - Default 1 if not specified
  - Memory `--mem=`
    - Default 2GB if not specified
  - GPU `--gres=`
    - Default none
  - Features `--constraint=`
    - GPU types
    - OS version
    - CPU architecture
    - Instruction Set
    - Default none
  - X Window `--x11=first`
    - Default none	

**Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 2GB of RAM, with X11 forwarding for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).**

```bash
[tutln01@login-prod-01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=2g --x11=first --pty bash
srun: job 296794 queued and waiting for resources
srun: job 296794 has been allocated resources
[tutln01@cc1gpu001 ~]$ 
```

Note: If you are requesting X11 forwarding with `srun`, `-XC` or`-YC` or `-XYC` must be used upon login with `ssh`.

**Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 4GB of RAM, with 1 A100 GPU for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).**

```bash
[tutln01@login-prod-01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=4g --gres=gpu:a100:1 --pty bash
```

Once your resource is allocated on a compute node, use `nvidia-smi` to check GPU info.



## Batch Job

Write a batch submission script e.g. **myjob.sh**

```bash
#!/bin/bash -l
#SBATCH -J My_Job_Name   #job name
#SBATCH --time=00-00:20:00  #requested time (DD-HH:MM:SS)
#SBATCH -p batch,preempt    #running on "batch" or "preempt" partition, wherever resource is available first
#SBATCH -N 1    #1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2   #2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --mem=2g  #requesting 2GB of RAM total for the number of cpus you requested
##SBATCH --gres=gpu:a100:1	#requesting 1 A100 GPU, in this case, the "-p" needs to be switched to a partition has the requested GPU resources
#SBATCH --output=MyJob.%j.%N.out  #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=MyJob.%j.%N.err   #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL    #email optitions
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]
# unload all mo
module purge
# for example, running a python script 
# load the module so the correct version python is available to you
module load anaconda/2021.05
# If you have a conda env that you would like to use, activate it here using "source activate xxx". DO NOT USE "conda activate"
source activate [target_env]
# run python script
python myscript.py #make sure myscript.py exists in the current directory
# make sure you save all plots, data, outputs generated to files in your script
# Don't forget to deactivate your conda env if you are using one
conda deactivate
```

**Submit** the job using the following command from command line interface:

`$ sbatch myjob.sh`

**Sample Scripts** including R, conda, matlab, gaussian

`/cluster/tufts/hpc/tools/slurm_scripts`



## Job Status

- Checking your **active** jobs

```bash
[tutln01@cc1gpu001 ~]$ squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
            296794   preempt     bash tutln01  R       5:12      1 cc1gpu001 
[tutln01@cc1gpu001 ~]$ squeue -u tutln01
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON) 
            296794   preempt     bash tutln01  R       5:21      1 cc1gpu001 
```

To check your active jobs in the queue:

`$ squeue -u $USER` or `$ squeue -u your_utln`

To cancel a specific job:

`$ scancel JOBID`

To cancel all of your jobs:

`$ scancel -u $USER` or `$ scancel -u your_utln`

To check details of your **active jobs** (running or pending):

`$ scontrol show jobid -dd JOBID`

```bash
[tutln01@cc1gpu001 ~]$ scontrol show jobid -dd 296794
JobId=296794 JobName=bash
   UserId=tutln01(31003) GroupId=tutln01(5343) MCS_label=N/A
   Priority=10833 Nice=0 Account=(null) QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=0 Restarts=0 BatchFlag=0 Reboot=0 ExitCode=0:0
   DerivedExitCode=0:0
   RunTime=00:10:33 TimeLimit=1-02:30:00 TimeMin=N/A
   SubmitTime=2021-03-22T22:18:50 EligibleTime=2021-03-22T22:18:50
   AccrueTime=2021-03-22T22:18:50
   StartTime=2021-03-22T22:18:55 EndTime=2021-03-24T00:48:55 Deadline=N/A
   PreemptEligibleTime=2021-03-22T22:18:55 PreemptTime=None
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2021-03-22T22:18:55
   Partition=preempt AllocNode:Sid=login-prod-01:34458
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=cc1gpu001
   BatchHost=cc1gpu001
   NumNodes=1 NumCPUs=2 NumTasks=2 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=2,mem=2G,node=1,billing=2
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   JOB_GRES=(null)
     Nodes=cc1gpu001 CPU_IDs=30-31 Mem=2048 GRES=
   MinCPUsNode=1 MinMemoryNode=2G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=bash
   WorkDir=/cluster/home/tutln01
   Power=
   MailUser=tutln01 MailType=NONE
```

- Checking your **finished** jobs

*You can no longer see these jobs in `squeue` command output.*

Querying finished jobs helps users make better decisions on requesting resources for future jobs. 

**Check job CPU, memory usage, and efficiency:**

`$ seff JOBID`

```bash
[tutln01@login-prod-01 ~]$ seff 296794
Job ID: 296794
Cluster: pax
Use of uninitialized value $user in concatenation (.) or string at /usr/bin/seff line 154, <DATA> line 602.
User/Group: /tutln01
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 2
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:22:12 core-walltime
Job Wall-clock time: 00:11:06
Memory Utilized: 1.16 MB (estimated maximum)
Memory Efficiency: 0.06% of 2.00 GB (2.00 GB/node)
```

**Check job detailed accounting data:**

`$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j JOBID`

```bash
[tutln01@login-prod-01 ~]$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j  296794
 Partition      State  Timelimit               Start                 End    Elapsed     MaxRSS     ReqMem  MaxVMSize   NNodes      NCPUS        NodeList 
---------- ---------- ---------- ------------------- ------------------- ---------- ---------- ---------- ---------- -------- ---------- --------------- 
   preempt  COMPLETED 1-02:30:00 2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06                   2Gn                   1          2       cc1gpu001 
           OUT_OF_ME+            2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06         8K        2Gn    135100K        1          2       cc1gpu001 
            COMPLETED            2021-03-22T22:18:56 2021-03-22T22:30:01   00:11:05       592K        2Gn    351672K        1          2       cc1gpu001 
```

NOTE: there are more format options, see [sacct](https://slurm.schedmd.com/sacct.html)

**OR** utilize `hpctools` module on the cluster to make things a little easier

```[tutln01@login-prod-01 ~]$ module load hpctools
[tutln01@login-prod-01 ~]$ module load hpctools
	 command: hpctools
[tutln01@login-prod-01 ~]$ hpctools
 Please select from the following options:

  1. Checking Free Resources On Each Node in Given Partition(s)

  2. Checking Free GPU Resources On Each Node in Given Partition(s)

  3. Checking tutln01 Past Completed Jobs in Given Time Period

  4. Checking tutln01 Active Job informantion

  5. Checking Project Space Storage Quota Informantion

  6. Checking Any Directory Storage Usage Informantion

  Press q to quit

Your Selection:

```

Then follow the onscreen instructions to get desired information.



**Useful Link**

[https://tufts.box.com/v/HPC-Table-of-Contents](https://tufts.box.com/v/HPC-Table-of-Contents)
