# Slurm Job Arrays

Job arrays offer a mechanism for submitting and managing collections of similar jobs quickly and easily, saving both time and computational resources.

## Use cases

- I have 1000 samples and they all need to run the same workflow.
- I need to run a simulation 1000 times with a different set of parameters.


**Why not use serial jobs?** 

A common approach is to use bash loops to submit jobs one by one, but this is not efficient for large numbers of tasks. For example:

```
for fq in *.fastq.gz; do 
  fastqc -t 4 $fq
done
```
<img src="../images/serial_job.png" alt="serial_job" style="zoom:40%;" />



Using bash loops works but often results in jobs taking much longer. **Instead, using SLURM job arrays can streamline this process.**





## Slurm arrays

### Basic Syntax

Job arrays are only supported for batch jobs, and the array index values are specified using the `--array` or `-a` option of the `sbatch` command or `#SBATCH` inisde job script. 

```
--array=<indices>
```

- You can specify the array indices in different ways:

  - `--array=0-100`: Runs jobs with indices from 0 to 100.

  - `--array=2,4,6,8,10`: Runs jobs with specific indices (2, 4, 6, 8, and 10).

  - `--array=2-1000:2`: Runs jobs with a step size, in this case, every 2nd job from 2 to 1000.

- You can limit the number of array jobs which are allowed to run at once by using the `%` character when specifying indices.

  - `1-16%2` Create 16 jobs, but only allow two to run at a time


### Job ID and Environment Variables 

#### SLURM_ARRAY_JOB_ID

- This environment variable represents the job ID of the entire job array.

- It is the same for all tasks within that job array.

- If you submit a job array with 10 tasks, each of those tasks will have the same `SLURM_ARRAY_JOB_ID`.

**Example**:

If you submit a job array with `sbatch --array=1-10 script.sh`, and the job array is assigned the job ID **12345**, then

- `SLURM_ARRAY_JOB_ID` for all tasks will be **12345**.

#### **SLURM_ARRAY_TASK_ID**

- This environment variable represents the unique identifier of each task within the job array.

- It differentiates each task in the array and usually corresponds to the index you specified when submitting the job array.

- **This is the variable you use to handle task-specific operations within the script**.

Example:

If you submit a job array with `sbatch --array=1-10 script.sh`, and the job array is assigned the job ID **12345**, then:

- Task 1 will have SLURM_ARRAY_TASK_ID=1.

- Task 2 will have SLURM_ARRAY_TASK_ID=2.

- And so on, up to SLURM_ARRAY_TASK_ID=10 for the last task.

In a simple case, you can directly use the `$SLURM_ARRAY_TASK_ID` variable in your script to set up your job array. 

For instance, if you have a fasta file for each sample like: sample1.fa, sample2.fa, sample3.fa ... sample10.fa, and you want each of the 10 Slurm array tasks to handle a separate sample file, you can replace the line specifying the sample filename with `sample${SLURM_ARRAY_TASK_ID}.fa`. 

This means that for array task 1, the script will run sample1.fa, for array task 2 it will run sample2.fa, and so on.

### Monitor and cancel jobs

You can cancel a particular array task using the respective JOBID in the first column, e.g. `scancel 7456478_2`, or you can cancel all array tasks in the array job by just specifying the main job ID, e.g. `scancel 7456478`.

```
[yzhang85@login-prod-01 array]$ sbatch fastqc_array.sh 
Submitted batch job 7456347
[yzhang85@login-prod-01 array]$ squeue --me
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         7456478_1   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
         7456478_2   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
         7456478_3   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
         7456478_4   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
         7456478_5   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
         7456478_6   preempt   fastqc yzhang85  R       1:30      1 s1cmp004
```

### Limiting the number of tasks to run simultaneously
By default, if sufficient resources are available, all tasks in a job array will run simultaneously. However, if you wish to limit the number of tasks running at once, you can use the `%N` parameter with the `--array` option (where N specifies the maximum number of tasks to execute concurrently). 

In the following example, I used `--array=0-999%10`, which creates a total of 1000 tasks. By appending `%10`, I limited the number of tasks that can run concurrently to 10, meaning that instead of all 1000 tasks running simultaneously, only 10 tasks will be executed at any given time. This helps manage resource usage on the cluster.

```
         JOBID       PARTITION  NAME     USER  ST       TIME  NODES NODELIST(REASON)
7689847_[10-999%10   preempt array_te yzhang85 PD       0:00      1 (JobArrayTaskLimit)
         7689847_0   preempt array_te yzhang85  R       0:03      1 p1cmp029
         7689847_1   preempt array_te yzhang85  R       0:03      1 p1cmp029
         7689847_2   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_3   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_4   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_5   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_6   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_7   preempt array_te yzhang85  R       0:03      1 p1cmp043
         7689847_8   preempt array_te yzhang85  R       0:03      1 p1cmp045
         7689847_9   preempt array_te yzhang85  R       0:03      1 p1cmp045
```

### Example job scripts

In the following example, I have many `fastq.gz` files in the folder `fastq`. I need to run `fastqc`  to check the quality of each of these `fastq.gz` files. 

```
$ ls -1 fastq/*.gz
fastq/SRX1693951_1.fastq.gz
fastq/SRX1693951_2.fastq.gz
fastq/SRX1693952_1.fastq.gz
fastq/SRX1693952_2.fastq.gz
fastq/SRX1693953_1.fastq.gz
fastq/SRX1693953_2.fastq.gz
fastq/SRX1693954_1.fastq.gz
fastq/SRX1693954_2.fastq.gz
fastq/SRX1693955_1.fastq.gz
fastq/SRX1693955_2.fastq.gz
fastq/SRX1693956_1.fastq.gz
fastq/SRX1693956_2.fastq.gz
```

For each `fastq.gz` file, I want to submit a separate slurm job to our cluster. This can be achieved with slurm job array. 

We can use `fastq/SRX169395${SLURM_ARRAY_TASK_ID}_1.fastq.gz` and `fastq/SRX169395${SLURM_ARRAY_TASK_ID}_2.fastq.gz` to represent pairs of `fastq.gz` files.

```
#!/bin/bash
#SBATCH -p preempt  # batch, gpu, preempt, mpi or your group's own partition
#SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
#SBATCH -N 1   # Number of nodes
#SBATCH -n 1   # Number of tasks per node
#SBATCH -c 4   # Number of CPU cores per task
#SBATCH --mem=8G       # Memory required per node
#SBATCH --array=1-6		# An array of 6 jobs
#SBATCH --job-name=fastqc      # Job name
#SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
#SBATCH --mail-user=yzhang85@tufts.edu       # Email address for notifications
#SBATCH --error=%x-%A_%a.err   # Standard error file: <job_name>-<job_id>-<taskid>.err
#SBATCH --output=%x-%A_%a.out  # Standard output file: <job_name>-<job_id>-<taskid>.out

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load fastqc/0.12.1
fastqc -t 4 fastq/SRX169395${SLURM_ARRAY_TASK_ID}_1.fastq.gz fastq/SRX169395${SLURM_ARRAY_TASK_ID}_2.fastq.gz -o fastqcOut
```



Output logs

```
[yzhang85@login-prod-01 array]$ ls -hl
total 13K
drwxrws--- 2 yzhang85 workshop 4.0K Aug 30 11:51 fastq/
drwxrws--- 2 yzhang85 workshop 4.0K Aug 30 11:39 fastqcOut/
-rw-rw---- 1 yzhang85 workshop 1.2K Aug 30 11:54 fastqc-7456478_1.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_1.out
-rw-rw---- 1 yzhang85 workshop 1.1K Aug 30 11:54 fastqc-7456478_2.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_2.out
-rw-rw---- 1 yzhang85 workshop 1.1K Aug 30 11:54 fastqc-7456478_3.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_3.out
-rw-rw---- 1 yzhang85 workshop 1.1K Aug 30 11:54 fastqc-7456478_4.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_4.out
-rw-rw---- 1 yzhang85 workshop 1.1K Aug 30 11:54 fastqc-7456478_5.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_5.out
-rw-rw---- 1 yzhang85 workshop 1.1K Aug 30 11:54 fastqc-7456478_6.err
-rw-rw---- 1 yzhang85 workshop  110 Aug 30 11:52 fastqc-7456478_6.out
-rw-rw---- 1 yzhang85 workshop  918 Aug 30 11:48 fastqc_array.sh
```

## Limits

HPC is a valuable shared resource that allows many users to perform complex calculations simultaneously. To ensure a productive and fair environment for everyone, we have implemented policies and practices that promote equitable access to our computational resources.

There are several limits for array jobs. 
If you submit too many array jobs and exceed the limits, you will get the below error message:
```
$ sbatch array.sub 
sbatch: error: AssocMaxSubmitJobLimit
sbatch: error: Batch job submission failed: Job violates accounting/QOS policy (job submit limit, user's size and/or time limits)
```

#### MaxArraySize

The array index should be smaller than `MaxArraySize`. 

```
scontrol show conf | grep MaxArraySize
$ scontrol show config | grep -i array
MaxArraySize            = 2000
```

Since `MaxArraySize` is set as 2000, the maximum array index you can use is **1999**. So "1000-1999" is valid, but "1001-2000" is invalid. 

#### MaxSubmit
Our cluster does not allow users to submit > 1000 jobs. As a result, the maximum array size is 1000. So "0-999" and "1-1000" is valid, but "1-1001" or "0-1000" is invalid. 

#### CPUs, RAM and GPUs

```
Public Partitions (batch+mpi+largemem+gpu)
CPU: 1000 cores
RAM: 4000 GB
GPU: 10

Preempt Partition (preempt)
CPU: 2000 cores
RAM: 8000 GB
GPU: 20
```

Please note that the above limits are subject to change in the future. To ensure optimal resource allocation, the limit value is dynamic and may change as we evaluate system demands.



## Using parameter file to manage array tasks


In most cases, your script will loop through different input parameters, which are usually not number 1-10 or 1-100. 

In this situation, we would like to use a parameter file with input parameters for each job.  

### Required files

1. **Parameter File:** A file containing the parameters that your array job will iterate through. This file could include different variables or data that each array task will process individually.
3. **Slurm Script:** A simple shell script that sends your jobs to the SLURM scheduler. This script makes it easy to run multiple tasks automatically, with each task using different parameters from the parameter file.



### Fastqc Command 

Here is an example of fastqc command that generates html report for each pair of fastq files

```r
fastqc ${fastq_R1} ${fastq_R2} -o ${output_folder}
```



### Parameter File

Here’s an example of the parameter file `id_sample.tsv` used in the job array. Each row includes a sample ID, along with the corresponding forward and reverse read FASTQ files.

```
1     SRX1693953_SRR3362663_1.fastq.gz SRX1693953_SRR3362663_2.fastq.gz
2     SRX1693956_SRR3362666_1.fastq.gz SRX1693956_SRR3362666_2.fastq.gz
3     SRX1693952_SRR3362662_1.fastq.gz SRX1693952_SRR3362662_2.fastq.gz
4     SRX1693955_SRR3362665_1.fastq.gz SRX1693955_SRR3362665_2.fastq.gz
5     SRX1693951_SRR3362661_1.fastq.gz SRX1693951_SRR3362661_2.fastq.gz
6     SRX1693954_SRR3362664_1.fastq.gz SRX1693954_SRR3362664_2.fastq.gz
```


### Slurm Script

The following shell script submits the jobs to the SLURM scheduler as an array of tasks. Each task processes a sample. 

```bash
#!/bin/bash
#SBATCH -p preempt  # batch, gpu, preempt, mpi or your group's partition
#SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
#SBATCH -N 1        # Number of nodes
#SBATCH -n 1        # Number of tasks per node
#SBATCH -c 4        # Number of CPU cores per task
#SBATCH --mem=8G    # Memory required per node
#SBATCH --array=1-6 # An array of 10 jobs
#SBATCH --job-name=fastqc
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=utln@tufts.edu
#SBATCH --error=%x-%A_%a.err   # Standard error file: <job_name>-<job_id>-<taskid>.err
#SBATCH --output=%x-%A_%a.out  # Standard output file: <job_name>-<job_id>-<taskid>.out

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load fastqc/0.12.1
ID=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" id_sample.tsv) 
fastq_R1=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" id_sample.tsv) 
fastq_R2=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" id_sample.tsv)

echo $ID ${fastq_R1} ${fastq_R2}
fastqc ${fastq_R1} ${fastq_R2} -o fastqcOut

```

#### Script Details

- `SBATCH --array=1-6` tells SLURM to run jobs for rows 1 to 6 of the parameter file.
- The `awk` commands extract the `ID`, `fastq_R1` and `fastq_R2` values from the specified row and columns (1st, 2nd and 3rd).
- The script submits 6 jobs, each running the FASTQC command with different `fastq_R1` and `fastq_R2` values.



After the array jobs are submitted, we can see that 6 separate jobs are running, with ***\*SLURM_ARRAY_TASK_ID\**** from 1-6.

```
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         7975940_1   preempt   fastqc     utln  R       4:35      1 d1cmp024
         7975940_2   preempt   fastqc     utln  R       4:35      1 d1cmp032
         7975940_3   preempt   fastqc     utln  R       4:35      1 d1cmp031
         7975940_4   preempt   fastqc     utln  R       4:35      1 d1cmp031
         7975940_5   preempt   fastqc     utln  R       4:35      1 d1cmp022
         7975940_6   preempt   fastqc     utln  R       4:35      1 d1cmp022
```



#### Output files

```
-rw-rw---- 1 yzhang85 workshop 386K Oct 22 14:40 SRX1693951_SRR3362661_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 590K Oct 22 14:40 SRX1693951_SRR3362661_1_fastqc.html
-rw-rw---- 1 yzhang85 workshop 386K Oct 22 14:40 SRX1693954_SRR3362664_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 588K Oct 22 14:40 SRX1693954_SRR3362664_1_fastqc.html
-rw-rw---- 1 yzhang85 workshop 379K Oct 22 14:40 SRX1693953_SRR3362663_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 585K Oct 22 14:40 SRX1693953_SRR3362663_1_fastqc.html
-rw-rw---- 1 yzhang85 workshop 385K Oct 22 14:41 SRX1693952_SRR3362662_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 590K Oct 22 14:41 SRX1693952_SRR3362662_1_fastqc.html
-rw-rw---- 1 yzhang85 workshop 383K Oct 22 14:41 SRX1693956_SRR3362666_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 588K Oct 22 14:41 SRX1693956_SRR3362666_1_fastqc.html
-rw-rw---- 1 yzhang85 workshop 386K Oct 22 14:41 SRX1693955_SRR3362665_1_fastqc.zip
-rw-rw---- 1 yzhang85 workshop 589K Oct 22 14:41 SRX1693955_SRR3362665_1_fastqc.html 
```



### Customizing the Array

You can adjust the `--array` option to change the range of jobs. For example, to run jobs for every other line from 2 to 1000, you can specify:

```
#SBATCH --array=2-1000:2
```

This would submit jobs for rows 2, 4, 6, ..., up to 1000.



## Useful links:  

https://blog.ronin.cloud/slurm-job-arrays/ 

