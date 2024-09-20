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


<img src="/Users/yucheng/Documents/GitHub/Tufts_2024Fall_Training/images/serial_job.png" alt="serial_job" style="zoom:40%;" />



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

Submitting too many jobs 

#### MaxArraySize

To query `MaxArraySize` , you can use 

```
scontrol show conf | grep MaxArraySize
$ scontrol show config | grep -i array
MaxArraySize            = 2000
```




Public Partitions (batch+mpi+largemem+gpu)
CPU: 1000 cores
RAM: 4000 GB
GPU: 10
Preempt Partition (preempt)
CPU: 2000 cores
RAM: 8000 GB
GPU: 20





# Create a contig file for your array tasks (Change the title)
Use R script as an example. 
In most cases, your script will loop through different input parameters, which are usually not number 1-10, 1-100. In this situation, we would like to a config file with input parameters for each job.  (Will revise later)


## Required files

1. **Parameter File:** A file containing the parameters that your array job will iterate through. This file could include different variables or data that each array task will process individually.
2. **Script (R, Shell, Python, etc.):** The main script that will perform the analysis or visualization tasks. While the example here is in R, the same structure applies to other languages like shell, Python, or Perl. Adapt the script according to the specific tool or language you are using for the job.
3. **Wrapper Shell Script:** A simple shell script that sends your jobs to the SLURM scheduler. This script makes it easy to run multiple tasks automatically, with each task using different parameters from the parameter file.



## R Script Example

Here is n example of an R script that generates scatter plots of gene expression based on raw RNA-seq count data:

```r
# Load libraries
library(tidyverse)
library(ggrepel)

# Read in parameters
args <- commandArgs(trailingOnly = TRUE)
gene <- as.character(args[1])
padj <- as.numeric(args[2])

# Subset the gene of interest
dt <- read.table("salmon.merged.gene_counts.tsv", header=T)
d <- dt[match(gene, dt$gene_name),]
d <- gather(d, key = "condition", value = "expression", GFPkd_1:PRMT5kd_3)

# Reformat for ggplot
d_long <- separate(d, col = "condition", into = c("treatment", "replicate"), sep = "_")

# Ggplot to visualize
p <- ggplot(d_long, aes(treatment, expression)) +
        geom_point(size=5, color="steelblue", alpha=0.5) +
        geom_label_repel(aes(label=replicate)) +
        theme_classic() +
        xlab("Treatment") +
        ylab("Gene expression") +
        ggtitle(paste0(gene,": padj ", padj))

# Save plot to a pdf file
ggsave(plot=p, file=paste0(gene, ".pdf"), width=4, height=4)
```



### Script Purpose

This R script creates scatter plots for gene expression levels between control and treated groups from an RNA-seq analysis. It reads in two parameters from the command line: the gene name (`genename`) and the adjusted p-value (`padj`). The input data file is `salmon.merged.gene_counts.tsv`.



## Example Parameter File

Here’s an example of the parameter file (`table.tsv`) used in the job array. Each row contains gene expression information, and the R script will extract specific columns for each job.

```
gene_id baseMean        log2FoldChange  lfcSE   pvalue  padj    genename
ENSG00000078018 1126.709        -2.161184       0.05810824      1.578201e-304   2.054292e-301   MAP2
ENSG00000004799 1224.003        -2.199776       0.06003955      1.17799e-295    1.4154e-292     PDK4
ENSG00000272398 2064.696        1.615232        0.04513554      2.024618e-282   2.258895e-279   CD24
ENSG00000135046 12905.46        -0.8779134      0.02467955      1.349814e-278   1.405606e-275   ANXA1
```

The R script reads the `genename` from column 7 and `padj` from column 6 for each gene.

## Shell Wrapper Script

The following shell script submits the jobs to the SLURM scheduler as an array of tasks. Each task processes a different gene from the parameter file.

```bash
#!/bin/bash
#SBATCH -p preempt  # batch, gpu, preempt, mpi or your group's partition
#SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
#SBATCH -N 1        # Number of nodes
#SBATCH -n 1        # Number of tasks per node
#SBATCH -c 4        # Number of CPU cores per task
#SBATCH --mem=2G    # Memory required per node
#SBATCH --array=2-11 # An array of 10 jobs
#SBATCH --job-name=Rplot
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=xue.li37@tufts.edu
#SBATCH --error=%x-%A_%a.err   # Standard error file: <job_name>-<job_id>-<taskid>.err
#SBATCH --output=%x-%A_%a.out  # Standard output file: <job_name>-<job_id>-<taskid>.out

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

module load R/4.4.0
GENE=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$7}" table.tsv) 
Padj=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$6}" table.tsv)

echo $GENE $Padj
Rscript R_scatter_vis.r $GENE $Padj  
```



### Script Details

- `SBATCH --array=2-11` tells SLURM to run jobs for rows 2 to 11 of the parameter file.
- The `awk` commands extract the `GENE` and `Padj` values from the specified row and columns (7th and 6th).
- The script submits 10 jobs, each running the R script with different `GENE` and `Padj` values.



### Customizing the Array

You can adjust the `--array` option to change the range of jobs. For example, to run jobs for every other line from 2 to 1000, you can specify:

```
#SBATCH --array=2-1000:2
```

This would submit jobs for rows 2, 4, 6, ..., up to 1000.
