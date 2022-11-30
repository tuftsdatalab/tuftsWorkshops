## Log into the HPC cluster’s On Demand Interface

- Open a Chrome browser and go to [On Demand](https://ondemand.pax.tufts.edu/)
- Log in with your Tufts Credentials
- On the top menu bar choose `Clusters->Tufts HPC Shell Access`

![](./images/shell.png)

- You'll see a welcome message and a bash prompt, for example for user `tutln01`:

```
[tutln01@login001 ~]$
```

- This indicates you are logged in to the login node of the cluster. **DO NOT** run anything from this node as it is a shared login node. To complete jobs we either need to start an interactive session to get on to a compute node (more details can be found [here](https://tufts.box.com/v/Pax-User-Guide) about this option) or we can write a batch script to submit our job to a compute node.

## AlphaFold2 Batch Script

A batch script can be broken into two parts - the header section with information on how to run the job and a command section where we use UNIX commnads to do a job. Here is our script:

```
#!/bin/bash
#SBATCH -p ccgpu                            # partition we submit to
#SBATCH -n 8                                # The number of cpu cores we would like
#SBATCH --mem=64g                           # The amount of RAM we would like
#SBATCH --time=2-0:00:00                    # The time we think our job will take       
#SBATCH -o output.%j                        # The name of the output file
#SBATCH -e error.%j                         # The name of the error file
#SBATCH -N 1                                # The number of nodes we would like
#SBATCH --gres=gpu:1                        # The number of GPUs
#SBATCH --exclude=c1cmp[025-026]            # The nodes to exclude when using AlphaFold2

# Load the AlphaFold2 and NVIDIA modules
module load alphafold/2.1.1
nvidia-smi

# Make the results direcory
mkdir /path/to/your/home/directory/af2

# Specify where your output directory and raw data are
outputpath=/path/to/your/home/directory/af2
fastapath=/path/to/your/home/directory/data/1AXC.fasta

# Date to specify if you want to avoid using template
maxtemplatedate=2020-06-10

source activate alphafold2.1.1

# Running alphafold 2.1.1
runaf2 -o $outputpath -f $fastapath -t $maxtemplatedate -m multimer 
```

!!! example "What do these commands mean?"

    - `module load alphafold/2.1.1`: load the AlphaFold2 module
    - `nvidia-smi`: load the NVIDIA-SMI module
    - `mkdir /path/to/your/home/directory/af2`: make a directory for our outputs
    - `outputpath=/path/to/your/home/directory/af2`: specify where our outputs should go
    - `fastapath=/path/to/your/home/directory/data/1AXC.fasta`: specify our input FASTA file
    - `maxtemplatedate=2020-06-10`: specify our maximum template date
    - `source activate alphafold2.1.1`: activate our AlphaFold2 conda environment
    - `runaf2 -o $outputpath -f $fastapath -t $maxtemplatedate -m multimer`: run our AlphaFold2 program
        - `-o`: output path
        - `-f`: our input FASTA file path
        - `-t`: our maximum template date
        - `-m`: whether or not this is a multimeric protein
    
    
    
!!! info "What's up with the maxtemplatedate?"

    The `maxtemplatedate` option is a bit more complicated. If we ask AlphaFold to predict the structure of a protein with a structure **already** in the [Protein Data Bank (PDB)](https://www.rcsb.org/) - then we have the option of using that structure in the prediction. If we do not want AlphaFold 2 to use this structure in the prediction we need to specify a date **before** the release date of that structure.
    

??? question "What do the SBATCH commands mean?"

    |`command`|description|
    |-|-|
    |`#!/bin/bash`|specify our script is a bash script|
    |`#SBATCH -p ccgpu`| The partition we are requesting, and if you don't have ccgpu access, use "preempt"|
    |`#SBATCH -n 8` | The number of cpu cores we would like - here it is 8|
    |`#SBATCH --mem=64g ` |The amount of RAM we would like - here it is 64 Gigabytes|
    |`#SBATCH --time=2-0:00:00` | The time we think our job will take - here we say 2 days (days-hours:minutes:seconds)|
    |`#SBATCH -o output.%j` |The name of the output file - here it is "output.jobID"|
    |`#SBATCH -N 1`|The number of nodes we would like - here it is 1|
    |`#SBATCH --gres=gpu:1 `|The number of GPUs - here we ask for 1|
    |`#SBATCH --exclude=c1cmp[025-026]`|These are nodes to exclude when using AlphaFold2|

## Run the AlphaFold2 Batch Script

- First save your batch script, here we will save it as `runaf2.sh`. 
    - The bash script must end in `.sh`
- Now submit your script with the command:
    - `sbatch runaf2.sh`
- To check on the status of this script use the following command:
    - `squeue -u your_utln`
    - Be sure to swap out `your_utln` with your Tufts utln!
    

