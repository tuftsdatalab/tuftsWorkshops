# AlphaFold on Tufts HPC Cluster

### General Tufts HPC Cluster Access Info

Please review  https://tufts.box.com/v/Pax-User-Guide before proceeding forward.

### Login and Allocate Computing Resources

Tips:

1. Login. If you have a Mac, use the Terminal app. If you have a 
2. If you need to use GPU resources and don't have access to contrib node partitions, "preempt" is the best option

### Alphafold

The Alphafold script is available for everyone in `/cluster/tufts/hpc/tools/alphafold/2.2.0/runaf2test.sh`

Make a copy of the file to your own folder  (e.g. your home directory):

`$ cp /cluster/tufts/hpc/tools/alphafold/2.2.0/runaf2test.sh /your/own/directory`

Go to your own copy of the script:

```
#!/bin/bash
#SBATCH -p preempt  #if you DO have ccgpu access, use "ccgpu"
#SBATCH -n 8    # 8 cpu cores
#SBATCH --mem=64g       #64GB of RAM
#SBATCH --time=2-0      #run 2 days, up to 7 days "7-00:00:00"
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH -N 1
#SBATCH --gres=gpu:a100:1    # number of GPUs. please follow instructions in "Pax User Guide" when submit jobs to different partition and selecting different GPU architectures. 

module load alphafold/2.2.0
module list
nvidia-smi

module help alphafold/2.2.0 # this command will print out all input options for "runaf2" command

#Please use your own path/value for the following variables
#Make sure to specify the outputpath to a path that you have write permission
outputpath=/cluster/tufts/hpc/tools/alphafold/2.2.0/test
fastapath=/cluster/tufts/hpc/tools/alphafold/2.2.0/T1050.fasta
maxtemplatedate=2020-06-10

source activate alphafold2.2.0

#running alphafold 2.2.0

runaf2 -o $outputpath -f $fastapath -t $maxtemplatedate

```

Make sure you specify the `outputpath` to a path that you have write permission. 

Make sure you specify the `fastapath` to FASTA file containing the protein sequence for which you wish to predict the structure.

Make sure `maxtemplatedate` is set to be before the release date of the structure. 

### Module Help

Please see `$ module help alphafold/2.2.0` for additional input options ( __Required Parameters__ & __Optional Parameters__) for `runaf2` command.

```
----------- Module Specific Help for 'alphafold/2.2.0' ------------

	This module adds AlphaFold 2.2.0 to the PATH
Run AlphaFold 2.2.0 with:
runaf2 <Required parameters> <Optional Parameters>
Please make sure all REQUIRED parameters are given

Required Parameters:
-o <output_dir>       Path to a directory that will store the results, make sure the user has write permission to the directory.
-f <fasta_path>       Path to a FASTA file containing sequence. If a FASTA file contains multiple sequences, then it will be folded as a multimer
-t <max_template_date> Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets

Optional Parameters:
-g <use_gpu>          Enable NVIDIA runtime to run with GPUs (default: true)
-n <openmm_threads>   OpenMM threads (default: all available cores)
-a <gpu_devices>      Comma separated list of devices to pass to 'CUDA_VISIBLE_DEVICES' (default: 0)
-m <model_preset>     Choose preset model configuration - the monomer model, the monomer model with extra ensembling, monomer model with pTM head, or multimer model (default: 'monomer')
-c <db_preset>        Choose preset MSA database configuration - smaller genetic database config (reduced_dbs) or full genetic database config (full_dbs) (default: 'full_dbs')
-p <use_precomputed_msas> Whether to read MSAs that have been written to disk. WARNING: This will not check if the sequence, database or configuration have changed (default: 'false')
```

### Submit Job

To submit your job, go to the folder that contains `runaf2test.sh`

From command line, submit with `$ sbatch runaf2test.sh`

Then follow the instructions in  https://tufts.box.com/v/Pax-User-Guide to check your job status.

### AlphaFold output

The outputs will be in a subfolder of `output_dir` that you specified in `runaf2test.sh`. 

They include the computed MSAs, unrelaxed structures, relaxed structures, ranked structures, raw model outputs, prediction metadata, and section timings. The`output_dir` directory will have the following structure:

```
output_dir/

  features.pkl

  ranked_{0,1,2,3,4}.pdb

  ranking_debug.json

  relaxed_model_{1,2,3,4,5}.pdb

  result_model_{1,2,3,4,5}.pkl

  timings.json

  unrelaxed_model_{1,2,3,4,5}.pdb

  msas/

​    bfd_uniclust_hits.a3m

​    mgnify_hits.sto

​    uniref90_hits.sto
```

The contents of each output file are as follows:

- `features.pkl` – A `pickle` file containing the input feature Numpy arrays used by the models to produce the structures.

- `unrelaxed_model_*.pdb` – A PDB format text file containing the predicted structure, exactly as outputted by the model.

- `relaxed_model_*.pdb` – A PDB format text file containing the predicted structure, after performing an Amber relaxation procedure on the unrelaxed structure prediction, see Jumper et al. 2021, Suppl. Methods 1.8.6 for details.

- `ranked_*.pdb` – A PDB format text file containing the relaxed predicted structures, after reordering by model confidence. Here `ranked_0.pdb` should contain the prediction with the highest confidence, and `ranked_4.pdb` the prediction with the lowest confidence. To rank model confidence, we use predicted LDDT (pLDDT), see Jumper et al. 2021, Suppl. Methods 1.9.6 for details.

- `ranking_debug.json` – A JSON format text file containing the pLDDT values used to perform the model ranking, and a mapping back to the original model names.

- `timings.json` – A JSON format text file containing the times taken to run each section of the AlphaFold pipeline.

- `msas/` - A directory containing the files describing the various genetic tool hits that were used to construct the input MSA.

- `result_model_*.pkl` – A `pickle` file containing a nested dictionary of the various Numpy arrays directly produced by the model. In addition to the output of the structure module, this includes auxiliary outputs such as distograms and pLDDT scores. If using the pTM models then the pTM logits will also be contained in this file.

