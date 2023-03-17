## AlphaFold2 Accuracy Assessment

We can assess the accuracy of the AlphaFold prediction using:

- Predicted Local Distance Difference Test (pLDDT)
- Predicted Alignment Error

## Predicted Local Distance Difference Test (pLDDT)

- per-residue confidence metric  ranging from 0-100 (100 being the highest confidence)
- Regions below 50 could indicate disordered regions

![](images/plddt.png)

## Predicted Alignment Error (PAE)

- The Predicted Alignment Error (PAE) gives us an expected distance error based on each residue.
- If we are more confident that the distance between two residues is accurate, then the PAE is lower (darker green). If we are less confident that the distance between two residues is accurate, the PAE is higher (lighter green)

![](images/pae.png)

## Navigate To The Cluster

- Now that we have an idea of what these metrics mean, let's try generating these plots for our Cas12a mutants on the cluster. First navigate to:

[https://ondemand.pax.tufts.edu/](https://ondemand.pax.tufts.edu/)

- Log in with your Tufts credentials
- On the top menu bar choose `Clusters->Tufts HPC Shell Access`

![](images/shell.png)

- You'll see a welcome message and a bash prompt, for example for user `tutln01`:

```
[tutln01@login001 ~]$
```

- This indicates you are logged in to the login node of the cluster. Please **do not** run any program from the login node.

??? question "Are you logged into OnDemand?"

    - Yes
    - No
    
    
## Starting an Interactive Session

- To run our analyses we will need to move from the login node to a compute node. We can do this by entering:

```bash
srun -p batch --time=3:00:00 -n 2 --mem=4g --reservation=bioworkshop --pty bash
```

Where:

!!! example "Explanation of Commands"

    - `srun`: SLURM command to run a parallel job
    - `-p`: asking for a partition, here we are requesting the batch partition
    - `--time`: time we need here we request 3 hours
    - `-n`:  number of CPUs needed here we requested 2
    - `--mem`:  memory we need here we request 4 Gigabytes
    - `--reservation`: the reservation of compute resources to use here we use the `bioworkshop` reservation
    - `--pty`: get a pseudo bash terminal
    
!!! warning 
    
    The `bioworkshop` reservation will be unavailable after December 7th. This reservation is a temporary reservation for this class. 

- When you get a compute node you'll note that your prompt will no longer say login and instead say the name of the node:

```
[tutln01@c1cmp048 ~]$
```

## Set Up For Analysis

- To get our AlphaFold data we will enter:

```bash
cp -r /cluster/tufts/bio/tools/training/af2Workshop ./
```
