# Getting to Know Tufts HPC Cluster Resources



## Account & Storage Requests

Go to **[Tufts HPC website](https://it.tufts.edu/high-performance-computing)** for HPC Cluster Account and Storage Requests

For Faculties ONLY: **[HPC for Classes](https://tufts.qualtrics.com/jfe/form/SV_d7o0UZFgK1PFXnv)**



## Tufts HPC Cluster

- Moved to **[MGHPCC](https://www.mghpcc.org/)** in January 2024!



## Cluster Storage

* __Home Directory__

Every user has a home directory.

Be aware! Your Home Directory (**30GB**, fixed) should be `/cluster/home/your_utln`

If you are not sure how much storage you have used in your home directory, feel free to contact us and we can provide you with the number. 

For self-service, you can use the following commands from a shell terminal to find out your home directory usage:

`$ module load hpctools`

`$ hpctools` (from any node) 

OR 

`$ du -a -h --max-depth=1 ~ | sort -hr` from a **compute node** in your home directory to find out the detailed usage. 

* __Reserach Project Storage__

**[Request research project storage](https://it.tufts.edu/research-technology/)**

Created for research labs and classes. A Tufts HPC cluster research project storage share can **only** be owned by a **Tufts** **faculty**.

New Storage Policies can be found on [RT Announcements](https://it.tufts.edu/research-technology/announcements) page. - Tiered Storage

Your research projet storage (from **50GB**) path should be `/cluster/tufts/yourlabname/`, and each member of the lab group has a dedicated directory `/cluster/tufts/yourlabname/your_utln`

To see your **research project storage quota** by running the following command from **any node on the new cluster Pax**:

`$ df -h /cluster/tufts/yourlabname ` 

OR 

`$ module load hpctools`

`$ hpctools`

**NOTE:** Accessing your research project storage space for the __first time__ in your current session, please make sure you type out the __FULL PATH__ to the directory `/cluster/tufts/yourlabname/`.



## CPUs

Primarily Intel Xeon CPUs, from Broadwell to Emerald Rapids, with hyperthreading enabled*.

Compute nodes are grouped into **partitions** based on <u>functionality</u> and <u>priority</u> levels.

**Public Partitions :**

```
PARTITION       TIMELIMIT      
batch*          7-00:00:00          
gpu             7-00:00:00        
interactive     4:00:00        
largemem        7-00:00:00        
mpi             7-00:00:00         
preempt         7-00:00:00     
```

* **preempt** - Be aware, the `preempt` partition consists of most of the nodes on the cluster, including public nodes and **contrib nodes** from different research labs. When submitting jobs to preempt partition, you acknowledge that your jobs are taking the **risk** of being preempted by higher priority jobs. In that case, you will simply have to resubmit your jobs. 

  

## GPUs

__NVIDIA GPUs__ are available in `gpu` and `preempt` partitions

- **Request GPU resources with `--gres`. See details below.**

- If no specific architecture is required, GPU resources can be request with`--gres=gpu:1` (one GPU)

- You can request more than one type of GPUs with `constraint`, e.g.  

  `--gres=gpu:1 --constraint="t4|p100|v100"`

- Please **DO NOT** manually set `CUDA_VISIBLE_DEVICES`. 

- Users can ONLY see GPU devices that are assigned to them with `$ nvidia-smi`.

  If you submit batch jobs, it's recommended adding `nvidia-smi` in your slurm job submission script.

- `gpu` partition`-p gpu`:

  - NVIDIA A100
    - In "gpu" partition
    - Request with: `--gres=gpu:a100:1`(one A100 GPU, can request up to 8 on one node)
    - `--constraint="a100-80G"`
    - Each GPU comes with 80GB of DRAM
    - Driver supports upto CUDA 12.2
  - NVIDIA P100s
    - In "gpu" partition
    - Request with: `--gres=gpu:p100:1`(one P100 GPU, can request up to 4 on one node)
    - `--constraint="p100"`
    - Each GPU comes with 16GB of DRAM
    - Driver supports upto CUDA 12.2

- `preempt` partition `-p preempt`:

  - `a100`, `v100`, `p100`, ` rtx_6000`, `rtx_a6000`, `rtx_6000ada`, `rtx_a5000`, `h100`, `l40s`, `t4`

  - NVIDIA T4

    - In "preempt" partition
    - Request with: `--gres=gpu:t4:1`(one T4 GPU, can request up to 4 on one node)
    - `--constraint="t4"`
    - Each GPU comes with 16GB of DRAM
    - Driver supports upto CUDA 10.2

  - NVIDIA P100

    - In "preempt" partition
    - Request with: `--gres=gpu:p100:1`(one P100 GPU, can request up to 6 on one node)
    - `--constraint="p100"`
    - Each GPU comes with 16GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA rtx_6000

    - In "preempt" partition
    - Request with: `--gres=gpu:rtx_6000:1`(one RTX_6000 GPU, can request up to 8 on one node)
    - `--constraint="rtx_6000"`
    - Each GPU comes with 24GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA rtx_a6000

    - In "preempt" partition
    - Request with: `--gres=gpu:rtx_a6000:1`(one RTX_A6000 GPU, can request up to 8 on one node)
    - `--constraint="rtx_a6000"`
    - Each GPU comes with 48GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA rtx_6000ada

    - In "preempt" partition
    - Request with: `--gres=gpu:rtx_6000ada:1`(one RTX_6000Ada GPU, can request up to 4 on one node)
    - `--constraint="rtx_6000ada"`
    - Each GPU comes with 48GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA V100

    - In "preempt" partition
    - Request with: `--gres=gpu:v100:1`(one V100 GPU, can request up to 4 on one node)
    - `--constraint="v100"`
    - Each GPU comes with 16GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA A100

    - In "preempt" partition
    - Request with: `--gres=gpu:a100:1`(one A100 GPU, can request up to 8 on one node)
    - `--constraint="a100-80G"`
    - `--constraint="a100-40G"`
    - `--constraint="a100"`
    - Each GPU comes with 40GB of DRAM or 80GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA H100

    - In "preempt" partition
    - Request with: `--gres=gpu:h100:1`(one V100 GPU, can request up to 3 on one node)
    - `--constraint="h100"`
    - Each GPU comes with 80GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA L40s

    - In "preempt" partition
    - Request with: `--gres=gpu:l40s:1`(one L40s GPU, can request up to 4 on one node)
    - `--constraint="l40s"`
    - Each GPU comes with 48GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA RTX A5000

    - In "preempt" partition
    - Request with: `--gres=gpu:rtx_a5000:1`(one RTX A5000 GPU, can request up to 4 on one node)
    - `--constraint="l40s"`
    - Each GPU comes with 48GB of DRAM
    - Driver supports upto CUDA 12.2

  - NVIDIA L40

    - In "preempt" partition

    - Request with: `--gres=gpu:l40:1`(one L40 GPU, can request up to 4 on one node)

    - `--constraint="l40"`

    - Each GPU comes with 48GB of DRAM

    - Driver supports upto CUDA 12.2

      

## Cluster Resource Limit

* **Public Partitions** (batch+mpi+largemem+gpu) 

  * CPU: 500 cores

    RAM: 2000 GB

    GPU: 5

* **Preempt Partition** (preempt) 

  * CPU: 1000 cores

    RAM: 4000 GB

    GPU: 10

- **Maximun Number of Jobs submitted**: 1000

  

## Cluster Computing Resource Availability

`$ module load hpctools`

`$ hpctools` 

Then follow the on-screen instructions to extract the information you need. 


---

