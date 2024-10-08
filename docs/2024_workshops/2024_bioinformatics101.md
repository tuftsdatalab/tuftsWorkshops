# Introduction to Tufts HPC Cluster



>Date: 10/09/2024
>
>Delilah Maloney
>Sr. HPC Specialist
>
>###### TTS - Research Technology



Questions, Issues, Requests? 

**Contact us at: tts-research@tufts.edu**

**[Research Technology Website](https://it.tufts.edu/researchtechnology.tufts.edu)**



[TOC]

---

# Tufts HPC Terminologies



## What is "HPC"

- High-performance computing (HPC) is the ability to process data and perform complex calculations at high speeds.

- Everything is relative

## What is "Cluster"

- A computer cluster is a set of loosely or tightly **connected** **computers** that work together so that, in many respects, they can be viewed as a single system. Computer clusters have each node set to perform the same/similar tasks, controlled and scheduled by **[software](#SLURM)**. 

  <img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/Cluster.png" alt="Cluster" width=70%>



***To parallel, or not to parallel? That is the question.***

## CPU vs GPU

- CPU -- **Central Processing Unit** 

  - A CPU can never be fully replaced by a GPU

  - Can be thought of as the **taskmaster** of the entire system, coordinating a wide range of general-purpose computing tasks

- GPU -- **Graphics Processing Unit**

  - GPUs were originally designed to create images for computer graphics and video game consoles

  - GPGPU

  - Performing a narrower range of more specialized tasks

    <img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/CPUGPU.png" alt="CPU-GPU" width=60%>

    

## Cores vs Node

- A **node** is a single computer in the system, which has a number of computing units, **cores**. 

  <img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/CoreNode.png" alt="Core-Node" width=70%>

## Memory vs Storage

The central processing unit (CPU) of a computer is what manipulates data by performing computations.  In practice, almost all computers use a [storage hierarchy](https://en.wikipedia.org/wiki/Memory_hierarchy), which puts fast but expensive and small storage options close to the  CPU and slower but less expensive and larger options further away.  Generally the fast volatile technologies (which lose data when off  power) are referred to as "memory", while slower persistent technologies are referred to as "storage".

- Memory
  - Small, fast, expensive
  - Used to store information for immediate use
  - Volatile.
- Storage
  - Larger, slower, cheaper
  - Non-volatile (retaining data when its power is shut off)

<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/MemStorage.png" alt="Memory-Storage" width=60%>





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

If you are not sure how much storage you have used in your home directory, feel free to contact us and we can provide you the number. 

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



# **Tufts HPC Cluster Access**



---

> **VPN**
>
> - Off-campus Non-Tufts Network please connect to [Tufts VPN](https://access.tufts.edu/vpn)

> __2FA__
>
> * DUO is needed on Tufts Network (not needed for [**OnDemand**](https://ondemand.pax.tufts.edu), https://ondemand.pax.tufts.edu)
> * DUO is needed when using FastX11 from  [__OnDemand__](https://ondemand.pax.tufts.edu)

> **SSH**
>
> - The SSH protocol aka **Secure Shell** is a method for secure remote login from one computer to another. 

> **X Window System (X11)**
>
> - The X Window System (X11) is an open source, cross platform,  client-server computer software system that provides a **GUI** in a  distributed network environment.

> **Login Hostname**
>
> - login.cluster.tufts.edu = login.pax.tufts.edu
>
> **Cluster New User Guides**
>
> - [Table of Contents](https://tufts.box.com/v/HPC-Table-of-Contents)



## OnDemand

Preferred browser: **Chrome or FireFox**

Go to [**OnDemand**](https://ondemand.pax.tufts.edu), https://ondemand.pax.tufts.edu



<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_login.png" alt="Core-Node" width=70%>

Use your **Tufts UTLN** (Tufts username, lower case!) and **password** to login. 

<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_home.png" alt="Core-Node" width=70%>



<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_menu.png" alt="Core-Node" width=70%>

__`Clusters`__, you can start a shell access to the HPC cluster. 

**`Tufts HPC Shell Access`** = `$ ssh your_utln@login.cluster.tufts.edu `= `$ ssh your_utln@login.pax.tufts.edu `

OR

Use the `>_Open in Terminal` button in `Files` to open a terminal in whichever directory you navigated to.

If you need X11 access through OnDemand to display any GUI applications, please use our [OnDemand](https://ondemand.pax.tufts.edu) **`Clusters`** for this option:

**`Tufts HPC FastX11 Shell Access`** = `$ ssh -XYC your_utln@login.cluster.tufts.edu` (with X11 for GUI applications)

[FastX Web/Desktop Client Setup Instructions](https://tufts.box.com/s/s1vig4km289dzx8qkq4mbhlp4es0oxu1)

OR 

You also have the option to use the `Xfce Terminal` under new  [OnDemand](https://ondemand.pax.tufts.edu) `Interactive Apps` with limited computing resources.



## Shell Access

> **Hostname**: `login.cluster.tufts.edu` or `login.pax.tufts.edu`

### Mac OSX & Linux

- **Terminal** 

  - Shell environment (default: bash):

    `$ ssh your_utln@login.cluster.tufts.edu`

    `$ ssh your_utln@login.cluster.tufts.edu`

    With GUI (Graphical User Interface):

    `$ ssh -XC your_utln@login.cluster.tufts.edu`

    or

    `$ ssh -YC your_utln@login.cluster.tufts.edu`

    X Window System need to be locally installed.

    Now you are on a **Login Node** of the cluster (login-prod-[01-03]) and in your **Home Directory** (~). 

    `$ [your_utln@login-prod-03 ~]`

  * Setting up [SSH keyless access](_https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/_) 
    * Be sure your `~/.ssh` permission is correct! Otherwise, SSH won't work properly.
    * `. ssh` directory: 700 ( drwx------ )
    * public key ( `. pub` file): 644 ( -rw-r--r-- )
    * private key (` id_rsa` ): 600 ( -rw------- )

### Windows

- **[PowerShell](https://learn.microsoft.com/en-us/powershell/scripting/learn/remoting/ssh-remoting-in-powershell-core?view=powershell-7.2)**
- **[WSL - Windows Subsystem for linux](https://learn.microsoft.com/en-us/windows/wsl/install)**
- **[PuTTY](https://www.putty.org/)**  
- **[Cygwin](https://www.cygwin.com/)** 

Need Assistance? Contact us at tts-research@tufts.edu



## File Transfers To/From Tufts HPC Cluster

> **File Transfer Hostname**
>
> - xfer.cluster.tufts.edu = xfer.pax.tufts.edu
>
> **File Transfer Protocol**
>
> - SCP
> - SFTP - Use this for NCBI uploads
> - rsync over SSH
>
> **Globus** * - Coming Soon!
>
> - Tufts HPC Cluster
> - Tufts Box
> - Tufts Sharepoint

### File Transfer Clients

- Windows Only - **[WinSCP](https://winscp.net/eng/index.php)**
- **[FileZilla](https://filezilla-project.org/)**    
- **[Cyberduck](https://cyberduck.io/)**

### OnDemand

- **[OnDemand](https://ondemand.pax.tufts.edu)** (Single file size up to 976MB)

  Only for transfering files size **less than 976MB per file.**

  Go to OnDemand:

  **[https://ondemand.pax.tufts.edu/]( https://ondemand.pax.tufts.edu/)** 

  Under **`Files`**, using the **`Upload`** or **`Download`** buttons to transfer. Make sure you navigate to the destination/source directory on cluster using the **`Go To`** button before transfering files.

  <img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_homedir.png" alt="Core-Node" width=70%>

### Terminal

- **Hostname for file transfer: xfer.cluster.tufts.edu**

  > NOTE:
  >
  > * __Local_Path__ is the path to your files or directory on your local computer
  > * __Cluster_Path__ is the path to your files or directory on the cluster
  >   * Home Directory: */cluster/home/your_utln/your_folder*
  >   * Research Project Storage Space Directory: */cluster/tufts/yourlabname/your_utln/your_folder*

  ***Execute from your local machine terminal.* **

  General Format:

  `$ scp From_Path To_Path`

  `$ rsync From_Path To_Pat`

  ***NOTE: If you are transfering very large files that could take hours to finish, we would suggest using `rsync` as it has ability to restart from where it left if interrupted.***

  **File** Transfer with `scp`or `rsync`:

  * Download from cluster

  `$ scp your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

  `$ rsync your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

  * Upload to cluster

  `$ scp Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

  `$ rsync Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

  **Directory** Transfer with `scp` or `rsync`:

  * Download from cluster

  `$ scp -r your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

  `$ rsync -azP your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

  * Upload to cluster

  `$ scp -r Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

  `$ rsync -azP Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

### Globus - Coming Soon!

[Globus](https://www.globus.org/) is a research cyberinfrastructure, developed and operated as a not-for-profit service by the University of Chicago.

Local Computer: [Globus Connect Personal](https://docs.globus.org/globus-connect-personal/)

Collections: 

- Tufts HPC Cluster Storage (Home and Project)

- Tufts Box

- Tufts Sharepoint

  

## Software, Modules, Packages

### What are modules?

- A tool that **simplify** shell initialization and lets users easily modify their environment during the session with modulefiles
- Each modulefile contains the **information** needed to configure the shell for an application. (PATH, LD_LIBRARY_PATH, CPATH, etc.). Without modules, these environment variables need to be set manually in every new session where the application is needed. 
- Modules are useful in managing **different versions** of applications. 
- Modules can also be bundled into metamodules that will load an entire **set of different applications**. 

### What's NEW?

- Switched from TCL to [Lmod](https://lmod.readthedocs.io/en/latest/010_user.html)
- All previous module commands work as they should + more
- Allows module usage tracking

### How to use modules?

> **Cheat Sheet**
>
> `module av` - check available modules on the MODULEPATH
>
> `module av <software>` - check if a specific software is available as a module
>
> `module spider <keyword>` * - lists all possible modules and not just the modules that can be seen in the current MODULEPATH (such as private modules)
>
> `module --raw show <module_name>` * - printing the modulefile
>
> `module list` - check loaded modules
>
> `module load <software>` - load a specific module
>
> `module unload <software>` - unload a specific module
>
> `module swap <loaded_software> <new_software>` - switch a loaded module for a new one
>
> `module purge` - unload all loaded modules

To **check available modules** installed on the cluster, this may take a few minutes as there are a lot of modules, and be sure to browse the entire list as there are several module file locations:

```
[tutln01@login-prod-01 ~]$ module av
```

Upon login, environment variable **`PATH`** is set for the system to search **executables** along these paths:

```
[tutln01@login-prod-01 ~]$ echo $PATH
/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/home/tutln01/.local/bin
```

For example, I would like to use `gcc` compiler, to check what versions of gcc compiler is available, load the version I would like to use, and use it:

```
[tutln01@login-prod-01 ~]$ module av gcc

----------------------------------------------------------- /opt/shared/Modules/modulefiles-rhel6 ------------------------------------------------------------
gcc/4.7.0 gcc/4.9.2 gcc/5.3.0 gcc/7.3.0

-------------------------------------------------------------- /cluster/tufts/hpc/tools/module ---------------------------------------------------------------
gcc/8.4.0 gcc/9.3.0 gcc/11.2.0
```

Use `module list` to **check loaded modules** in current environment:

```
[tutln01@login-prod-01 ~]$ module load gcc/7.3.0
[tutln01@login-prod-01 ~]$ module list
Currently Loaded Modulefiles:
  1) use.own     2) gcc/7.3.0
```

```
[tutln01@login-prod-01 ~]$ which gcc
/opt/shared/gcc/7.3.0/bin/gcc
[tutln01@login-prod-01 ~]$ echo $PATH
/opt/shared/gcc/7.3.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/home/tutln01/.local/bin
[tutln01@login-prod-01 ~]$ gcc --version
gcc (GCC) 7.3.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

**swap a module for another** (doesn't have to be the same software):

```
[tutln01@login-prod-01 ~]$ module swap gcc/7.3.0 gcc/9.3.0 
[tutln01@login-prod-01 ~]$ module list
Currently Loaded Modulefiles:
  1) use.own     2) gcc/9.3.0
```

 **unload loaded modules**:

```
[tutln01@login-prod-01 ~]$ module unload gcc
[tutln01@login-prod-01 ~]$ echo $PATH
/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/home/tutln01/.local/bin
```

**unload ALL** of the loaded modules in the current environment:

```
[tutln01@login-prod-01 ~]$ module purge
```

### Install Software/Packages

- [R](https://tufts.box.com/s/qximkv5ke2y4k0vbg6m04m6fc6exh88h) (R command line recommanded)
  - Packages need to be <u>reinstalled</u> for each version of R
  - OnDemand 
    - Interactive Apps - RStudio Pax
    - Bioinformatics Apps - RStudio for bioinformatics, RStudio for scRNA-Seq
- [Python](https://tuftsdatalab.github.io/tuftsWorkshops/2024_workshops/2024_bioinformatics301/03_conda/) (Conda env recommanded)
  - Not recommended: anaconda/3 or anaconda/2 (older versions), only use this version when it's absolutely necessary
- Other software compiled from **source**
  - gcc
  - cmake
  - Autotools - automake, autoconf, autogen, .etc
  - ... any **dependencies**, load if available, install if not. Some environment variables may need to be set manually.
  - **Follow instructions** (read it through first)
  - Use **"--prefix="** to install in non-standard locations
  - Modify the environment variables !!! (such as PATH, LD_LIBRARY_PATH, CPATH, .etc)
  - You can make a private module for your locally installed software. Here is [HOW](https://tufts.box.com/s/fewt978g2hmmmhskm1n5dg4bj6xt0cdd)
  - OR you can submit a request to tts-research@tufts.edu for the software to be installed globally on HPC cluster to share with the community. 

---



# Getting Work Done on Tufts HPC Cluster with SLURM

<img src="https://ormstonhouse.com/wp-content/uploads/2019/09/tetris.jpg" alt="teris" width=50%>



**ALL** work **MUST** to be performed on **compute nodes**!

If you see prompt like this `[your_utln@login-prod-01]`, `[your_utln@login-prod-02]`, `[your_utln@login-prod-03]`, **DON'T** run any programs! **Get resource allocation first**!

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
# unload all modules to ensure a clean start.
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
