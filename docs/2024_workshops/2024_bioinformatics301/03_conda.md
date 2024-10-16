# Package Installation Using Conda Environment

Author/Presenter: [Yucheng Zhang](https://github.com/zhan4429), Bioinformatics Engineer, TTS Research Technology        






![Python Environment](https://imgs.xkcd.com/comics/python_environment.png)

## Conda

Conda is an open-source package manager and virtual environment manager for installing packages. 


## conda install is not the first step
`conda install` is a versatile command that simplifies the installation and management of packages.
If you run `conda install` directly without creating and activating a `conda environment` first,  conda will install packages into conda's `root` environment, where you do not have write permissions. 

```
EnvironmentNotWritableError: The current user does not have write permissions to the target environment.
  environment location: /cluster/tufts/software/anaconda/Anaconda3-2020.02
  uid: 33874
  gid: 7593
```
## Conda modules
```
$ module avail anaconda

--------------------- /opt/shared/Modules/modulefiles-rhel6 ------------------------------------------
   anaconda/2  anaconda/3

-------------------- /cluster/tufts/hpc/tools/module --------------------------------------------------
anaconda/bio35  anaconda/2020.02   anaconda/2021.05   anaconda/2021.11   anaconda/2023.07.tuftsai   anaconda/2024.06-py312 (D)
```

`Libmamba` became the default solver for Conda starting with version **23.10.0**, released in November 2023. The new versions of anaconda modules(2024.06-py312) are significantly faster in solving dependencies thanks to the integration of `libmamba`, a highly efficient package management library that outperforms Conda’s traditional solver. **Not recommend to use old anaconda modules**.


!!! note "Anaconda has updated its terms of service([TOS](https://www.anaconda.com/blog/anaconda-commercial-edition-faq)), as below: "

      "We clarified our definition of commercial usage in our Terms of Service in an update on Sept. 30, 2020. The new language states that use by individual hobbyists, students, universities, non-profit organizations, or businesses with **less than 200 employees** is allowed, and all other usage is considered commercial and thus requires a business relationship with Anaconda."

Due to this updated TOS, it's likely we will have to uninstall anaconda from Tufts HPC and other Tufts-owned computers, and migrate to miniforge. **[Miniforge](https://github.com/conda-forge/miniforge) provides a similar repository of packages and software to enable your install and use of Python for your research workflows. Right now, we are waiting to see what Anaconda Inc. decides. In the meantime, you can still use Anaconda, but we recommend new users to use miniforge instead.**

![conda-env-mod-workflow](images/conda_vs_anaconda.png)


| Feature               | Anaconda                                      | Miniconda                                     | Miniforge                                      |
|-----------------------|-----------------------------------------------|-----------------------------------------------|------------------------------------------------|
| **Size**              | ~3 GB                                         | ~400 MB                                       | ~400 MB                                        |
| **Package Manager**    | Conda                                         | Conda                                         | Conda      |
| **Default Channels**   | Anaconda's `defaults`  | Anaconda's `defaults`  | `conda-forge` (community-driven open-source)   |
| **Included Packages**  | Over 300 scientific packages pre-installed  | Minimal installation (essential tools only)   | Minimal installation (essential tools only)    |
| **Target Audience**    | Data scientists needing full environment      | Developers who want a lightweight version     | Developers focused on open-source solutions    |
| **Ease of Use**        | Easy with pre-installed packages              | Requires manual installation of most packages | Similar to Miniconda, but optimized for `mamba`|
| **Solver**             | Mamba                                         | Mamba                                         | Mamba   |


## miniforge modules
```
$ module avail miniforge
 ----------------------/cluster/tufts/hpc/tools/module-----------------------
   miniforge/24.3.0-py310    miniforge/24.7.1-py312 (D)
```

## conda channels

- Conda channels are the locations where packages are stored. 

- Conda search or download packages from channels. 

- For anaconda and miniconda, the default set of channels is called **defaults**. 

- To install a package that is not in **defaults**, you need to tell conda which channel contains the package. 

```
$ conda install -c pytorch pytorch 
$ conda install -c conda-forge r-base
$ conda install -c bioconda samtools
```
[Bioconda](https://bioconda.github.io/#) is a community-enabled conda channel of thousands of bioinformatics packages.
![bioconda](https://bioconda.github.io/_images/bioconda.png)
### .condarc

By default, Conda stores packages in your $HOME directory. Since each user has a storage limit of 30GB in $HOME, we recommend avoiding package installations there to prevent exceeding your quota. If you are a member of the XXXXlab group on the cluster, we suggest using the group’s research storage space for this purpose.

To change the default location to your group's project directory, use a text editor to create a hidden file called `.condarc`. The path to that file should be `$HOME/.condarc`. 
```
$ cd 
$ ls -a # -a will reveal hidden files starting with .
$ touch .condarc # If you haven't created .condarc before, you can create it using the touch command
```
Create two directories in your group research storage space (one for storing the `envs`, one for storing the `pkgs`, for example: condaenv, condapkg)

```
$ mkdir /cluster/tufts/XXXXlab/$USER/condaenv/
$ mkdir /cluster/tufts/XXXXlab/$USER/condapkg/
```

Now add the following 4 lines to the `.condarc` file in your home directory (modify according to your real path to the directories):

```
envs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condaenv/
pkgs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condapkg/
```

**OR** you can do so from command line with the following commands:

```
$ module load miniforge # or module load anaconda
$ conda config --append envs_dirs /cluster/tufts/XXXXlab/$USER/condaenv/
$ conda config --append pkgs_dirs /cluster/tufts/XXXXlab/$USER/condapkg/
```

Another common use of `.condarc` is setting the channels that will be used to search for packages to install. For example, to add the `bioconda` and `conda-forge` channels, you can run :

```
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

**According to Anaconda's updated terms of service, you cannot download and install packages from Anaconda's `defaults` channel without a paid license**. You'd better remove this `defaults` channel from your `.condarc` by:
```
$ conda config --remove channels defaults 
```

After this, your `.condarc`file should look like this

```
$ cat ~/.condarc
envs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condaenv/
pkgs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condapkg/
channels:
  - bioconda
  - conda-forge
```

## Conda environment

With Conda, you can create and use environments that have different versions of python and other packages. Switching between conda environments is called activating the environment. 

Conda environments are a useful way to install packages with specific versions and dependencies. 

**It is recommended to always use conda environments**. By default, conda environments are located in `./conda/envs` of $HOME. If you have created `.condarc` following the above steps, conda environment will be stored in other directories.

### conda create -n
To create a new Conda environment named `myenv`, you can use the following command:
```
$ conda create -n myenv
```
You can also use this way:
```shell-session
$ conda create --name myenv
```
This will create a new environment called `myenv` without installing any specific packages. If you want to create the environment with specific packages (e.g., Python, NumPy, etc.), you can specify them as follows:
```
$ conda create -n myenv python=3.11 numpy
```

### conda create -p

You can control where a conda environment lives by providing a path to a target directory when creating the environment. For example, you want to create a conda environment in your lab's project folder instead of $HOME.

```
$ conda create --prefix /cluster/tufts/zhanglab/condaenv/myenv2
or 
$ conda create -p /cluster/tufts/zhanglab/condaenv/myenv2
```

### source activate and deactivate

```
$ source activate myenv1
$ source deactivate
$ source activate /cluster/tufts/zhanglab/condaenv/myenv2
$ source deactivate
```


!!! note "**conda activate** is supported"

  Since `anaconda/2024.06-py312`, our anaconda and miniforge modules will support both `conda activate` and `source activate`: 

  ```
  $ conda activate myenv1
  $ conda deactivate
  ```
  You still need to use **source activate** in old anaconda modules. 

### Do not run `conda init`
```
CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'.
To initialize your shell, run

    $ conda init <SHELL_NAME>

Currently supported shells are:
  - bash
  - fish
  - tcsh
  - xonsh
  - zsh
  - powershell

See 'conda init --help' for more information and options.

IMPORTANT: You may need to close and restart your shell after running 'conda init'.
```

You can run the below command to avoid `conda init`
```
$ eval "$(conda shell.bash hook)"
```

We do not recommend using `conda init`, it permanently alters your `.bashrc` so that the `base` environment is always activated in new shells. 

```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/cluster/tufts/software/anaconda/Anaconda3-2020.02/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/cluster/tufts/software/anaconda/Anaconda3-2020.02/etc/profile.d/conda.sh" ]; then
        . "/cluster/tufts/software/anaconda/Anaconda3-2020.02/etc/profile.d/conda.sh"
    else
        export PATH="/cluster/tufts/software/anaconda/Anaconda3-2020.02/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

```
(base) $ which conda
/cluster/tufts/software/anaconda/Anaconda3-2020.02/bin/conda
```

If you ran `conda init` previously, it is recommended to delete the `conda init` section from your `.bashrc`.

## Pip

pip stands for "pip installs packages". It is a package manager for Python packages. pip installs packages that are hosted on the Python Package Index or [PyPI](https://pypi.org/).

### conda install vs pip install

| Feature                  | `conda install`                                              | `pip install`                                                |
| ------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **Package Management**   | Manages both python packages and non-python libraries, dependencies, and environments. | Manages only python packages and dependencies.               |
| **Environment Creation** | Can create and manage isolated environments with both python and non-python dependencies. | Requires `venv` or similar tools for environment creation.   |
| **Dependency Handling**  | Handles complex dependencies more robustly by using curated, pre-built binaries. | Resolves dependencies but might lead to version conflicts (dependency hell). |

### conda install or pip install?
If your package or versions exist supports both `conda install` and `pip install`.  You should almost always favor conda over pip. Because conda uses a more robust dependency resolution mechanism. The packages available in the Conda repositories are pre-built and tested to work together, reducing the likelihood of version conflicts (sometimes called **dependency hell**).

### pip venv
Similar to conda environment, it is also recommended to create isolated python environments to avoid/reduce conflicts. Python3 provides a `venv` module that can be used to lightweight virtual environment. 

A virtual environment can be created by running "`python3 -m venv <DIR>`", where `<DIR>` is the directory of the created environment. 

**`venv` is very useful to manage difffernt python projects. However, we would like to reduce users' confusion, and hope users to stick to conda environment for both python and non-python applications**.

### pip install --user is not recommended
In user mode installation, packages are installed under the user’s home directory (in `$HOME/.local/` and carry the risk of being *contaminated* over time. For example, during minor upgrade of the python version older packages may no longer work. Similarly, if users install multiple packages with conflicting dependencies over time, some packages would break. This can be avoided by creating a virtual environment (such as an conda environment), which provides an isolated location for installing a group of packages. 

### conda create --name myenv python=3.xx
- Creates a new conda environment named myenv with python version 3.xx pre-installed.
- **Python modules installed using pip within this environment will be installed in the myenv environment's specific directory.** This ensures isolation and prevents conflicts with other environments or the base conda environment.

### conda create --name myenv: **This is not recommended**
Creates a new conda environment named myenv without specifying a python version.

**By default, this environment will likely use the base Python installation on your system.** If you install python modules using pip, their location will depend on the default python installation's configuration (usually it is user’s home directory: `$HOME/.local/`), which can lead to potential conflicts and inconsistencies. 


## Using Conda environments in Jupyter

If you would like to use JupyterNotebook or JupyterLab from OnDemand, you can follow the instructions below and **run your conda env as a kernel in Jupyter**.

- Make sure with python 3.7+ and make sure you load cluster’s anaconda module (this only works with py3.7+)
- Create conda environment with its own python

  ```
  module load anaconda/2024.06-py312
  conda create --name myenv python=3.xx
  ```

- Activate the conda environment

  ```
  $ conda activate myenv
  ```

- Install the ipykernel python module 

  ```
  $ pip install ipykernel
  ```

- Create the kernel file

  ```
  $ python -m ipykernel install --user --name=mykernelname ## give a meaningful kernel name
  ```

- Restart Jupyter from OnDemand

- :)

## Recommended workflow
```
$ module load miniforge/24.7.1-py312  
$ conda config --append envs_dirs /cluster/tufts/XXXXlab/$USER/condaenv/ # Once
$ conda config --append pkgs_dirs /cluster/tufts/XXXXlab/$USER/condapkg/ $ Once
$ conda config --remove channels defaults  # Once
$ conda config --add channels bioconda   # Once
$ conda config --add channels conda-forge # Once
$ conda create -n myenv python=3.12
$ conda activate myenv
$ conda install xxx
$ pip install xxx
$ conda deactivate # After your job/analysis is complete
```

## Caveats
- Do not run `conda init`
- Do not install python modules with `pip install --user`
- Watch for disk usage in your $HOME
  ```
   $ module load hpctools
   $ hpctools
   Please select from the following options:
  1. Checking Free Resources On Each Node in Given Partition(s)
  2. Checking Free GPU Resources On Each Node in Given Partition(s)
  3. Checking Past Completed Jobs in Given Time Period
  4. Checking Active Job informantion
  5. Checking Project Space Storage Quota Informantion
  6. Checking Any Directory Storage Usage Informantion
  Press q to quit
  ## You can select 6 to check your $HOME disk usage
  ```

[Previous: R](02_Rpackage.md)

[Next: conda-env-mod](04_conda-env-mod.md)