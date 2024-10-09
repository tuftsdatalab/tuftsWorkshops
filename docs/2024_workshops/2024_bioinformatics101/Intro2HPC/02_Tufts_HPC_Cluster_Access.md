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
> - The X Window System (X11) is an open source, cross platform,  client-server computer software system that provides a Graphical User Interface (**GUI**) in a  distributed network environment.

> **Login Hostname**
>
> - login.cluster.tufts.edu = login.pax.tufts.edu
>
> **Cluster New User Guides**
>
> - [Table of Contents](https://tufts.box.com/v/HPC-Table-of-Contents) - New documentation site coming soon! Look out for our emails!



## OnDemand

> Preferred browser: **Chrome or FireFox**

Go to [**OnDemand**](https://ondemand.pax.tufts.edu), https://ondemand.pax.tufts.edu



<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_login.png" alt="Core-Node" width=70%>

Use your **Tufts UTLN** (Tufts username, lower case!) and **password** to login. 

<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_home.png" alt="Core-Node" width=70%>



<img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_menu.png" alt="Core-Node" width=70%>

__`Clusters`__
- Start a shell access to the HPC cluster.

  **`Tufts HPC Shell Access`** = `$ ssh your_utln@login.cluster.tufts.edu `= `$ ssh your_utln@login.pax.tufts.edu `

  > Note: the `>_Open in Terminal` button in `Files` also opens a terminal in whichever directory you navigated to.

- If you need X11 access through OnDemand to display any GUI applications, use [OnDemand](https://ondemand.pax.tufts.edu)
  
  **`Tufts HPC FastX11 Shell Access`** = `$ ssh -YC your_utln@login.cluster.tufts.edu` (with X11 for GUI applications)

  [FastX Web/Desktop Client Setup Instructions](https://tufts.box.com/s/s1vig4km289dzx8qkq4mbhlp4es0oxu1)


## Shell Access

> **Hostname**: `login.cluster.tufts.edu` or `login.pax.tufts.edu`

### Mac OSX & Linux

**Terminal** 

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

  - Setting up [SSH keyless access](_https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/_) 
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

