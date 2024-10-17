# Installation from source codes

![permission](images/permission.png)



## Understanding Make, CMake, and Software Installation on HPC

This guide will help you understand how **GNU Make** and **CMake** are used in software installation, particularly for bioinformatics applications. These tools are essential for managing the building and compiling of programs from source code. Additionally, you will learn about specific software installation instructions for bioinformatics tools like BWA, HMMER, and RegTools on a High-Performance Computing (HPC) environment.



## What is Make?

[GNU Make](https://www.gnu.org/software/make/) is a program often used for compiling software. It uses a plain text file named **makefile** or **Makefile**, which lists each of the non-source files and how to compute it from other files.

`make` and `Makefile` are also widely used in building reproducible workflows. This [ariticle](http://www.bioinformaticszen.com/post/makefiles/) is a good introduction.



!!! note "What is complier? "

    A compiler is a program that translates source code written in a high-level programming language (such as C, C++, or Java) into machine code or bytecode that a computer's processor can understand and execute. This process involves several steps, including lexical analysis, syntax analysis, optimization, and code generation.

### Steps for software installation using make

1. Unpack the source code archive. 

2. **Configure** the package.                                # Some packages do not have the **configure** file

3. Run **make** to build the programs. 

4. Run **make install** to install the package.  # Optional

   ❌ Do not run **sudo** **make install**

**Tip:** By default, **make install** will install applications into `/usr/local`, but regular users do not have permission to write into `/usr/local`. 

The best way is to install applications into your home directory or your group's shared directory by passing the option `--prefix=TargetDirName` to `./configure`. 



!!! note "What is configuration?"

    **Configuration** refers to the arrangement or setup of various components and settings within a system, software, or device to achieve a specific behavior or function. It involves specifying options and parameters that control how the system or software operates.



!!! note "`make` and `make install`"

    **make**: Compiles the source code and creates binaries, typically in the current directory.              
    **make install**: Installs the compiled program into system-wide directories, so it can be run from anywhere on the system. This step usually follows after make. 



### Installing bwa using make

[BWA](https://bio-bwa.sourceforge.net/) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. **Installation guide from the developer can be found [HERE](https://github.com/lh3/bwa).** 

#### Installation on Tufts HPC 

Step1: Go to the folder where you would like to install the tools. Ex: We create `apps` folder under $HOME to install the tools.

```
cd $HOME
mkdir apps
cd apps
```

Step2: Load the compiler. **GCC** (GNU Compiler Collection) is required for installing **BWA**. Since BWA is written in C, you need a C compiler like GCC to compile the source code. When you run the `make` command to compile BWA, it invokes the GCC compiler to build the binaries from the source code. 

```
module av gcc          # check which version of gcc is available. 
module load gcc/11.2.0 # Recommend to load the newest version of gcc
```

Step3: Clone the bwa repository

```git clone https://github.com/lh3/bwa.git```

Step4: Configure and build

```
cd bwa
make
```

Step5: Add the build directory to your PATH

```export PATH=$PATH:$HOME/apps/bwa```

Step6: Now bwa is read to use

```
$ bwa
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.18-r1243-dirty
Contact: Heng Li <hli@ds.dfci.harvard.edu>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
         pemerge       merge overlapping paired ends (EXPERIMENTAL)
         aln           gapped/ungapped alignment
         samse         generate alignment (single ended)
         sampe         generate alignment (paired ended)
         bwasw         BWA-SW for long queries (DEPRECATED)

         shm           manage indices in shared memory
         fa2pac        convert FASTA to PAC format
         pac2bwt       generate BWT from PAC
         pac2bwtgen    alternative algorithm for generating BWT
         bwtupdate     update .bwt to the new format
         bwt2sa        generate SA from BWT and Occ

Note: To use BWA, you need to first index the genome with `bwa index'.
      There are three alignment algorithms in BWA: `mem', `bwasw', and
      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
      first. Please `man ./bwa.1' for the manual.
```



### Installing hmmer using make
[HMMER](http://hmmer.org) is used for searching sequence databases for sequence homologs, and for making sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs). **Installation guide from the developer can be found [HERE](https://github.com/EddyRivasLab/hmmer).**

#### Installation on Tufts HPC 

Make sure you are on compute mode by running the command below

```
srun -p interactive -n 1 --time=02:00:00 --mem 4g --pty bash
```



Step1: Go to the folder where you would like to install the tools. Ex: We create `apps` folder under $HOME to install the tools.

```
cd $HOME
mkdir apps
cd apps
```

Step2: Download and unpack the source code

```
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar -xvf hmmer-3.4.tar.gz
cd hmmer-3.4
```

Step3: Configure and build the software

```
./configure --prefix=$HOME/apps   # replace /your/install/path with what you want
make
make check                        # optional: run automated tests
make install                      # Install HMMER programs, a bin folder will be created under $HOME/apps 
```

Step4: Add HMMER to your PATH

```
export PATH=$PATH:$HOME/apps/bin
```

Step5: Now HMMER is read to use

```
$ hmmsearch -h
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# Copyright (C) 2023 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: hmmsearch [options] <hmmfile> <seqdb>
```



**Tip**: Check what has been created after each step to understand how the software is built and installed. 

```
# After "tar -xvf hmmer-3.4.tar.gz"
ls -lhtr $HOME/apps/hmmer-3.4
# After "make"
ls -lhtr $HOME/apps/hmmer-3.4
ls -lhtr $HOME/apps/
# After "make install"
ls -lhtr $HOME/apps/hmmer-3.4
ls -lhtr $HOME/apps/
```





## What is CMake?
**CMake** is an open-source, cross-platform family of tools designed to build, test, and package software. It controls the software compilation process by generating native build scripts (like Makefiles or project files) for a wide variety of platforms and compilers.

Installation of some bioinformatics applications requires both **make** and **cmake**.

### Installing RegTools Using Make and CMake
[RegTools](https://github.com/griffithlab/regtools) integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context. Installation guide from the developer can be found [HERE](https://regtools.readthedocs.io/en/latest/)

#### Installation on Tufts HPC

```
module avail cmake         # Always use the latest version

--------------------- /opt/shared/Modules/modulefiles-rhel6 ------------------------------
   cmake/2.8    cmake/2.8.11.2    cmake/3.2.1    cmake/3.4.3

--------------------- /cluster/tufts/hpc/tools/module ------------------------------------
   cmake/3.18    cmake/3.23_gui (D)

```

```
# Load modules
module load gcc/11.2.0
module load cmake/3.23_gui ## Recommand to use the latest version of cmake

# Go to the folder where you would like to install tools
cd $HOME/apps

# Clone the github repo
git clone https://github.com/griffithlab/regtools

# Create a folder called `build`. 
# cmake is designed to work well with out-of-source builds, where the build directory contains all the generated files (e.g., Makefiles, binaries, configuration files). 
# It's a common practice to create build folder

cd regtools/
mkdir build
cd build/

# Run cmake first and then make
cmake ..
make
```

```
export PATH=$PATH:$HOME/apps/regtools/build

# The tool is now successfully installed. 
regtools --help
```

### DCMAKE_INSTALL_PREFIX

Some applications' installation also has install stage, which will have `make intall` as the last step. For these installations, we have to include `-DCMAKE_INSTALL_PREFIX` in the `cmake ..` step. Below are the common steps for such installations:

```
module load gcc/11.2.0
module load cmake/3.23_gui
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make
make install
```

## Summary

In this guide, we covered the use of **Make** and **CMake** in software installation, particularly in the context of HPC environments. We also provided specific examples for installing bioinformatics tools such as **BWA**, **HMMER**, and **RegTools**. These examples emphasize the importance of configuring installation paths to avoid system permission issues and ensure software is installed in user-accessible directories.



[Previous: Intro](00_introduction.md)

[Next: R](02_Rpackage.md)