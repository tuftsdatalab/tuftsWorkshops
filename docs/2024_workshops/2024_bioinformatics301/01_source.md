# Installation from source codes

![permission](images/permission.png)

## make
GNU Make is a program often used for compiling software. It uses a plain text file named **makefile** or **Makefile**.

`make` and `Makefile` are also widely used in building reproducible workflows. This [ariticle](http://www.bioinformaticszen.com/post/makefiles/) is a good introduction.

### Steps

1. Unpack the source code archive. 

2. **Configure** the package. ## Some packages do not have the **configure** file

3. Run **make** to build the programs. 

4. Run **make install** to install the package. # Optional

   ❌ Do not run **sudo** **make install**

By default, **make install** will install applications into `/usr/local`, but regular users do not have permission to write into `/usr/lobcal`. 

The best way is to install applications into your home directory or your group's shared directory by passing the option `--prefix=TargetDirName` to `./configure`. 


### bwa
[BWA](https://bio-bwa.sourceforge.net/) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
#### Installation guide from [README](https://github.com/lh3/bwa)
```
git clone https://github.com/lh3/bwa.git
cd bwa; make
./bwa index ref.fa
./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz
```

```
## Load the compiler 
$ module av gcc
$ module load gcc/11.2.0 ## Recommend to load the newest version of gcc
$ cd $HOME
$ mkdir apps 
$ cd apps
$ git clone https://github.com/lh3/bwa.git
$ cd bwa; make
```

```
# Add the bwa directory to $PATH
$ export PATH=$PATH:$HOME/apps/bwa
# bwa is ready to use
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

### hmmer
[HMMER](http://hmmer.org) is used for searching sequence databases for sequence homologs, and for making sequence alignments. It implements methods using probabilistic models called profile hidden Markov models (profile HMMs).
#### Installation guide from [README](https://github.com/EddyRivasLab/hmmer)
```
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.4
./configure --prefix /your/install/path   # replace /your/install/path with what you want, obv 
make
make check                                # optional: run automated tests
make install                              # optional: install HMMER programs, man pages
```

```
$ cd $HOME
$ cd apps
# Download source code
$ wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
# Unpack source code archive 
$ tar -xvf hmmer-3.4.tar.gz
$ cd hmmer-3.4
# compile the code
./configure --prefix=$HOME/apps 
make
make check
make install
```

```
# Add the bwa directory to $PATH
$ export PATH=$PATH:$HOME/apps/bin
# hmmer is ready to use
$ hmmsearch -h
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# Copyright (C) 2023 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: hmmsearch [options] <hmmfile> <seqdb>
```

## cmake
**CMake** is an open-source, cross-platform family of tools designed to build, test, and package software. It controls the software compilation process by generating native build scripts (like Makefiles or project files) for a wide variety of platforms and compilers.

Installation of some bioinformatics applications requires both **make** and **cmake**.

### RegTools
[RegTools](https://github.com/griffithlab/regtools) integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context.

#### Installation [guide](https://regtools.readthedocs.io/en/latest/) from the developer 

```
git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
```

#### Installation on HPC

```
 $ module avail cmake

--------------------- /opt/shared/Modules/modulefiles-rhel6 ------------------------------
   cmake/2.8    cmake/2.8.11.2    cmake/3.2.1    cmake/3.4.3

--------------------- /cluster/tufts/hpc/tools/module ------------------------------------
   cmake/3.18    cmake/3.23_gui (D)

```

```
$ module load gcc/11.2.0
$ module load cmake/3.23_gui ## Recommand to use the latest version of cmake
$ cd $HOME/apps
$ git clone https://github.com/griffithlab/regtools
$ cd regtools/
$ mkdir build
$ cd build/
$ cmake ..
$ make
```

```
$ export PATH=$PATH:$HOME/apps/regtools/build
```

### DCMAKE_INSTALL_PREFIX
Some applications' installation also has install stage, which will have `make intall` as the last step. For these installations, we have to include `-DCMAKE_INSTALL_PREFIX` in the `cmake ..` step. Below are the common steps for such installations:

```
module load gcc
module load cmake
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make
make install
```

