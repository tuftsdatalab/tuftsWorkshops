# Introduction To Unix

# Intro to Command Line
=======================================

* TOC
{:toc}

This short workshop provides some basic training on bash and shell scripting on the command line on the Linux-based Tufts HPC cluster.

This course is not meant to be comprehensive, but provides some insights into how the command line works as well as some strategic resources for studying and understanding command line on the HPC cluster.

#### Helpful Tip: Some Vocabulary
===================================

"Command line" is a more general term to indicate that you are using text commands on a terminal (linux bash shell or similar). Command line differs from "Graphical User Interface (GUI)" because all commands are texts instead of drag-and-drop or interactive formats such as the Windows or Mac Operating Sytems provide.

"HPC" stands for High Performance Computing, "cluster" refers to a shared computer resource to enable more powerful computation than regularly available on an individual machine.

"Linux" can refer to any of the free open source version of "Unix" from AT&T Bell labs who pioneered the language in 1965. There are a number of Linux operating systems installed on HPC clusters (Ubuntu, Debian, RedHat Enterprise License (RHEL), CentOs, Fedora, etc.) Each of these systems have slight differences that may impact the commands demoed here. Tufts University Research Cluster is currently using RHEL7.

"Bash" is one type of languages used in a "shell", the text interface on the Linux system. This lesson introduces a few objectives to help users understand how to use bash commands on the Linux RHEL shell of our HPC. Other shell languages have slight differences that affect how commands are run (e.g. new MacOSX ship with "zsh" as the default shell language on their installed terminal programs).


## Learning Objectives
-----------------------

- What is the shell?
- How do you access it?
- How do you use it and what is it good for?

  * Running commands
  * File Directory Structure
  * Manipulating files
  * Simple Bash Scripts

## What is the shell?
------------------

The **shell** is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

There are many reasons to learn about the shell.  A few specific ones:

* For most bioinformatics tools, you have to use the shell. There is no
  graphical interface. If you want to work in metagenomics or genomics you're
  going to need to use the shell.

* The shell gives you **power**. The command line gives you the power to
  do your work more efficiently and more quickly. Shell allows users to automate repetitive tasks.

* To use remote computers or cloud computing, you need to use the shell.


### Knowing Shell Increases Speed and Efficiency Through Automation

The most important reason to learn the shell is to learn about
**automation**.  Any time you find yourself doing roughly the same
computational task more than few times, it may be worth automating it;
the shell is often the best way to automate anything to do with files.

In this lesson, we're going to go through how to access Unix/Linux and some of the basic
shell commands. We will finish with a demonstration of how to run programs interactively as well by submitting a job to SLURM (https://it.tufts.edu/sites/default/files/uploaded-files/2020-03/QuickStart%20for%20Slurm.pdf). Slurm is a scalable cluster management and job scheduling system for Linux clusters. Other job scheduling systems you may be familiar with from other universities are "PBS" and "SGE_Batch".

### Where to learn shell commands
-----------------------------------------

The challenge with bash for the command line is that it's not particularly simple - it's a
power tool, with its own deep internal logic with lots of details.

Practice is the best way to learn, but here are some helpful shell command resources:

* [Fun With Unix Cheat Sheet](https://files.fosswire.com/2007/08/fwunixref.pdf)
* [Shell Cheatsheet - Software Carpentry](https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md)
* [Explain shell](http://explainshell.com) - a web site where you can see what the different
components of a shell command are doing.

### How to Access the Shell at Tufts
------------------------------------

[Log In through OnDemand](https://ondemand.pax.tufts.edu/pun/sys/dashboard){:target="_blank" rel="noopener"}

<img width="692" alt="Ondemand_Shell" src="https://user-images.githubusercontent.com/8632603/179539946-5d4fa52d-95ae-4215-ab16-24c912879aeb.png">

Mac
---

On Mac the shell is available through Terminal  
Applications -> Utilities -> Terminal  
Go ahead and drag the Terminal application to your Dock for easy access.
Note: newer versions of MacOSX ship with "zsh" as the default shell language in their terminal. It is possible to change the preference to "bash". However, if you are only using the terminal to log into the Tufts cluster, you don't necessarily need to do this, because you will be using "bash" once you are on the cluster. "zsh" would only impact scripts and commands run locally on your own machine.

Windows
-------

For Windows, an easy one to install and use right away is  gitbash.  
Download and install [gitbash](https://gitforwindows.org/)
Open up the program.

Other options: 

* [Microsoft Window Terminal](https://docs.microsoft.com/en-us/windows/terminal/install)
* [Conemu](https://conemu.github.io/)

Linux
-----

You probably already know how to find the shell prompt.


[Return to TOC](#intro-to-command-line)

---------------------------

## Starting with the Shell
---------------------------

We will spend most of our time learning about the basics of the shell by manipulating some experimental data that we download from the internet.

Open up the shell through a terminal (OnDemand or on your laptop) and type the command::

```
whoami
```

and then hit ENTER 

(This is a good question for Mondays ....)

When you are on the Tufts cluster, this will return your username according to the cluster. This username is attached to you wherever you are in the cluster and creates a home where your files can be kept, regardless of which machine you are on in the cluster. [If you are on your laptop or personal computer, the answer to this may be different before you log in.]

---
---
#### For Attendees Using Terminal Programs to Access the Cluster (instead of the Web Browser "OnDemand")

If you are using a terminal on your home machine to connect to the tufts cluster, you will first need to log in by sending a simple command. **Ignore this if you are using the web browser login tool.**

Replace "username01" with your tufts username.

```
ssh username01@login.pax.tufts.edu
```

Your username will have been created when your account was set up. If you do not have a cluster account, you can still follow this tutorial from your laptop or personal computer, except that the file structure will be different from what is described.

The login will ask you for your Tufts password.

If you are not on the Tufts network, you will need to set up the Tufts VPN (Virtual Private Network) before logging in:

https://it.tufts.edu/guides/vpn-virtual-private-network/anyconnect-desktop-application

---
---


## Running Commands
--------------------

When you login, you will see some letters and characters at the beginning of the line.

This is called the 'command line prompt.'

```
[username01@login-prod-02 ~]$
```

This information helps orient you to who you are (`username01`) and which computer you are currently on(`login-prod-02`). In this case, it is a server that is intended only for "login", no big programs should be run from this computer, but it is fine for practicing a few bash commands.

#### Helpful Tip
===================

The name of the computer you are on is important informatiom when troubleshooting the cluster. `login` machines will reject large commands and output an error. Make sure to switch machines before running jobs. This is explained in the HPC portion of the lesson.

----------------------


The `$` at the end of the line is where you start typing your commands. The `$` (on the Mac it is a `%`) is not part of the command.

The outputs from commands will not have that piece of information or `$` at the beginning of the line.

Let's try some simple commands.

Much like text shortcuts, shell commands often use abbreviations to get their point across.

For example, the command *pwd* is short for "pass working directory."

Now type the command

```
pwd
```

You should see something similar to this:

```
/cluster/home/username01/
```

Try this command

```
ls
```

It may be empty for the moment, or it may not if this is not your first time using the shell.


### Takeaways
=================

`pwd` and `ls` are examples of commands - programs you run at the shell
prompt that do stuff. 

* `pwd` stands for 'print working directory', while
* `ls` stands for 'list files'. 

It is similar to the abbreviations used in texting, it takes less time to get the point across (lol, tbh, imho, afaik, ftw -- you're saying them outloud in your head, right now, correct?)


## Parameters for Bash Commands
--------------------------------

Many bash commands have special **parameters**, sometimes referred to as **flags** that open up a lot more possibilities.

Let's start by going to your home directory (you choose the command)


As you start using bash more and more, you will find a mix of files and directories/folders. If we want to know which is which, we can type::

```
ls -F
```

Anything with a "/" after it is a directory.  Things with a `*` after
them are programs.  It there's nothing there it's an otherwise
unremarkable file (e.g. a data file).

Depending on which terminal you are using, some of the file types may have different colors. 

In our ondemand shell:

Files are white
Directories are blue
Programs are green
Compressed files are red (e.g. files that end in .zip or .gzip or .tar)


You can also use the command::

    ls -l

to see whether items in a directory are files or directories. `ls -l`
gives a lot more information too, such as the size of the file.


It also shows the permissions of who can read, write or execute a file.


```

drwxrwx--- 2 username05 username05     4096 Jul 18 09:57 JulyWorkshop


```

The first 10 letters in this line indicates the permission settings.


<img width="523" alt="File_Permissions" src="https://user-images.githubusercontent.com/8632603/179539739-75f4edf9-5f5d-4de9-b20c-97abc7869be6.png">


#### Helpful Tip
================

There are an overwhelming number of possibilities with some of these shell commands, so knowing how to find help on demand is important.

For example, `ls` has a lot of flags that can be used.

```
ls --help
```

This outputs a list of all the ways that `ls` can be altered to find information about your files.

Parameters can be added together in some cases.

```
ls -ltr
```

This command strings together three flags.

`ls -l` is list with details
`ls -t` is sort the list by creation time
`ls -r` is sort the list in reverse

For very full directories, this is helpful because it outputs the most recent set of files as the last in the list.


Another way to get help is to use the `man` command. Not every unix installation has this installed, but the Tufts cluster does.

```
man ls
```

This opens up the manual on the `ls` command. It spells out the meaning of all the parameters in detail.

Most common bash commands have a `man` page that explains it (I wish they had this for emojis....).

Many programs have a help function built in, try adding `--help` or `-h` to see if some helpful information pops up. Sometimes just running the command without any arguments or parameters leads to some usage information or describes the correct command to get help.

For example, if I want to understand the command `tr`

```
tr -h
```

shell outputs

```
tr: invalid option -- 'h'
Try 'tr --help' for more information
```


[Return to TOC](#intro-to-command-line)

## Navigating in the Shell
===============================

We are going to make a place to work for this workshop.

The following command makes a new directory.

```
mkdir JulyWorkshop

```
----------------------

#### Helpful Tip
===================

When nameing files and directories, avoid spaces and special characters except underscores ("_").

**Spelling** and **Capitalization** are literal in unix, be careful when making and using files to be consistent in your process. This will make it easier to find files later.

---------------------

You can check that the new directory was created by repeating the list command.

```
ls
```

A directory is like a desk drawer. We create them to store files that relate to each other mostly.

When creating directories and filenames it is helpful to put some information about the project and the date of activity.


<img width="786" alt="File_Folder_Structure" src="https://user-images.githubusercontent.com/8632603/179539866-ecd6e880-f468-4151-bbaa-149f52c328b4.png">

## Absolute and Relative Paths
-------------------------------

Let's go into our directory and look around.

Another command you'll find yourself using a lot is `cd`, which stands
for 'change directory'.  Try typing::

```
cd JulyWorkshop

```
and then

```
pwd
```

You should now see something like this:

```
/cluster/home/username01/JulyWorkshop
```

This is an example of an **Absolute Path**.

It gives an address for where you are located on the cluster, much like a postal address that defines where you are in several layers (e.g. /country/state/city/street/specific_house.

<img width="821" alt="tufts_root_path" src="https://user-images.githubusercontent.com/8632603/179759502-549b38b0-4957-4105-aee3-8bca4271bf7b.png">

If you want to go back to the directory that is in the level above our current file, another common shortcut used in bahs is `..`.


```
cd ..
```

`..` is a reference to a **RELATIVE PATH**

```
pwd
```

You should be back in your home directory.

```
/cluster/home/username01/
```

If you want to go back to the directory that you just left, type this command.

```
cd -
```
Then find your location.

```
pwd
```

You should be back in the directory you came from.

```
/cluster/home/username01/JulyWorkshop
```

A **RELATIVE PATH* means that the command only works from the relative location that you are in.

`cd ..` and `cd -` are examples of relative path commands.

This can get confusing if you are moving around a lot in your directories or sending commands to SLURM, so the alternative method to navigating around the cluster is using an **ABSOLUTE PATH**.


Note: If you ever type `cd` without a word behind it, it will send you back to your home directory.

Your home directory is not all the way back at the root, it is set within the cluster as `/cluster/home/username01/`.

You can make sure that you are in the right directory by using the command `cd` with the absolute path.

```
cd /cluster/home/username01/JulyWorkshop
```

#### Helpful Tip
=================
Many commands in bash can be used with the ABSOLUTE PATH.

```
ls /cluster/home/username01/JulyWorkshop
```

Using an absolute path to find files in a directory is helpful for checking for outputs from SLURM jobs when they are running.


## Creating and Manipulating Files
----------------------------------

Let's make a file here using a common command "echo" to start creating our file structure.

```
echo "Hello World " > helloworld.txt

```

The `>` in this command tells the command to place the output into the place it is pointing. In this case, it creates the file `helloworld.txt` and puts the phrase `Hello World` into the file. 

#### Helpful Tip
===================

Be careful with redirect.

When using `>` to redirect content into a file, if the filename already exists, it will **overwrite** the file. This means that the original file is gone, and there is no undo in shell.

If you want to add to a file (for example if you are running the same command on several files and extracting a piece of information that you want to put together at the end) you can use another form of redirect `>>`. Using the double redirect will **add** to the file instead of overwriting it.

Which one is used depends on your process. If you are only running a command once, or have an intermediate file in a process that does not need to be retained at the end, then `>` is okay to use.

-------------------------

### Reading File Contents
---------------------------

Let's look inside the file. We have several methods of viewing the content of files that we have created.

A helpful command is `cat`.

```
cat helloworld.txt
```

"cat" will open the entire file, so this is not the best command for long files.

In that case "head" is a good option. Head pulls the top ten lines of the file and prints them to the screen.

```
head helloworld.txt
```

It does not look any different from cat in this case because there is only one line in the file.

A third way to check file contents is by using a program called "less" (or "more").

"less" will open the file interactively, then you can scroll through it and when you are done, push "q" on your keyboard to close the file.

```
less helloworld.txt
```

Press "q" to close the file.

There are many versions of these tools on command line, but "cat", "head" and "less" are very common.


### Copying Files
------------------

Sometimes we have a file that we want to reuse.

When copying within the same directory, make sure to change the name of the file, or the original will be **overwritten**.

When copying to a new directory, the name can stay the same.

This command copies the file within the same directory with a new name. Both files are kept.

```
cp helloworld.txt helloworld1.txt
```
Check this with `ls`

These commands make a new directory, and then copies the file into the new directory with the same name.

```
mkdir helloworld
cp helloworld.txt helloworld
```

Check this with `ls helloworld` (lists the contents of the directory).

### Moving Files
---------------------------

`mv` is an option for renaming files, but also has the potential to **overwrite** existing files.

For example, this command changes the name of the file and removes the original file. If `helloworld2.txt` already existed, it would be replaced.

```
mv helloworld1.txt helloworld2.txt
```

Check this with `ls`

### Removing Files
---------------------------

`rm` and `rmdir` are permanent in shell, so make sure you are ready to delete files.

```
rm helloworld/helloworld.txt
```

Once the directory is empty, we can remove the directory.

```
rmdir helloworld
```

It will throw an error if the directory is not empty.

If you are positive that you want to remove a directory and all the files within it, then add two flags, `-r` for recursive and `-f` for force.

Both commands above could have been replaced with one remove command: `rm -rf helloworld`

Until you are confident with file structure and bash commands, it is a good idea to copy instead of move and to set the interactive flag `-i` on the commands `rm` and `mv` to set up a question that you answer `y` or `n` to before removing.

```
rm -i helloworld/helloworld.txt
```

Generates this question
```
rm: remove regular file ???helloworld/helloworld.txt????
```

`mv -i` only generates a question if you are in danger of overwriting an existing file.

For example:

```
cp helloworld.txt helloworld1.txt
mv -i helloworld.txt helloworld1.txt
```
Generates the question:

```
mv: overwrite ???helloworld1.txt????
```


## Going Home
---------------------------


Sometimes we get lost, so it is useful to know a few ways to get back to where you started.

```
cd

```

This command returns you to your home directory. Check by typing this command.

```
pwd
```

Other options for going back to your home directory:

```
cd ~
```

```
cd $HOME
```

When lost in the file structure, going home is a good place to start.

### Find and Tree

File structures can get complicated quickly.

Two tools to understand where your files are that can help are `find` and `tree`.

From your home directory, you can find your file named `helloworld.txt` by typing the following:

```
cd
find . -name helloworld.txt
```

From the home directory, the answer should look like:

```
./JulyWorkshop/helloworld.txt
```

`.` in this case is another RELATIVE path direction that indicates "from this directoy that I am in currently". Note that the answer is given in the RELATIVE path format, starting with `.` = here.

It is also possible to provide an ABSOLUTE path to this command.

```
find /cluster/home/username01 -name helloworld.txt
```
This command will work from anywhere in the cluster. Note that the answer is given in the ABSOLUTE path format.

```
/cluster/home/arhode05/JulyWorkshop/helloworld.txt
```

Another helpful bash command for finding files is `tree`.

```
tree
```

This outputs your directory structure with lines that indicate the tree-like branches of your file structure.

This could be very messy if you already have a lot of files in a directory, so limit the level by adding a flag.

```
tree -L 2
```

This just shows the top two levels of the file structure.

[Return to TOC](#intro-to-command-line)

## Running Programs Interactively
----------------------------------

### HPC Etiquette
------------------

Try not to use the login computers for programs or large file management jobs. Looking things up and small commands such as `cat` or `head` are fine, but running programs may block others from logging in to the cluster.

### Switch to an Interactive Session
-------------------------------------

Do this first before running programs or testing your code.

```
srun -p batch --reservation=bioworkshop -n 2 --mem=8g -t 1-0 --pty bash
```

This command only works during the workshop on July 20 and 21. To use this command after the workshop is over or if you are working on your own, just remove the reservation flag.

* `-p` which partition to use, outside of class it is okay to use `interactive` or `batch`, your group may have it's own partition. You can read more about what is available by going to the OnDemand dropdown menu for "Misc" and look at "Scheduler Info" to find all the partition names.

* `--reservation` only applies during a specific class or workshop
* `-n` is the number of cpus to request, 2 or 4 is sufficient for most tests.
* `--mem` specifies the memory requested, 8g is usually sufficient for small jobs, consult the documentation for a program to find out if a minimum memory requirement is needed.
* `-t` indicates the time `1-0` means one day, so for tomorrow's session you will need to rerun this command.
* `--pty bash` just indicates that the shell opens in `bash`, meaning that all the commands that we learned today will work.


Let's go back into the JulyWorkshop directory, but this time use your ABSOLUTE path by changing `username01` to your username. If you forget your username, try `whoami`.

```
cd /cluster/home/username01/JulyWorkshop
```

---

#### Helpful Tip
================

There are some keyboard shortcuts that can help when writing complex commands and running programs interactively.

* Control-C will terminate a running process
* Control-A will put your cursor at the beginning of the line
* Control-E will put your cursor at the end of the line
* Up and down arrows will scroll through recent commands - If you make a mistake, just hit up to reveal the command and work on the part that was a mistake instead of retyping the whole thing.

**Extra Tip**: When trouble shooting a command using tickets, screen shots of error messages are a good option. (On Macs, Command-Shift-4)

---

### What is BLAST?

BLAST is the **B**asic **L**ocal **A**lignment **S**earch **T**ool.
It uses an index to rapdily search large sequence databases;
it starts by finding small matches between the two sequences and extending those matches.


![blst_0602](https://user-images.githubusercontent.com/8632603/180028827-f0980933-95f4-4c3b-b3d8-fece6e93eae8.gif)


For more information on how BLAST works and the different BLAST functionality,
check out the summary on [Wikipedia](https://en.wikipedia.org/wiki/BLAST) or
the NCBI's list of [BLAST resources](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs).

BLAST can be helpful for identifying the source of a sequence,
or finding a similar sequence in another organism.
In this lesson, we will use BLAST to find zebrafish proteins that
are similar to a small set of mouse proteins.

### Why use the command line?
BLAST has a very nice graphical interface for searching sequences in NCBI's database.
However, running BLAST through the commmand line has many benefits:
  * It's much easier to run many BLAST queries using the command line than the GUI
  * Running BLAST with the command line is reproducible and can be documented in a script
  * The results can be saved in a machine-readable format that can be analyzed later on
  * You can create your own databases to search rather than using NCBI's pre-built databases
  * It allows the queries to be automated
  * It allows you to use a remote computer to run the BLAST queries
  
Later on in the workshop we will talk more about these advantages and have a more in-depth explanation of the shell.

----


### Loading Modules

Many common programs are pre-loaded into the Tufts HPC using a system called "modules".

To see what versions of blast are available as a module, try running this command.

```
module av blast
```

As of July 2022, these are the modules you might see displayed.

<img width="711" alt="Blast_modules" src="https://user-images.githubusercontent.com/8632603/179539551-1d0c8933-30f2-43d5-957c-f4216d849ca6.png">


Choose the latest blast-plus version of the module and load it. 

```
module load blast-plus/2.11.0
```

When there is only one version of a module, the full version does not need to be provided, but it is always best to inclue the version as we are loading and updating versions of programs all of the time.

Confirm that the module is loaded.

```
module list
```

### Bringing in Files from the Internet

We need some data!  Let's grab the mouse and zebrafish RefSeq protein data sets from NCBI, and put them in our home directory. (this example is adapted from a lesson from Titus Brown's summer institute: https://angus.readthedocs.io/en/2019/running-command-line-blast.html)

For genomics projects, the files are often stored in pubic repositories and we must go and get those files before proceeding. These files originally came from the
[NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot), a copy has been placed in our github directory for future reference.


Now, we'll use `curl` to download the files from a Web site onto our computer. You will need to be connected to the internet for these commands to work.

* `-o` indicates this is the name we are assigning to our files in our own directory
* `-L` provides the full path for the download


It is possible to copy and paste both conmands to your terminal, they will run in sequence if there is not an error.

```
curl -o mouse.1.protein.faa.gz -L https://tuftsdatalab.github.io/Research_Technology_Bioinformatics/workshops/hpcForLifeSciences_July2022/IntroToLinux/mouse.1.protein.faa.gz

curl -o zebrafish.1.protein.faa.gz -L https://tuftsdatalab.github.io/Research_Technology_Bioinformatics/workshops/hpcForLifeSciences_July2022/IntroToLinux/zebrafish.1.protein.faa.gz
```

Another method for pulling files from the internet is `wget`, which will be demoed tomorrow. `curl` can pull more file types than `wget`, but in this simple case, either can be used.


If you look at the files in the current directory:

```
ls -l
```

You should now see these 3 files with details on who has permissions and when the files were created (notice that the dates are not today).

```
total 29908
-rw-rw-r-- 1 username01 username01 12553742 Jun 29 08:41 mouse.1.protein.faa.gz
-rw-rw-r-- 1 username01 username01 13963093 Jun 29 08:42 zebrafish.1.protein.faa.gz
```

The three files you just downloaded are the last three on the list - the
`.faa.gz` files.

All three of the files are FASTA protein files (that's what the .faa
suggests) that are compressed with `gzip` (that's what the .gz means).

Uncompress them:

```
gunzip *.faa.gz
```

and let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

These are protein sequences in FASTA format.  FASTA format is something
many of you have probably seen in one form or another -- it's pretty
ubiquitous.  It's a text file, containing records; each record
starts with a line beginning with a '>', and then contains one or more
lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll
do this using output redirection with '>', which says "take
all the output and put it into this file here."

`-n` flag for `head` specifies a number of lines to pull.

The first 11 lines contain two protein sequences. Let's extract those for blasting to test that our process is working.

```
head -n 11 mouse.1.protein.faa > mm-first.faa
```

```
cat mm-first.faa
```

Should produce:

```
>YP_220550.1 NADH dehydrogenase subunit 1 (mitochondrion) [Mus musculus domesticus]
MFFINILTLLVPILIAMAFLTLVERKILGYMQLRKGPNIVGPYGILQPFADAMKLFMKEPMRPLTTSMSLFIIAPTLSLT
LALSLWVPLPMPHPLINLNLGILFILATSSLSVYSILWSGWASNSKYSLFGALRAVAQTISYEVTMAIILLSVLLMNGSY
SLQTLITTQEHMWLLLPAWPMAMMWFISTLAETNRAPFDLTEGESELVSGFNVEYAAGPFALFFMAEYTNIILMNALTTI
IFLGPLYYINLPELYSTNFMMEALLLSSTFLWIRASYPRFRYDQLMHLLWKNFLPLTLALCMWHISLPIFTAGVPPYM
>YP_220551.1 NADH dehydrogenase subunit 2 (mitochondrion) [Mus musculus domesticus]
MNPITLAIIYFTIFLGPVITMSSTNLMLMWVGLEFSLLAIIPMLINKKNPRSTEAATKYFVTQATASMIILLAIVLNYKQ
LGTWMFQQQTNGLILNMTLMALSMKLGLAPFHFWLPEVTQGIPLHMGLILLTWQKIAPLSILIQIYPLLNSTIILMLAIT
SIFMGAWGGLNQTQMRKIMAYSSIAHMGWMLAILPYNPSLTLLNLMIYIILTAPMFMALMLNNSMTINSISLLWNKTPAM
LTMISLMLLSLGGLPPLTGFLPKWIIITELMKNNCLIMATLMAMMALLNLFFYTRLIYSTSLTMFPTNNNSKMMTHQTKT
KPNLMFSTLAIMSTMTLPLAPQLIT
```


Now let's BLAST these two sequences against the entire zebrafish
protein data set. First, we need to tell BLAST that the zebrafish
sequences are (a) a database, and (b) a protein database.  That's done
by calling 'makeblastdb':

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

Next, we call BLAST to do the search:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa
```

This should run pretty quickly, but you're going to get a lot of output!!
To save it to a file instead of watching it go past on the screen,
ask BLAST to save the output to a file that we'll name `mm-first.x.zebrafish.txt`:

```
blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
```

and then you can 'page' through this file at your leisure by typing:

```
less mm-first.x.zebrafish.txt
```

(Type spacebar to move down, and 'q' to get out of paging mode.)

-----



[Return to TOC](#intro-to-command-line)


## Writing a BASH Script and Running it as "Batch" 
--------------------------------------------------

In this example, we'll repeat the blast command above but refine it by outputting a table which summarizes each blast hit on one line. 

First, let's add more sequences to our query file. This will extract the first 186 sequences.  


```
head -n 999 mouse.1.protein.faa > mm-second.faa
```

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) for a description of the possible BLAST output formats.

In order to do this, we need to open a text editor.

### Opening a Text Editor

The easiest text editor to use on command line is `nano`, but there are many other types of command line text editors (`vi`,`emacs`,`vim`, etc.)

Nano is nice because it puts the instructions at the bottom of the editor in case you forget.

Open nano

```
nano

```

Hit Control-X to exit, say no and no. Nothing is saved, because we did not type into the file.

Let's reopen and copy and paste our script into the file.

Sometimes it is good to give a file name, so let's nano with a filename for our script.

```
nano sbatch.sh
```

Before closing, let's put some text into the file.

```
#!/bin/bash

#SBATCH --job-name=job
#SBATCH --nodes=1
#SBATCH -n 2
#SBATCH --partition=batch
#SBATCH --reservation=bioworkshop
#SBATCH --mem=8Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load blast-plus/2.11.0
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6


```

Control -X to close and save and use the same file name (sbatch.sh)

Because it is going to one or several virtual locations in the cluster, we need to reload the module as part of the script before running the script. This will make the command recognizable to the machine where the job is running.

```
cat sbatch.sh
```

Does it have all the elements?

If it does, a simple way to run it is by telling shell that it is a program.

```
sbatcb sbatch.sh
```



Delilah will explain the contents of this file, but let's go ahead and run it from the workshop directory.

Because we did not add any ABSOLUTE paths, then the sbatch command will look for the files where the program is running.

The results will also show up in that directory.

You can look at the output file with `less -S`, the flag allows scrolling from left to right instead of wrapping text or cutting it off:

```
less -S mm-second.x.zebrafish.tsv
```

(and again, type 'q' to get out of paging mode.)


The command line may move stuff around slightly, but it is a tab delimited file that can be downloaded to your computer and loaded into your spreadsheet program of choice.

`blastp` is a versatile tool for finding similar sequences, to see all the options, type `blastp -help`

#### Helpful Tip
=================

If writing the script on your laptop before copying and pasting, make sure to use a compatible text editor.

Even though you can't see it, popular word processors will add hidden symbols and change punctuation to your code.

There are several free tools available to avoid these errors.

* [Notepad+](https://notepad-plus-plus.org/downloads/) is free to download and use.
* [BBEdit](http://www.barebones.com/products/bbedit/) has a free version.

Other options are Subline and PyCharm, which have some features to help edit files.

------------------

## This concludes our Intro to Unix Lesson and we now Return to our Regularly Scheduled Slurm

[Return to TOC](#intro-to-command-line)

## Resources for Further Training in Command Line

* Udemy
* Coursera
* LinkedIn Learning

What are your favorites?








