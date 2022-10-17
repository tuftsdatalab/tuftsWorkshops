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

This command only works during the October-November 2022 workshops. To use this command after the workshop is over or if you are working on your own, just remove the reservation flag.

* `-p` which partition to use, outside of class it is okay to use `interactive` or `batch`, your group may have it's own partition. You can read more about what is available by going to the OnDemand dropdown menu for "Misc" and look at "Scheduler Info" to find all the partition names.

* `--reservation` only applies during a specific class or workshop
* `-n` is the number of cpus to request, 2 or 4 is sufficient for most tests.
* `--mem` specifies the memory requested, 8g is usually sufficient for small jobs, consult the documentation for a program to find out if a minimum memory requirement is needed.
* `-t` indicates the time `1-0` means one day, so for tomorrow's session you will need to rerun this command.
* `--pty bash` just indicates that the shell opens in `bash`, meaning that all the commands that we learned today will work.



### Finding your files interactively

When you request a computer using an `srun` command, the beginning of your command line should change to indicate that you are no longer on a `login` node and instead are on a `compute` node. It will tell you which node you are on.

!!! question

 What compute node are you on? Type it into the chat box.
 

Your files will `mount` to the new node. This means that you can be on any computer in the Tufts HPC and it will recognize your home directory structure.


Let's go back into the Oct22Workshop directory, but this time use your **ABSOLUTE** path by changing `username01` to your username. If you forget your username, try `whoami`.

```
cd /cluster/home/username01/Oct22Workshop
```

### Find and Tree

File structures can get complicated quickly.

Two tools to understand where your files are that can help are `find` and `tree`.

From your home directory, you can find your file named `helloworld.txt` by typing the following:

```
cd
find . -name helloworld.txt
```

`find` is a bash command
`.` means look from this location and into any subdirectories to this location
`name` is the file that you are looking for




From the home directory, the answer should look like:

```
./Oct22Workshop/helloworld.txt
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

!!! tip

 * `find` using the parameter `iname` allows the search to be insensitive to case `find . -iname "helloworld.txt"` finds `Helloworld.txt` AND `helloworld.txt`
 * wildcards allow for partial searches of many filenames that may match.
 * The easiest wildcard is `*` which means any number of characters can match, such as `find . -iname "hello*.txt"` finds any file that begins with `hello` and has any number of characters before `.txt` finds `helloworld1.txt` AND `helloworld2.txt`. 
 * `*` can be used anywhere in the pattern: `hello*.txt`, `hello*.txt`, `*.txt`


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


---

!!! tip

 There are some keyboard shortcuts that can help when writing complex commands and running programs interactively.

 * Control-C will terminate a running process
 * Control-A will put your cursor at the beginning of the line
 * Control-E will put your cursor at the end of the line
 * Up and down arrows will scroll through recent commands - If you make a mistake, just hit up to reveal the command and work on the part that was a mistake instead of retyping the whole thing.

!!! note

 When trouble shooting a command using tickets, screen shots of error messages are a good option. (On Macs, Command-Shift-4)

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
  
We are next going to write a script that we will send to SLURM which will demonstrate these advantages.

----
