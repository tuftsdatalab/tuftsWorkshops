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
