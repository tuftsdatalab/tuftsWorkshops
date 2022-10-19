## Writing a BASH Script and Running it as "Batch" 
--------------------------------------------------

In this example, we'll repeat the blast command above but refine it by outputting a table which summarizes each blast hit on one line. 

### Return to the Workshop Directory

cd ~/Oct22Workshop


First, let's add more sequences to our query file. This will extract the first 186 sequences.  


```
head -n 999 mouse.1.protein.faa > mm-second.faa
```

See [this link](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6){:target="_blank" rel="noopener"} for a description of the possible BLAST output formats.

In order to do this, we need to open a text editor.


## Opening a Text Editor

The easiest text editor to use on command line for beginners is `nano`, but there are many other types of command line text editors (`vi`,`emacs`,`vim`, etc.)

Nano is nice because it puts the instructions at the bottom of the editor in case you forget.

Open nano

```
nano
```

**Control-X** to exit, say **no** and **no**. Nothing is saved, because we did not type into the file.

Let's reopen and copy and paste our script into the file.

Sometimes it is good to give a file name, so let's nano with a filename for our script.

```
nano blast_sbatch.sh
```

Before closing, let's put some text into the file. 

Make sure to change the email address to your own email.

```
#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH -n 2
#SBATCH --partition=batch
#SBATCH --reservation=bioworkshop
#SBATCH --mem=8Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-user=youremail@tufts.edu

module load blast-plus/2.11.0
blastp -query mm-second.faa -db zebrafish.1.protein.faa -out mm-second.x.zebrafish.tsv -outfmt 6


```

**Control -X** to close and save and use the same file name (blast_sbatch.sh)

Because it is going to one or several virtual locations in the cluster, we need to reload the module as part of the script before running the script. This will make the command recognizable to the machine where the job is running.

```
cat blast_sbatch.sh
```

Does it have all the elements?

If it does, a simple way to run it is by telling shell that it is a program to run on SLURM.

```
sbatch blast_sbatch.sh
```



Let's go ahead and run it from the workshop directory where you copied your data to.

Because we did not add any ABSOLUTE paths, then the sbatch command will look for the files where the program is running.

The results will also show up in that directory.

You can look at the output file with `less -S`, the flag allows scrolling from left to right instead of wrapping text or cutting it off:

```
less -S mm-second.x.zebrafish.tsv
```

(and again, type 'q' to get out of paging mode)

The command line may move stuff around slightly, but it is a tab delimited file that can be downloaded to your computer and loaded into your spreadsheet program of choice.

`blastp` is a versatile tool for finding similar sequences, to see all the options, type `blastp -help`

!!! tip

  If writing the script on your laptop before copying and pasting, make sure to use a compatible text editor.

  Even though you can't see it, popular word processors will add hidden symbols and change punctuation to your code.

  There are several free tools available to avoid these errors.

  * [Notepad+](https://notepad-plus-plus.org/downloads/) is free to download and use.
  * [BBEdit](http://www.barebones.com/products/bbedit/) has a free version.

  Other options are Sublime and PyCharm, which have some features to help edit files.

------------------

## This concludes our Intro to Unix Lesson and we now Return to our Regularly Scheduled Slurm

[Return to TOC](#intro-to-command-line)

## Resources for Further Training in Command Line

* Udemy (free to the Tufts community)
* Coursera
* LinkedIn Learning

What are your favorites?
