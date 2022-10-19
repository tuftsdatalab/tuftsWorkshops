### Return to the workshop directory

```
cd ~/Oct22Workshop
```

### Loading Modules

Many common programs are pre-loaded into the Tufts HPC using a system called "modules".

To see what versions of blast are available as a module, try running this command. 

!!! tip
  
  You can use the first part of the program name to check if there is a module.

```
module av blast
```

As of October 2022, these are the modules you might see displayed.

<img width="728" alt="moduleavblast" src="https://user-images.githubusercontent.com/8632603/196236539-c2308c9f-bbf8-44e2-9679-731f7299cf9b.png">


Choose the latest blast-plus version of the module and load it. 

```
module load blast-plus/2.11.0
```

When there is only one version of a module, the full version does not need to be provided, but it is always best to inclue the version as we are loading and updating versions of programs all of the time.

Confirm that the module is loaded.

```
module list
```

tmux and blast should be listed.

If other programs are loaded with the module, they may also show up with this command.



### Bringing in Files from the Internet

We need some data!  Let's grab the mouse and zebrafish RefSeq protein data sets from NCBI, and put them in our home directory. (this example is adapted from a lesson from [Titus Brown's summer institute](https://angus.readthedocs.io/en/2019/running-command-line-blast.html){:target="_blank" rel="noopener"}. These lessons contain a lot of command line examples.

!!! note
  
  `curl` and `wget` are the two most common tools used to bring in files that are available from a url. We are going to use `curl` because that command works well for files coming from an `ftp://` url.

??? note "Copying files over from NCBI"
  For genomics projects, the files are often stored in pubic repositories and we must go and get those files before proceeding. These files originally came from the
  [NCBI FTP site](ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot){:target="_blank" rel="noopener"}, a copy has been placed in our github directory for future reference.
  
 
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
Compressed files may have a different color when you use the `ls` command.

Uncompress the files.

```
gunzip *.faa.gz
```

Because both files follow a very similar pattern, and we want to decompress all our .gz files, we can use the `*` wildcard (filenames that have a pattern that matches and number of missing letters before the part of the file name that is the same

!!! note "Regular Expressions"
  
  `*` and other wildcards are useful to save on typing scripts, because many actions can be combined in one request.
  
  Regular Expressions are a set of special characters combined with unix commands.
  
  Here is a [link](https://www.grymoire.com/Unix/Regular.html) that explains the basic syntax){:target="_blank" rel="noopener"}. 
  

## Checking the contents of a File

We've already used `cat` and `less` to look at the content of our helloworld.txt files. Some files are very large and we may only want to check the first few lines to reassure ourselves that the download worked correctly.

Let's look at the first few sequences in the file:

```
head mouse.1.protein.faa 
```

!!! note "FASTA format

  These are protein sequences in FASTA format.  FASTA format is something
  many of you have probably seen in one form or another -- it's pretty
  ubiquitous.  It's a text file, containing records; each record
  starts with a line beginning with a '>', and then contains one or more
  lines of sequence text.

Let's take those first two sequences and save them to a file.  We'll do this using output redirection with '>', which says "take
all the output and put it into this file here."

```
head -n 11 mouse.1.protein.faa > mm-first.faa
```

`-n` flag for `head` specifies a number of lines to pull.

The first 11 lines contain two protein sequences. Let's extract those for blasting to test that our process is working.


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
by calling `makeblastdb`

```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```

`makeblastdb` is a program that was loaded using the `module` command. If you unload the module, this command may not work.


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

!!! question "What are your questions?"

!!! note
  
  This command was an example of `interactive` shell scripting because we are typing in the commands manually and waiting for the results. If we walk away from our machine and the session times out, then the program may be interrupted. `tmux` allows us to keep running the program even if we take a break.
  
  The next session demonstrates how to combine all of these commands into a script that runs on **SLURM**. 
  
  **SLURM** differs from `interactive` computing because you activate the script instead of manually writing the commands one at a time.


