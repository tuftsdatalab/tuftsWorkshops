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