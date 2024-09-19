Shirley Li

Date: 2024-10-

xue.li37@tufts.edu



[TOC]



# Objectives



- The primary goal of this tutorial is to introduce participants to bioinformatics on a high-performance computing (HPC) cluster.
-  By the end of this tutorial, participants will: Understand basic command-line operations.  
- Be familiar with common bioinformatics file formats (e.g., FASTA, FASTQ, SAM, BAM). - Gain hands-on experience with bioinformatics tools and workflows.



# Getting Started



## Prerequisites

- Basic understanding of biology and bioinformatics 
- Familiarity with the command line 
- Access to an HPC cluster (e.g., login credentials, necessary software installations)



## Setup

1. **Connecting to the Cluster** 
   1. Instructions for accessing the cluster (e.g., SSH login)
   2. Overview of the cluster environment

2. **Setting Up the Environment**

   1. Start an interactive session, go from log in mode to compute mode. 
   2. Creating and managing directories for the tutorial

3. **Copying Sample Data**

   1. Go to your directory (Do not use Home directory). 

      You can use your lab storage 

      ```
      cd /cluster/tufts/XXlab/utln/
      ```

      or the dedicated space for workshops (Warning: This directory is accessible for many users and files will be deleted if keep untouched for 2 months.)

      ```
      cd /cluster/tufts/workshop/yourutln/
      ```

      Contact us at tts-research@tufts.edu if your lab doesn't have storage or you don't belong to any research lab, we will add you to `/cluster/tufts/workshop/`.

   2. Copy example data to your directory. 

      ```
      cp /cluster/tufts/workshop/exampledata/ ./
      ```

      or 

      ```
      cp -r /cluster/tufts/workshop/demo/bio_2024/week1/ /cluster/tufts/workshop/yourutln/
      ```

# Overview of the input files



## FASTA file

- FASTA is a text-based format for representing **nucleotide** sequences or **protein** sequences. It is widely used in bioinformatics for sequence data storage and analysis. Each sequence in a FASTA file is represented by a header line starting with a '>', followed by lines of sequence data.

- A typical FASTA file has the following structure: 
  - **Header Line**: This line starts with a '>' character, followed by a description or identifier of the sequence. 
  - **Sequence Data**: The nucleotide or amino acid sequence, which can span multiple lines.

1. Example

   ```
   cat other/rcsb_pdb_5XUS.fasta
   ```

   What does this file contain?

```
>5XUS_1|Chain A|LbCpf1|Lachnospiraceae bacterium ND2006 (1410628)
GSHMSKLEKFTNCYSLSKTLRFKAIPVGKTQENIDNKRLLVEDEKRAEDYKGVKKLLDRYYLSFINDVLHSIKLKNLNNYISLFRKKTRTEKENKELENLEINLRKEIAKAFKGNEGYKSLFKKDIIETILPEFLDDKDEIALVNSFNGFTTAFTGFFDNRENMFSEEAKSTSIAFRCINENLTRYISNMDIFEKVDAIFDKHEVQEIKEKILNSDYDVEDFFEGEFFNFVLTQEGIDVYNAIIGGFVTESGEKIKGLNEYINLYNQKTKQKLPKFKPLYKQVLSDRESLSFYGEGYTSDEEVLEVFRNTLNKNSEIFSSIKKLEKLFKNFDEYSSAGIFVKNGPAISTISKDIFGEWNVIRDKWNAEYDDIHLKKKAVVTEKYEDDRRKSFKKIGSFSLEQLQEYADADLSVVEKLKEIIIQKVDEIYKVYGSSEKLFDADFVLEKSLKKNDAVVAIMKDLLDSVKSFENYIKAFFGEGKETNRDESFYGDFVLAYDILLKVDHIYDAIRNYVTQKPYSKDKFKLYFQNPQFMGGWDKDKETDYRATILRYGSKYYLAIMDKKYAKCLQKIDKDDVNGNYEKINYKLLPGPNKMLPKVFFSKKWMAYYNPSEDIQKIYKNGTFKKGDMFNLNDCHKLIDFFKDSISRYPKWSNAYDFNFSETEKYKDIAGFYREVEEQGYKVSFESASKKEVDKLVEEGKLYMFQIYNKDFSDKSHGTPNLHTMYFKLLFDENNHGQIRLSGGAELFMRRASLKKEELVVHPANSPIANKNPDNPKKTTTLSYDVYKDKRFSEDQYELHIPIAINKCPKNIFKINTEVRVLLKHDDNPYVIGIDRGERNLLYIVVVDGKGNIVEQYSLNEIINNFNGIRIKTDYHSLLDKKEKERFEARQNWTSIENIKELKAGYISQVVHKICELVEKYDAVIALEDLNSGFKNSRVKVEKQVYQKFEKMLIDKLNYMVDKKSNPCATGGALKGYQITNKFESFKSMSTQNGFIFYIPAWLTSKIDPSTGFVNLLKTKYTSIADSKKFISSFDRIMYVPEEDLFEFALDYKNFSRTDADYIKKWKLYSYGNRIRIFRNPKKNNVFDWEEVCLTSAYKELFNKYGINYQQGDIRALLCEQSDKAFYSSFMALMSLMLQMRNSITGRTDVDFLISPVKNSDGIFYDSRNYEAQENAILPKNADANGAYNIARKVLWAIGQFKKAEDEKLDKVKIAISNKEWLEYAQTSVKH

>5XUS_2|Chain B|crRNA|synthetic construct (32630)
AAUUUCUACUAAGUGUAGAUGGAAAUUAGGUGCGCUUGGC

>5XUS_3|Chain C|DNA (29-MER)|synthetic construct (32630)
GCCAAGCGCACCTAATTTCCTAAAGGACG

>5XUS_4|Chain D|DNA (5'-D(*CP*GP*TP*CP*CP*TP*TP*TP*A)-3')|synthetic construct (32630)
CGTCCTTTA
```

Counting the number of sequencies in a FASTA file:

```
grep -c "^>" sequences.fasta
```

## FASTQ



- FASTQ is a text-based format used to store both nucleotide sequences and their corresponding quality scores. It is widely used in bioinformatics, particularly for storing data from high-throughput sequencing technologies.
- Structure of a FASTQ File A FASTQ file consists of a series of entries, each representing a single read. Each entry has four lines: 
  - **Header Line**: Starts with '@' followed by a sequence identifier and an optional description.  
  - **Sequence Line**: The raw sequence of nucleotides (A, T, C, G).  
  - **Separator Line**: Starts with a '+' character and can optionally be followed by the same sequence identifier and description as in the header line. 
  - **Quality Line**: Encodes the quality scores for the sequence in the sequence line, using ASCII characters.

Example:

```
head raw_fastq/Irrel_kd_1.subset.fq
```

Warning: Do not use `cat` to view fq files, as they are often very large. 

```
@HWI-ST330:304:H045HADXX:1:1214:19169:32660
CTTTTTTGCTGGAACTTTAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTACCCTTAGTACCAGGTGCTTTTTTGGGAGAAGCT
+
CCCFFFFFGGHGHJJJJJJJIEHIIIIIGIIIJJHIJIIJJJJJGGGIJJJJJJJJJJGHGHHFBCEFFEEEEE;@@CDDD5>ACDDDDDDBB?BDDDCA
@HWI-ST330:304:H045HADXX:2:1212:14280:80867
CCCTAATGATGATATATGGATCAAAAGTCTTCTTTGTAGTACAAACAGTCATGCTGCCTTCGATCAGGTCCAGGGTTGCATTAACATGATGTTCATTTAA
+
:?@DDDDFDHHGHJIGBFGIGIIIGDE:CFIHIIJE1?>CEFGCEGFG?FGI?F@FHCCGIHIIJGHC@FHHGEB@9AEEHFFFEDFCCACEEEEDDDDD
@HWI-ST330:304:H045HADXX:2:2111:14933:89958
GGCATCCATGTTCTTGCCCAAAACCTTGGTTACAGCAATCTGATACTTCTTTTGTGTGGGCTGGCATAGGTCAATGAGGCAGATCGGAAGAGCACACGTC
```

You can also use `less` to view the file and press `Q` key to exit.  



## GTF file



GTF (Gene Transfer Format) is a file format used to hold information about gene structure. It is widely used in genomics to store annotations of genomic features such as genes, exons, and regulatory elements. GTF is similar to GFF (General Feature Format) but includes additional standardized attributes. 

Structure of a GTF File 

A GTF file is a tab-delimited text file with one line per feature. Each line consists of nine fields: 

1. **seqname**: Name of the sequence (e.g., chromosome). 
2. **source**: Name of the program that generated the feature. 
3. **feature**: Type of feature (e.g., gene, exon, CDS). 
4. **start**: Start position of the feature. 
5. **end**: End position of the feature. 
6. **score**: Score value (can be a dot if not used). 
7. **strand**: Strand of the feature (+ or -). 
8. **frame**: Frame for coding sequences (0, 1, or 2). 
9. **attribute**: A semicolon-separated list of key-value pairs describing the feature.

```
less reference_data/chr1-hg19_genes.gtf 
```





```
chr1    unknown exon    14362   14829   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    14970   15038   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    15796   15947   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    16607   16765   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    16858   17055   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    17233   17368   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    17606   17742   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    17915   18061   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    18268   18366   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    24738   24891   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    29321   29370   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    34611   35174   .       -       .       gene_id "FAM138F"; gene_name "FAM138F"; transcript_id "NR_026820"; tss_id "TSS8099";
chr1    unknown exon    34611   35174   .       -       .       gene_id "FAM138A"; gene_name "FAM138A"; transcript_id "NR_026818"; tss_id "TSS8099";
chr1    unknown exon    35277   35481   .       -       .       gene_id "FAM138F"; gene_name "FAM138F"; transcript_id "NR_026820"; tss_id "TSS8099";
chr1    unknown exon    35277   35481   .       -       .       gene_id "FAM138A"; gene_name "FAM138A"; transcript_id "NR_026818"; tss_id "TSS8099";
chr1    unknown exon    35721   36081   .       -       .       gene_id "FAM138F"; gene_name "FAM138F"; transcript_id "NR_026820"; tss_id "TSS8099";
chr1    unknown exon    35721   36081   .       -       .       gene_id "FAM138A"; gene_name "FAM138A"; transcript_id "NR_026818"; tss_id "TSS8099";
chr1    unknown CDS     69091   70005   .       +       0       gene_id "OR4F5"; gene_name "OR4F5"; p_id "P9488"; transcript_id "NM_001005484"; tss_id "TSS13903";
chr1    unknown exon    69091   70008   .       +       .       gene_id "OR4F5"; gene_name "OR4F5"; p_id "P9488"; transcript_id "NM_001005484"; tss_id "TSS13903";
chr1    unknown start_codon     69091   69093   .       +       .       gene_id "OR4F5"; gene_name "OR4F5"; p_id "P9488"; transcript_id "NM_001005484"; tss_id "TSS13903";
chr1    unknown stop_codon      70006   70008   .       +       .       gene_id "OR4F5"; gene_name "OR4F5"; p_id "P9488"; transcript_id "NM_001005484"; tss_id "TSS13903";
ch
```

GTF files are essential for common bioinformatics analyses such as RNA-Seq analysis, genome annotation, and differential expression analysis. They provide the necessary annotation information to map reads to genomic features and to understand the functional elements of the genome.

You can download human gtf file from here `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/`



# Toy analysis

## Alignment with bowtie2

### Make sure you are in compute node

```
```

If not, do the following

```
```



### Load necessary modules

```
module load bowtie2/
```





## Output 

### SAM/BAM file

SAM (Sequence Alignment/Map) is a text-based format for storing biological sequences aligned to a reference sequence. It is widely used in bioinformatics for storing large-scale sequencing data, such as that from next-generation sequencing technologies.

A SAM file consists of a header section and an alignment section. 

1. **Header Section**: Optional. Contains metadata about the alignments and reference sequences. Header lines start with the '@' character. 
2. **Alignment Section**: Contains the aligned sequences and their corresponding information. Each line represents a single alignment and consists of 11 mandatory fields and optional fields.



### Manipulate the file





**Materials adapted from https://hbctraining.github.io/Intro-to-shell-flipped/lessons/03_working_with_files.html**

