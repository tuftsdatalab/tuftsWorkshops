# Advanced Linux/Unix Tools
## Awk
Awk is a powerful text-processing tool in Unix/Linux that allows you to manipulate and analyze text files and streams. It’s named after its creators (Aho, Weinberger, and Kernighan) and is commonly used for pattern scanning, processing, and reporting.

#### Syntax
```
awk 'pattern { action }' file
```

- **pattern**: A condition or regular expression that decides which lines are selected for processing. Patterns are  similar to  `if` statements in other languages: if the pattern’s expression evaluates to true or the regular expression matches, the statements inside **action** will run. **If we omit the pattern, Awk will run the action on all records**.

- **action**: Commands to be executed on the selected lines. **If we omit the action but specify a pattern, Awk will print all records that match the pattern**. 

- **file**: The text file to be processed.

  Awk processes input data a record at a time. Each record is composed of fields, separate chunks that awk automatically separates. Because awk was designed to work with `tabular data` each record is a line, and each field is a column’s entry for that record. The clever part about awk is that it automatically assigns the entire record to the variable `$0`, and field one’s value is assigned to `$1`, field two’s value is assigned to `$2`, field three’s value is assigned to `$3`, and so forth.	

#### Example
```
$ head -n 6 
gene_id	baseMean	log2FoldChange	lfcSE	pvalue	padj
ENSG00000000003	782.8404	-0.06662793	0.05691688	0.2130615	0.3800412
ENSG00000000419	746.4319	-0.1265914	0.06250695	0.02873058	0.07811519
ENSG00000000457	100.9565	0.1006396	0.14143475	0.3202419	0.5071726
ENSG00000000460	307.3198	0.1057537	0.08803984	0.1649367	0.3151837
ENSG00000000971	263.8175	-1.865821	0.11722751	6.524737e-59	4.470018e-57

$ awk 'BEGIN { FS="\t"; OFS="\t" } $6 < 0.05 && $3 > 1 { print $1, $3, $6 }' deseq2.results.tsv 
```

- **BEGIN { FS="\t"; OFS="\t" }**: Sets the input (`FS`) and output (`OFS`) field separators to tab (`\t`) since your file is a TSV (tab-separated values).
- **$6 < 0.05 && $3 > 1**: Filters rows where the padj (6th column) is less than 0.05 and the log2FoldChange (3rd column) is greater than 1.
- **{ print $1, $3, $6 }**: Prints the gene_id (1st column), log2FoldChange (3rd column), and padj (6th column).

## GNU Parallel
GNU Parallel is a command-line tool designed to execute shell commands or scripts in parallel on a local or remote system. It is especially useful for bioinformatics, data processing, and other fields that involve repetitive command execution, as it can significantly speed up tasks by utilizing multiple CPU cores.

### Basic syntax
#### Triple colon:::
```
$ parallel [options] command ::: arguments
```

- **::: arguments**: The list of arguments passed to the command. Each argument is passed to the command in parallel.


#### Quad colon::::
```
$ parallel [options] command :::: input_file
```

`-a` is alternative syntax to quadruple colon.
```
$ parallel [options] command -a input_file
```

Run `command` in parallel for each line in input_file

#### Pipe
```
$ command1 | parallel [options] command2
```
standard output from command1 as argument


#### Example
```
parallel -j N "fastqc {}" ::: *.fastq.gz
```

- **-j N**: Specifies the number of parallel jobs to run (replace N with the desired number, considering available CPU cores).
- **"fastqc {}"**: The FastQC command to execute in parallel, with {} representing each input file.
- **:::**: Separates the command from the list of files.
- ***.fastq.gz**: Wildcard pattern to match all FASTQ files with the **.fastq.gz** extension in the current directory. Modify as needed for different file extensions or locations.

Highly recommeded to read this article written by the developer Ole Tange in [Biostars](https://www.biostars.org/p/63816/). 
 <img src="http://i.stack.imgur.com/17FsG.png" alt="GNU parallel" style="height:500px;" />
