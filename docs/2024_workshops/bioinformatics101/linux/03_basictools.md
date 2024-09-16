# Must-known Linux/Unix tools

##  File and Directory Management
Linux provides powerful tools for managing files and file systems. Here we will introduce a few essential commands. 

### pwd: print the current working directory
#### Usage
```
$ pwd
/cluster/home/yzhang85
$ cd /cluster/tufts/rt/yzhang85/
$ pwd
/cluster/tufts/rt/yzhang85
```

### cd: change directory
#### Usage 
```
cd [directory]
```
If a directory is not supplied as an argument, it will default to your **home** directory. 
```
$ pwd
/cluster/tufts/rt/yzhang85
$ cd ..
$ pwd
/cluster/tufts/rt
$ cd  
$ pwd
/cluster/home/yzhang85
```

### touch: create new files and update timestamps
touch is used to create new files or to update the timestamps (access and modification times) of existing files. 
#### Create new file
```
touch newfile.txt
```
#### Update Timestamps of Existing Files
```
touch existingfile.txt
```

### mkdir: create new directory
#### Usage 
```
mkdir [options] dir_name
```

#### Common option
​- **-p**: Creates parent directories if they don't exist.

```
$ mkdir -p rnaseq/output 
```
This will create output folder as well as its parent folder rnaseq if it doesn't exist.

### mv: move a file/directory to a new location or rename it

#### Usage
```
mv [options] source destination
```
#### Common option
- **-i**: Prompts for confirmation before overwriting an existing file. Useful to avoid accidental data loss.- 
- **-f**: Forces the operation without prompting, even if an existing file would be overwritten. Use with caution!

### cp: copy a file/directory
#### Usage
```
cp [options] source destination
```
#### Common option
- **-r**:  To copy directory

### rm: remove files/directories
#### Usage 
```
rm [options] file/directory
```
#### Common option
- **-r:** Deletes recursively any file and subdirectories contained within the given directory

## Text processing
Linux command-line tools are invaluable for bioinformatics text processing due to their efficiency and flexibility. They allow for rapid manipulation and analysis of large biological datasets, such as DNA sequences, protein structures, and gene expression data. Commands like grep, sed, awk, and cut are essential for filtering, extracting, and reformatting text-based biological information.

### cat: Catenate files (joins their contents)**
#### Usage
```
cat [options] file1 file2 …
```
#### Common option
- **-n:** tells cat to **number each line of the output**. This is helpful for debugging scripts.

### head/tail: Displays the beginning/end of a file
#### Usage
```
head/tail [options] file
```
#### Common option
- **-n** [number]: Specifies the number of lines to display (default: 10).

### less/more: View the content of a file page by page
#### Usage
```
less largefile.txt
more largefile.txt
```

### cut: Extract sections from each line of files
cut is a text-processing utility used to extract specific sections of lines from a file or standard input. It is commonly used to split lines of text based on delimiters, extract columns of data, and work with fixed-width fields. The cut command is helpful for manipulating text files like TSVs and CSVs.

#### Usage
```
cut [OPTIONS] [FILE...]
```
#### Common options
- **-d**: Specifies the delimiter for field extraction. Default is TAB. Example: `-d,` (use a comma as the delimiter).
- **-f**: Selects specific fields, used with a delimiter (default is TAB). Example: `-f 1,3` (extract fields 1 and 3).

#### Example
```
cut -f1,3 -d, file.csv ##(Extract columns 1 and 3 from a comma-separated file)
```

### sort: Sort lines of text files
 `sort` is designed to sort plain-text data with columns. Running `sort` without any arguments simply sorts a file alphabetically.
 By default, sort treats blank characters (like tab or spaces) as field delimiters. If your file uses another delimiter (such as a comma for CSV files), you can specify the field separator with `-t` (e.g., `-t","`). 

#### Usage
##### Sort lines alphabetically
```
sort file.txt
```
##### Sort lines numerically
```
sort -n file.txt
```

##### Sort by a specific column
```
sort -k 2 file.txt  # Sort by the second column
```

##### Sort by multiple columns
```
sort -k 1,2 file.txt  # Sort by the first column, then the second
```

### uniq: Report or filter out repeated lines in a file.
#### Usage
```
sort file.txt | uniq
```

### grep:Extracting lines matching (not matching) a pattern**
#### Usage
 ```
 grep [options] PATTERN file
```
#### Common option
- **-i**: ignore cases
- **-v**: select non-matching lines.
- **-A NUM:** Print **NUM** lines of trailing context after matching lines.
- **-B NUM:** Print **NUM** lines of leading context before matching lines.

### sed: Stream editor for modifying file content
sed (short for stream editor) is a powerful text-processing tool in Bash that allows you to parse and transform text in files or streams. It is commonly used to perform basic text manipulations like search and replace, insert and delete lines, and apply regular expressions on text data.
#### Substitution (Search and Replace)
Replace the first occurrence of **old** with **new** in each line:
```
sed 's/old/new/' filename.txt
```
Replace **all** occurrences of **old** with **new** in each line:
```
sed 's/old/new/g' filename.txt
```

#### In-Place Substitution
```
sed -i 's/old/new/g' filename.txt
```
Warning: Use this command with caution as it directly modifies the original file. To create a backup, use `-i.bak`:
```
sed -i.bak 's/old/new/g' filename.txt
```
#### Delete Lines
```
sed '/pattern/d' filename.txt
```

`sed` and `awk` that we will introduce later are very powerful bash commands. If you have interest in learning more about their usage, below is a very good book.
![se and awk](https://m.media-amazon.com/images/I/916OGXYyYVL._SL1500_.jpg)

## Data Compression and Archiving
When working with files on Linux, compressing them to save space and bundling multiple files into a single archive is a common practice. The commands gzip, gunzip, and tar are essential tools for file compression and archiving in Bash.
### gzip: Compress files
gzip is used to compress files in Linux. It reduces file sizes using the **DEFLATE** algorithm, resulting in files with the `.gz` extension.
#### Usage
##### Compress file/files
```
gzip file.txt ## Compress file.txt and replace it with file.txt.gz.
gzip file1.txt file2.txt file3.txt ## Each file will be compressed and replaced with a .gz version.
```

##### Keep the Original File
```
gzip -k file.txt ## Create file.txt.gz while preserving the original file.txt.
```
##### Compress the folder
```
gzip -r directory_name
```
### gunzip: decompress .gz files
##### Decompress files
```
gunzip file.txt.gz
gunzip file1.txt.gz file2.txt.gz
```
##### Keep the Compressed File
```
gunzip -k file.txt.gz
```


### tar: Archive multiple files into one or extract them.
#### Usage
```
tar -cvf archive.tar directory/ ## Create an archive
tar -xvf archive.tar ## Extract an archive
```

 
## Other useful tools

### Redirection: >, >>, <

- `>`: Overwrites the contents of a file with the command's output

 	`cat file1 file2 > files`

- `>>`: Appends the output to the end of an existing file

​        `cat file3 >> files`

- `<`: Uses the contents of a file as input to a command

​	`sort < names.txt`

### Pipe: |

Pipes in Linux are a powerful feature that allows you to connect the output of one command directly as the input to another command. This is a key concept in Unix/Linux philosophy, which promotes the use of small, modular tools that can be combined to perform complex tasks.

A pipe is represented by the `|` symbol. When you place a pipe between two commands, the standard output (`stdout`) of the command on the left of the pipe becomes the standard input (`stdin`) for the command on the right.

#### Usage

```
command1 | command2
```


```
sort file.txt | uniq
```

​	•	sort file.txt: Sorts the lines in file.txt.

​	•	uniq: Removes duplicate lines from the sorted output.


### Wildcards: selecting multiple files/directories based on patterns

- **\***: Represents zero or more characters.

 		*.fastq.gz  matches all fastq.gz files

- **?**: Represents a single character.

 		file?.txt matches "file1.txt", "fileA.txt", but not "file12.txt".

- **[]**: Represents a single character within a specified range or set.

​		 [abc]at matches "bat", "cat", or "aat”.

- **[0-9]** matches any single digit.

### Alias

An alias in Linux is a custom shortcut or abbreviation for a command or a series of commands. Once defined, you can use the alias in place of the original command.

#### Creating an alias

To create an alias, use the alias command followed by the name you want to give the alias and the command it should execute.

```
alias alias_name='command'
```

#### Example

```
alias ll='ls -l'
alias la='ls -a'
alias mav='module avail'
alias ml='module load'
```

#### Using alias

```
[yzhang85@login-prod-03 MPI]$ ll
total 31
-rwxrwx--- 1 yzhang85 yzhang85 16256 May 16 07:57 hello_mpi*
-rw-rw---- 1 yzhang85 yzhang85   731 May 16 07:56 hello_mpi.c
-rwxrwxr-x 1 yzhang85 yzhang85  8392 Dec 23  2023 hello_world*
-rwxrwxr-x 1 yzhang85 yzhang85   731 Oct 31  2023 hello_world.c*
-rwxrwxr-x 1 yzhang85 yzhang85  4053 Oct 31  2023 mpi_hello_world.c*

[yzhang85@login-prod-03 MPI]$ mav blast
-------------------------------- /opt/shared/Modules/modulefiles-rhel6 --------------------------
   blast/2.2.24    blast/2.2.31    blast/2.3.0    blast/2.8.1

------------------------------- /cluster/tufts/hpc/tools/module ---------------------------------
   blast-plus/2.11.0    ncbi-magicblast/1.5.0

------------------------------- /cluster/tufts/biocontainers/modules -----------------------------
   blast/2.15.0 (D)

  Where:
   D:  Default Module
```



### ln -s: softlink
A `soft link` (also known as a symbolic link or symlink) is a type of file in Linux that points to another file or directory. It’s essentially a shortcut that references the location of another file, allowing you to access it from a different location in the filesystem.
#### Usage
To create a soft link, you use the `ln` command with the `-s` option:
```
ln -s target_file link_name
```
- target_file: The file or directory you want to link to.
- link_name: The name of the symlink that will point to the target.


[Previous: Files and File system](02_files.md)                                                                 

[Next: Loops](04_loops.md)