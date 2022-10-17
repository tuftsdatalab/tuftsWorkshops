## Setting Parameters for Bash Commands
--------------------------------

Many bash commands have special **parameters**, sometimes referred to as **flags** that open up a lot more possibilities.

Let's start by going to your home directory (you choose the command)



??? note "What if I don't have any files to list? 

    You can add an empty file and an empty directory by typing these two commands:

    The command `touch` makes an empty text file. The command `mkdir` makes a new directory.

    ```
    touch emptyfile.txt
    ```

    ```
    mkdir emptydir
    ```



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
