## Making files and directories
--------------------------------

Many bash commands have special **parameters**, sometimes referred to as **flags** that open up a lot more possibilities.

Let's start by going to your home directory (you choose the command)

For new users, this may not return any content (besides `privatemodules` if you loaded `tmux` at the beginning.)

Let's make a workshop directory and put a file into it.

1.) Go to your home directory 

```
cd
```

`cd` not only changes directory, it allows you to go home by typing the command all by itself without a directory name.

2.) Create a new directory for the workshop

```
mkdir Oct22Workshop
```
!!! note

    `mkdir` is a specific command that allows you to make a directory.
    `rmdir` is a command that allows you to remove a directory (but only if it is empty)
     When nameing files and directories, avoid spaces and special characters except underscores ("_") and hyphens ("-").


!!! note "Important"

    **Spelling** and **Capitalization** are literal in unix.
    Be careful when making and using files to be consistent in your process. 
    This will make it easier to find files later.



3.) Let's go into the directory using a very common command `cd` --> `change directory`

```
cd Oct22Workshop
```

4.) Make a new file that is empty

```
touch emptyfile.txt
```

`touch` is a bash command that creates an empty file.

??? question "Why would you want an empty file?"

    Some programs require some pre-existing file names to be created.

5.) Make a new file that contains "Hello World"

```
echo "Hello World" > helloworld.txt
```

`echo` is a command that prints the content to the terminal window (sometimes refered to as `print to screen`

5.) Return to your home directory and run `ls`

```
cd
```

```
ls
```

!!! tip
    If you want to speed up the execution of commands, you can copy and paste multiple commands at the same time.

    ```
    cd
    ls
    ```

!!! question
    Please put a green checkmark in your box if you see the new directory when you type `ls` from your home directory).
    
   
## Setting Parameters for Bash Commands


As you start using bash more and more, you will find a mix of files and directories/folders. 
If we want to know which is which, we can add a `parameter` (sometimes referred to as a `flag`)

This is an example of adding a `parameter` without an `argument`.

```
ls -F
```

## Adding Arguments to Bash Commands

An `argument` is a file name or other data that is provided to a command.


```
ls -F Oct22Workshop
```

It is possible to list the files and see their types inside a specific directory by adding the `argument` of the directory name to the `ls` command.

??? question "What do you see when you run the two commands above?"

    Anything with a "/" after it is a directory.  
    Anything with a `*` after it are programs. (we will make a program later)  
    If there's nothing there it's an otherwise unremarkable file (e.g. a data file or picture).

Depending on which terminal you are using, some of the file types may have different colors. 

In our ondemand shell:

Files are white
Directories are blue
Programs (also called `executables`) are green
Compressed files are red (e.g. files that end in .zip or .gzip or .tar)


## Other Useful Parameters for `ls`

Show hidden files

```
ls -a
```

You should see a file called `.bashrc` here. This may be a file we need for troubleshooting your work or where you can make shortcuts or add paths to your login.


Show the `long form` of the list command

```
ls -l
```


To see whether items in a directory are files or directories. `ls -l`
gives a lot more information too, such as the size of the file.


It also shows the permissions of who can read, write or execute a file.


```
drwxrwx--- 2 username05 username05     4096 Jul 18 09:57 JulyWorkshop

```

The first 10 letters in this line indicates the permission settings.


<img width="523" alt="File_Permissions" src="https://user-images.githubusercontent.com/8632603/179539739-75f4edf9-5f5d-4de9-b20c-97abc7869be6.png">


## Getting Help on the Command Line

There are an overwhelming number of possibilities with some of these shell commands, so knowing how to find help on demand is important.

For example, `ls` has a lot of flags that can be used.

```
ls --help
```

This outputs a list of all the ways that `ls` can be altered to find information about your files.


!!! tip "Parameters can be added together in some cases."

    ```
    ls -ltr
    ```

    This can replace `ls -l -t -r`
    `l` is for long form of the list (outputs the permission settings -- something we need to troubleshoot occasionally)
    `t` is to order the files chronologically
    `r` means to reverse the order of the files to put the newest file at the bottom

    This command strings together three flags.

    `ls -l` is list with details
    `ls -t` is sort the list by creation time
    `ls -r` is sort the list in reverse

    For very full directories, this is helpful because it outputs the most recent set of files as the last in the list.


Another way to get help is to use the `man` command. Not every unix installation has this installed, but the Tufts cluster does.

`man` is short for "manual"

!!! note "Navigating a `man` page"

    Use the `spacebar` to scroll through the document.
    Use `q` to leave the manual and go back to the command line prompt.

```
man ls
```

This opens up the manual on the `ls` command. It spells out the meaning of all the parameters in detail.

Most common bash commands have a `man` page that explains it (I wish they had this for emojis....).

Many programs have a help function built in, try adding `--help` or `-h` to see if some helpful information pops up. Sometimes just running the command without any arguments or parameters leads to some usage information or describes the correct command to get help.

For example, if I want to understand the command `tr` - which is used to change a word or character to a new value.

Most programs recognize when you ask for an incorrect parameter, and will tell you how to get more information, as in this example. 
To get help, type the command with the correct parameter.

!!! tip
    For some programs, the `help` function may be `-h`, `--help` 

```
tr -h
```

The shell outputs:

```
tr: invalid option -- 'h'
Try 'tr --help' for more information
```


```
tr --help
```

In this case, a `man` page does exist, so you can get even more direction by typing:

```
man tr
```





