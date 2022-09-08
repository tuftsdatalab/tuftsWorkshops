## Creating and Manipulating Files
----------------------------------

Let's make a file here using a common command "echo" to start creating our file structure.

```
echo "Hello World " > helloworld.txt

```

The `>` in this command tells the command to place the output into the place it is pointing. In this case, it creates the file `helloworld.txt` and puts the phrase `Hello World` into the file. 

#### Helpful Tip
===================

Be careful with redirect.

When using `>` to redirect content into a file, if the filename already exists, it will **overwrite** the file. This means that the original file is gone, and there is no undo in shell.

If you want to add to a file (for example if you are running the same command on several files and extracting a piece of information that you want to put together at the end) you can use another form of redirect `>>`. Using the double redirect will **add** to the file instead of overwriting it.

Which one is used depends on your process. If you are only running a command once, or have an intermediate file in a process that does not need to be retained at the end, then `>` is okay to use.

-------------------------

### Reading File Contents
---------------------------

Let's look inside the file. We have several methods of viewing the content of files that we have created.

A helpful command is `cat`.

```
cat helloworld.txt
```

"cat" will open the entire file, so this is not the best command for long files.

In that case "head" is a good option. Head pulls the top ten lines of the file and prints them to the screen.

```
head helloworld.txt
```

It does not look any different from cat in this case because there is only one line in the file.

A third way to check file contents is by using a program called "less" (or "more").

"less" will open the file interactively, then you can scroll through it and when you are done, push "q" on your keyboard to close the file.

```
less helloworld.txt
```

Press "q" to close the file.

There are many versions of these tools on command line, but "cat", "head" and "less" are very common.


### Copying Files
------------------

Sometimes we have a file that we want to reuse.

When copying within the same directory, make sure to change the name of the file, or the original will be **overwritten**.

When copying to a new directory, the name can stay the same.

This command copies the file within the same directory with a new name. Both files are kept.

```
cp helloworld.txt helloworld1.txt
```
Check this with `ls`

These commands make a new directory, and then copies the file into the new directory with the same name.

```
mkdir helloworld
cp helloworld.txt helloworld
```

Check this with `ls helloworld` (lists the contents of the directory).

### Moving Files
---------------------------

`mv` is an option for renaming files, but also has the potential to **overwrite** existing files.

For example, this command changes the name of the file and removes the original file. If `helloworld2.txt` already existed, it would be replaced.

```
mv helloworld1.txt helloworld2.txt
```

Check this with `ls`

### Removing Files
---------------------------

`rm` and `rmdir` are permanent in shell, so make sure you are ready to delete files.

```
rm helloworld/helloworld.txt
```

Once the directory is empty, we can remove the directory.

```
rmdir helloworld
```

It will throw an error if the directory is not empty.

If you are positive that you want to remove a directory and all the files within it, then add two flags, `-r` for recursive and `-f` for force.

Both commands above could have been replaced with one remove command: `rm -rf helloworld`

Until you are confident with file structure and bash commands, it is a good idea to copy instead of move and to set the interactive flag `-i` on the commands `rm` and `mv` to set up a question that you answer `y` or `n` to before removing.

```
rm -i helloworld/helloworld.txt
```

Generates this question
```
rm: remove regular file ‘helloworld/helloworld.txt’?
```

`mv -i` only generates a question if you are in danger of overwriting an existing file.

For example:

```
cp helloworld.txt helloworld1.txt
mv -i helloworld.txt helloworld1.txt
```
Generates the question:

```
mv: overwrite ‘helloworld1.txt’?
```
