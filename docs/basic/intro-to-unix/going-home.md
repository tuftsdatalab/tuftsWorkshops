
## Going Home
---------------------------


Sometimes we get lost, so it is useful to know a few ways to get back to where you started.

```
cd

```

This command returns you to your home directory. Check by typing this command.

```
pwd
```

Other options for going back to your home directory:

```
cd ~
```

```
cd $HOME
```

When lost in the file structure, going home is a good place to start.

### Find and Tree

File structures can get complicated quickly.

Two tools to understand where your files are that can help are `find` and `tree`.

From your home directory, you can find your file named `helloworld.txt` by typing the following:

```
cd
find . -name helloworld.txt
```

From the home directory, the answer should look like:

```
./JulyWorkshop/helloworld.txt
```

`.` in this case is another RELATIVE path direction that indicates "from this directoy that I am in currently". Note that the answer is given in the RELATIVE path format, starting with `.` = here.

It is also possible to provide an ABSOLUTE path to this command.

```
find /cluster/home/username01 -name helloworld.txt
```
This command will work from anywhere in the cluster. Note that the answer is given in the ABSOLUTE path format.

```
/cluster/home/arhode05/JulyWorkshop/helloworld.txt
```

Another helpful bash command for finding files is `tree`.

```
tree
```

This outputs your directory structure with lines that indicate the tree-like branches of your file structure.

This could be very messy if you already have a lot of files in a directory, so limit the level by adding a flag.

```
tree -L 2
```

This just shows the top two levels of the file structure.
