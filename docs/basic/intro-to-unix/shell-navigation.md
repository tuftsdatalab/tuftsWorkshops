## Navigating in the Shell
===============================

We are going to make a place to work for this workshop.

The following command makes a new directory.

```
mkdir JulyWorkshop

```
----------------------

#### Helpful Tip
===================

When nameing files and directories, avoid spaces and special characters except underscores ("_").

**Spelling** and **Capitalization** are literal in unix, be careful when making and using files to be consistent in your process. This will make it easier to find files later.

---------------------

You can check that the new directory was created by repeating the list command.

```
ls
```

A directory is like a desk drawer. We create them to store files that relate to each other mostly.

When creating directories and filenames it is helpful to put some information about the project and the date of activity.


<img width="786" alt="File_Folder_Structure" src="https://user-images.githubusercontent.com/8632603/179539866-ecd6e880-f468-4151-bbaa-149f52c328b4.png">

## Absolute and Relative Paths
-------------------------------

Let's go into our directory and look around.

Another command you'll find yourself using a lot is `cd`, which stands
for 'change directory'.  Try typing::

```
cd JulyWorkshop

```
and then

```
pwd
```

You should now see something like this:

```
/cluster/home/username01/JulyWorkshop
```

This is an example of an **Absolute Path**.

It gives an address for where you are located on the cluster, much like a postal address that defines where you are in several layers (e.g. /country/state/city/street/specific_house.

<img width="821" alt="tufts_root_path" src="https://user-images.githubusercontent.com/8632603/179759502-549b38b0-4957-4105-aee3-8bca4271bf7b.png">

If you want to go back to the directory that is in the level above our current file, another common shortcut used in bahs is `..`.


```
cd ..
```

`..` is a reference to a **RELATIVE PATH**

```
pwd
```

You should be back in your home directory.

```
/cluster/home/username01/
```

If you want to go back to the directory that you just left, type this command.

```
cd -
```
Then find your location.

```
pwd
```

You should be back in the directory you came from.

```
/cluster/home/username01/JulyWorkshop
```

A **RELATIVE PATH* means that the command only works from the relative location that you are in.

`cd ..` and `cd -` are examples of relative path commands.

This can get confusing if you are moving around a lot in your directories or sending commands to SLURM, so the alternative method to navigating around the cluster is using an **ABSOLUTE PATH**.


Note: If you ever type `cd` without a word behind it, it will send you back to your home directory.

Your home directory is not all the way back at the root, it is set within the cluster as `/cluster/home/username01/`.

You can make sure that you are in the right directory by using the command `cd` with the absolute path.

```
cd /cluster/home/username01/JulyWorkshop
```

#### Helpful Tip
=================
Many commands in bash can be used with the ABSOLUTE PATH.

```
ls /cluster/home/username01/JulyWorkshop
```

Using an absolute path to find files in a directory is helpful for checking for outputs from SLURM jobs when they are running.