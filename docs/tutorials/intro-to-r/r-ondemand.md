## Setup

Before getting started you will need:

- Account on [Tufts HPC](https://access.tufts.edu/research-cluster-account)
- [VPN](https://access.tufts.edu/vpn) if accessing the HPC from off campus

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `3`
- `Number of cores` : `1`
- `Amount of memory` : `32GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Default`
- `Load Supporting Modules`: `curl/7.47.1 gcc/7.3.0 hdf5/1.10.4 boost/1.63.0-python3 libpng/1.6.37 java/1.8.0_60 libxml2/2.9.10 libiconv/1.16 fftw/3.3.2 gsl/2.6`

Click `Lauch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

## Introduction To RStudio

RStudio is what is known as an Integrated Development Environment or IDE. Here you can write scripts, run R code, use R packages, view plots, and manage projects. This pane is broken up into three panels:

- **The Interactive R console/Terminal (left)**
- **Environment/History/Connections (upper right)**
- **Files/Plots/Packages/Help/Viewer (lower right)**

![](images/rstudio1.png)

# Project Management

Before we dive into R it is worth taking a moment to talk about project management. Often times data analysis is incremental and files build up over time resulting in messy directories:

![](images/messy.png)

Sifting through a non-organized file system can make it difficult to find files, share data/scripts, and identify different versions of scripts. To remedy this, It is reccomended to work within an R Project.

## R Project

To Create a new R project:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `new-project`)
5. `Create Project`

You will notice that your RStudio console switches to this project directory. When you log out of RStudio you can open this project again by clicking the `.Rproj` file in the project directory. 

!!! note
    The paths will be relative to this project directory as a safe guard against referencing data from outside sources. 

## File Organization

When working on a scientific project it is recommended that you put each project in its own directory and give it a name that is descriptive. Similarly, when naming scripts it is recommended that you also name these scripts after the function they are performing. When it comes to file structure within your project try following this folder structure:

- `doc` : folder for text documents associated with the project
- `data` : folder for your raw data/metadata
- `results` : folder for files generated from data/metadata
- `src` : folder for project's custom scripts/programs
- `bin` : folder for outside programs used in project

## Data Principles

When you deal with data treat it as read-only. Working with data files in something like excel can modify your original data source without any record of what was done to it. That being said, often times you will need to do some data cleaning. When you need to significantly modify your data source make a separate folder withing `data` for the `raw_data` and the `cleaned_data`. Also ensure that the scripts you used to clean the data are placed in a separate folder (e.g. `src/data_cleaning_scripts/`). Data that is generated from this raw data should be deposited in your `results` folder and should be treated as disposable. These files should be reproducible from your raw data using your scripts and are good candidate files to cut if you are getting low on storage.

## Script Management

When performing analyses you'll note that some code blocks are useful in multiple scenarios. It is a good idea to store these reusable chunks in a separate folder to use in other analysis scripts. 

## Setting The Directory

When we create a project all the paths are now relative to the project. This helps when you need to specify where a file is. So for example instead of:

```
cluster/home/user/new-project/data/file.txt
```

We only need to specify where things are with respect to our project:

```
data/file.txt
```

We can tell where we are using `getwd()`, so if we were in `new-project`:

```
getwd()
```

```
cluster/home/user/new-project/
```

If we want to specify a **new** base directory we can use `setwd()`:

```
setwd("cluster/home/user/new-project/data")
```

Here we set it to `data` so that if we were to pull that file again, the path would be:

```
file.txt
```

But let's set it back to the project directory for now:

```
setwd("cluster/home/user/new-project")
```
