# Introduction To RStudio OnDemand

!!! example "Prerequisites"
    - [Request an account](http://research.uit.tufts.edu/) on the Tufts HPC Cluster
    - Connect to the [VPN](https://access.tufts.edu/vpn)
    - Please be sure to have followed the instructions on the [setup page](../setup.md)

---

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `3`
- `Number of cores` : `1`
- `Amount of memory` : `4GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `GET THIS FROM DELILAH`
- `Load Supporting Modules`: `curl/7.47.1 gcc/7.3.0 hdf5/1.10.4 boost/1.63.0-python3 libpng/1.6.37 java/1.8.0_60 libxml2/2.9.10 libiconv/1.16 fftw/3.3.2 gsl/2.6`

Click `Launch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

??? question "Have you logged into the cluster?"
    - Yes (put up a green check mark in zoom)
    - No (raise hand in zoom)
    
---

## Introduction To RStudio

RStudio is what is known as an Integrated Development Environment or IDE. Here you can write scripts, run R code, use R packages, view plots, and manage projects. This pane is broken up into three panels:

- **The Interactive R console/Terminal (left)**
- **Environment/History/Connections (upper right)**
- **Files/Plots/Packages/Help/Viewer (lower right)**

![](images/rstudio1.png)

## Project Management

Before we dive into R it is worth taking a moment to talk about project management. Often times data analysis is incremental and files build up over time resulting in messy directories:

![](images/messy.png)

Sifting through a non-organized file system can make it difficult to find files, share data/scripts, and identify different versions of scripts. To remedy this, It is reccomended to work within an R Project.

## R Project

To Create a new R project:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `R-Practice`)
5. `Create Project`

You will notice that your RStudio console switches to this project directory. When you log out of RStudio you can open this project again by clicking the `.Rproj` file in the project directory. 

!!! note
    The paths will be relative to this project directory as a safe guard against referencing data from outside sources. 

??? question "Have you created the project?"
    - Yes (put up a green check mark in zoom)
    - No (raise hand in zoom)

--- 

## File Organization

Let's start by adding some directories to store the data and scripts we will use. In the `Terminal` tab enter the following commands:

```
mkdir data
mkdir scripts
mkdir results
```

Now we will need some data to play with so let's copy some over:

```
cp /cluster/tufts/bio/tools/training/intro-to-r/data/* data/
```
---

## Data Principles

- Treat data as read-only
- Store raw data separately from cleaned data if you do need to manipulate it
- Ensure scripts to clean data are kept in a separate `scripts` folder
- Treat reproducible results as disposable

!!! tip
    Result files are good candidate files to cut if you are getting low on storage.

---
