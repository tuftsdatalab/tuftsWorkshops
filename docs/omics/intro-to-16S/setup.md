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

## Project Setup

We are going to create a new project to begin:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `intro-to-16S`)
5. `Create Project`

In the `Terminal` tab we will set up our project space:

```bash
mkdir data
mkdir results
mkdir scripts
```

We will copy over our sample data:

```bash
cp -r /cluster/tufts/bio/tools/training/microbiome16S/subsampled/* ./data/
```
## Script Setup

Now to get started we will need to setup a script:

- Go to `File > New File > R script`
- Click the `New File` icon and select R script

In the R script, start by loading the libraries we need:

```R
#LIB='/cluster/tufts/bio/tools/R_libs/4.0.0'
LIB='/cluster/home/jlaird01/R/x86_64-pc-linux-gnu-library/4.0/'
.libPaths(c(LIB))
library(dada2)
library(phyloseq)
library(ggplot2)
library(DESeq2)
```


