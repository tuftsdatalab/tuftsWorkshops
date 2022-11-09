!!! example "Prerequisites"
    - [Request an account](http://research.uit.tufts.edu/) on the Tufts HPC Cluster
    - Connect to the [VPN](https://access.tufts.edu/vpn) if off campus

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/){:target="_blank" rel="noopener"} and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `4`
- `Number of cores` : `1`
- `Amount of memory` : `16GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Bioinformatics Workshop`
- `Load Supporting Modules`: `boost/1.63.0-python3 java/1.8.0_60 gsl/2.6`

Click `Lauch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

## Copy Over the Workshop Folder Into Your Home Directory

Copy this code chunk by clicking on the "copy" icon that is in the upper right corner.
Paste the code chunk into the console window in R studio (lower left window)

```R

file.copy(from="/cluster/tufts/bio/tools/training/metagenomics/Metagenomics2022",to="~/", recursive = TRUE)

```

## Project Setup

We are going to create a new project in the existing directory that we just copied over:

1. Go to `File` > `New Project`
2. Click "Don't Save" if it asks about saving the workspace.
3. `Existing Directory`
4. Choose the directory we just copied over (e.g. `~/Metagenomics2022`)

!!! question "What if my screen goes blank?"
    
    Don't worry if the screen goes blank for a moment. It just means that a shiny new workspace is being created so you can start on your project.


## Open the R Notebook

Now that we have our data and scripts copied, let's navigate to our scripts folder and open up `Metagenomics.Rmd`.



## What to Expect from Today

Here is a rendered version of the notebook that displays only the text and the code and not the output.

[Today's Notebook](http://htmlpreview.github.io/?https://github.com/tuftsdatalab/tuftsWorkshops/blob/main/docs/omics/intro-metagenomics/Metagenomics.nb.html)



