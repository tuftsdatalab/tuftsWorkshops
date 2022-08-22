## Setup

Before getting started you will need:

- Account on [Tufts HPC](https://access.tufts.edu/research-cluster-account)
- [VPN](https://access.tufts.edu/vpn) if accessing the HPC from off campus

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > JupyterLab` and you will see a form to fill out to request compute resources to use JupyterLab on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `3`
- `Number of cores` : `1`
- `Amount of memory` : `32GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Default`
- `Load Supporting Modules`: `curl/7.47.1 gcc/7.3.0 hdf5/1.10.4 boost/1.63.0-python3 libpng/1.6.37 java/1.8.0_60 libxml2/2.9.10 libiconv/1.16 fftw/3.3.2 gsl/2.6`

Click `Lauch` and wait until your session is ready. Click `Connect To JupyterLab Server`, and you will notice a new window will pop up with JupyterLab. 

## Introduction to JupyterLab

Jupyterlab is a web-based user interface to run Python code and is not a traditional Integrated Development Environment (IDE) where you create scripts via some text editor and then submit directly to command line. JupyterLab has several advantages, including being able to run code in chunks, annotating code with links, and displaying figures right next to code! For this reason, JupyterLab is a robust tool for script development/data analysis. When you open JupyterLab you will notice:

- **Left Sidebar**: containing your file browser, list of running kernels/terminals, table of contents, extension manager
- **Main Work Area**: containing options for file/window types to open (ipykernels, terminal environments, text files, markdown files, and python files)

![](images/jupyterlab.png)

We are going to start by opening up a `.ipynb` file by clicking `Notebook Python 3 (ipykernel)`. These are not python scripts, but notebook files that contain code but also text, links and images. These files can easily be converted to a python script (file ending in `.py`) by going to:

- `File`
- `Download as`
- `Python (.py)`

For now let's work in the Jupyter notebook (`.ipynb` file)!

## Code Vs. Markdown

You will notice when you open up your notebook that you are working in blocks:

![](images/blocks.png)

These blocks can either be:

- **raw blocks:** raw data that can be converted into HTML/Latex formats
- **code blocks:** python code that can be run in chunks
- **markdown blocks:** a plain text format that can render links, lists, and images like what you might find on a website

Here we will focus on code blocks to run chunks of python code, and markdown blocks which can add in images, links, etc. to annotate our code.

## Markdown Basics

**markdown code:**

```
- list item 1
- list item 2
```
**output:**
- list item 1
- list item 2

**markdown code:**

```
1. numbered list item 1
2. numbered list item 2
```
**output:**
1. numbered list item 1
2. numbered list item 2

**markdown code:**

```
# Level 1 Heading
## Level 2 Heading
```
**output:**
# Level 1 Heading
## Level 2 Heading

**markdown code:**

```
[google link](https://www.google.com/)
```
**output:**
[google link](https://www.google.com/)

