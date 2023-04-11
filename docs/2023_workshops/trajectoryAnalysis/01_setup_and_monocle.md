
!!! example "Prerequisites"
    - [Request an account](http://research.uit.tufts.edu/) on the Tufts HPC Cluster
        - Note if you signed up for the Introduction to Single-Cell RNA-Seq Time Series and Trajectory Analysis workshop this will have been already taken care of for you!
    - Connect to the [VPN](https://access.tufts.edu/vpn) if off campus
    

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/){:target="_blank" rel="noopener"} and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

!!! info "OnDemand Layout"

    ![](images/ondemand_layout_pic.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `5`
- `Number of cores` : `1`
- `Amount of memory` : `16GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Bioinformatics Workshops`
    - **NOTE: This reservation closed on April 26th 2023, use `Default` if running through the materials after that date.**
- `Load Supporting Modules`: `boost/1.63.0-python3 java/1.8.0_60 gsl/2.6`

Click `Launch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

??? question "Are you connected to RStudio?"
    - Yes (put up a green check mark in zoom)
    - No (raise hand in zoom)

## Today's Data

Today we will be working with data from  [Paulson et al. 2022](https://www.nature.com/articles/s41586-021-04358-6) which found cell-type-specific neurodevelopmental abnormalities that were shared across ASD risk genes. To this end they leveraged organoid single-cell RNA-seq data to investigate these abnormalities:

!!! info "[Paulson et al. 2022](https://www.nature.com/articles/s41586-021-04358-6)"

    ![](images/asd_figure_1.png)
    

## Data & Scripts

To copy over this data we will enter the following command into the console:

```R
file.copy(from="/cluster/tufts/bio/tools/training/trajectory_analysis",to="./", recursive = TRUE)
```

## Project Setup

Now that we have copied over the directory for today's workshop we are going to use this folder to create a new R project. R projects are great for managing analyses in a containerized way. To create an R project from within our `trajectory_analysis` directory we will:

- Go to `File` > `New Project`
- `Existing Directory`
- Browse for the `trajectory_analysis` folder
- Click `Create Project`

Now that we have our data and scripts copied, let's navigate to our `scripts` folder and open up "trajectory_analysis.Rmd".

## Loading Libraries and Data

Before we analyze our single-cell RNA-Seq data we will need to load the libraries needed for our analysis:

```R
# --- Load Libraries -----------------------------------------------------------
LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("",LIB))
.libPaths()
library(Seurat)
library(monocle3)
library(clusterProfiler)
library(patchwork)
library(tidyverse)

# --- Set color palette --------------------------------------------------------

cols = c("#41ae76","#ee8866","#bebada","#bbcc33","#fdb462",
         "#f768a1","#fa9fb5","#77aadd","darkgray",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#5C61FF","#B87ACF")
         
# --- Load Data ----------------------------------------------------------------

# start with the day 35 seurat object 
seur <- readRDS("./results/asd_organoids/suv420h1_mito210_d35_sub.rds")
```

## Monocle3 Cell Data Objects

- We will be working with single-cell RNA-seq data in R today. This data is often stored as a Seurat object if you are performing differential expression testing. To understand how to work with a Seurat object check out our [quick tutorial](0X_seurat_manipulation.md). Today, we will be performing trajectory analysis using the R package Monocle3. Monocle3 stores single-cell RNA-seq data as a cell data set object, which has the following structure:

!!! info "Monocle3's Cell Data Set Object"

    ![](images/cell_data_set.png)

- Let's start by loading the libraries we need to import and manipulate this object!

```R
# isolate the gene names
fd = data.frame("gene_short_name" = rownames(GetAssayData(seur)))
rownames(fd) = rownames(GetAssayData(seur))

# create a cell data cell object that monocle can use for trajectory
# analysis
cds <- new_cell_data_set(GetAssayData(seur),
                         cell_metadata = seur@meta.data,
                         gene_metadata = fd)
cds
```

!!! info "output"

    ```R
    class: cell_data_set 
    {==dim: 3233 27666==} 
    metadata(1): cds_version citations
    assays(1): counts
    {==rownames(3233): RP11-54O7.1 AURKAIP1 ... MT-ND5 MT-CYB==}
    rowData names(1): gene_short_name
    {==colnames(27666): 1_AAACCCAAGCACGGAT 1_AAACCCAAGTCATACC ... 4_TTTGGTTCAAGCCTGC 5_AGACAGGTCACAAGGG==}
    {==colData names(11): orig.ident nCount_RNA ... org Size_Factor==}
    {==reducedDimNames(0): ==}
    altExpNames(0):
    ```
    
- Here we will highlight that we have X rows and Y columns, our rownames are gene names, our column names are the cell names, we have one assay (`counts`), we have 11 columns of meta data under `colData`, and we have no dimension reductions under `reducedDimNames`.
- Let's investage a few helpful functions that can help access these data:

```R
# access the gene names
rownames(cds)[1:5]
```

!!! info "output"

    ```R
    [1] "RP11-54O7.1" "AURKAIP1"    "SSU72"       "C1orf233"    "MIB2" 
    ```

```R
# access the cell names
colnames(cds)[1:5]
```


!!! info "output"

    ```R
    [1] "1_AAACCCAAGCACGGAT" "1_AAACCCAAGTCATACC" "1_AAACCCAGTAGGAGTC" "1_AAACCCATCACTCTTA" "1_AAACCCATCTTGAACG"
    ```
    
```R
# access the feature data
head(rowData(cds))
```

!!! info "output"

    ```R
    DataFrame with 6 rows and 1 column
                gene_short_name
                    <character>
    RP11-54O7.1     RP11-54O7.1
    AURKAIP1           AURKAIP1
    SSU72                 SSU72
    C1orf233           C1orf233
    MIB2                   MIB2
    MMP23B               MMP23B
    ```

```R
# access the meta data
head(colData(cds))
```

!!! info "output"

    ```R
    [DataFrame with 6 rows and 11 columns
                       orig.ident nCount_RNA nFeature_RNA percent.mito percent.ribo            CellType       treat seurat_clusters         dataset      org Size_Factor
                         <factor>  <numeric>    <integer>    <numeric>    <numeric>         <character> <character>        <factor>     <character> <factor>   <numeric>
    1_AAACCCAAGCACGGAT          1       3086          544    0.1148179     0.474840                 aRG          wt              4  SUV_Mito210_d35        1    0.575785
    1_AAACCCAAGTCATACC          1       6499         1364    0.0745991     0.107398         Subcortical          wt              15 SUV_Mito210_d35        1    1.221490
    1_AAACCCAGTAGGAGTC          1      12168         1618    0.0678902     0.179736 Cycling Progenitors          wt              9  SUV_Mito210_d35        1    1.332651
    1_AAACCCATCACTCTTA          1       5393         1237    0.1456875     0.124290       Cajal-Retzius          wt              16 SUV_Mito210_d35        1    1.141904
    1_AAACCCATCTTGAACG          1       4951         1054    0.0581777     0.246685                 aRG          wt              0  SUV_Mito210_d35        1    0.988533
    1_AAACGAAGTGGCAACA          1       7081         1429    0.0966265     0.108369         Newborn PNs          wt              3  SUV_Mito210_d35        1    1.289798
    ```

```R
# access the assay data
head(assay(cds)[,1:3])
```

!!! info "output"

    ```R
    6 x 3 sparse Matrix of class "dgCMatrix"
                1_AAACCCAAGCACGGAT 1_AAACCCAAGTCATACC 1_AAACCCAGTAGGAGTC
    RP11-54O7.1                  .           .                  .       
    AURKAIP1                     .           5.337212           5.443008
    SSU72                        .           .                  4.069609
    C1orf233                     .           4.648863           .       
    MIB2                         .           4.648863           .       
    MMP23B                       .           .                  .       
    ```

## Matching Cell Types Between Conditions

Often times, you will be examining cell types from multiple conditions (i.e. control v. treatment). It is good practice to ensure you have similar numbers of cells between conditions as to not bias your analysis. In our data we have far more mutant cells than we do wild-type cells. So we will be downsampling our mutant cells to match the number of wild-type cells.

!!! info "Note about downsampling"

    If you want your work to be reproducible, it is recommended that you set a seed so that you are consistently sampling the same set of cells.
    
```R
# --- Downsampling -------------------------------------------------------------

#look at which one has fewer cells and downsample the other to that number
table(colData(cds)$treat) 
```


!!! info "output"

    ```R
    mut    wt 
    19480 13833 
    ```

```R
# the wild type cell type has the least number of cells
# let's grab all the wild type cells
wtCells = rownames(colData(cds)[colData(cds)$treat=="wt",])

# now we will need to downsample our mutant cells
# to match the number of wild-type cells
# set a seed to make sure the same cells are sampled
set.seed(123)
mutCells = sample(
  rownames(colData(cds)[colData(cds)$treat=="mut",]), 
  nrow(colData(cds)[colData(cds)$treat=="wt",]),replace = F)

# we will save our cell data set object with the matched
# numbers of cell types
cds = cds[,c(wtCells,mutCells)]

# now how do our cell type proportions look?
table(colData(cds)$treat)
```

!!! info "output"

    ```R
    mut    wt 
    13833 13833 
    ```
