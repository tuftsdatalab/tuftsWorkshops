## Setup and Seurat Objects

## Loading Libraries and Data

- We will be working with Single-Cell RNA-seq data in R today. This data is often stored as a Seurat object, which has the following structure:

![](images/single_cell_exp_obj.png)

- Let's start by loading the libraries we need to import and manipulate this object!

```R
# --- Load Libraries -----------------------------------------------------------
.libPaths(c("","/cluster/tufts/bio/tools/R/bio_sup"))
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)

# --- Set color palette --------------------------------------------------------

cols = c("#41ae76","#ee8866","#bebada","#bbcc33","#fdb462",
         "#f768a1","#fa9fb5","#77aadd","darkgray",
         "#cc6677","#882255","#225522","#aa4499","#332288",
         "#009988","#5C61FF","#B87ACF")
# --- Load Data ----------------------------------------------------------------

# start with the day 35 seurat object 
asd_d35 <- readRDS("./results/asd_organoids/asd_d35.rds")
```

## Seurat Objects

```R
# --- Explain the Seurat Object ------------------------------------------------

# what is in this Seurat Object?
asd_d35

# gene names
rownames(asd_d35)

# cell names
colnames(asd_d35)

# what assays do I have?
Seurat::Assays(asd_d35)

# how do I access these assays?
GetAssayData(object = asd_d35, 
             assay = "RNA",
             slot = "counts")

# how do I switch the default assay to be used?
DefaultAssay(asd_d35) <- "RNA"

# how do I access the meta data?
head(asd_d35@meta.data)

# how can I access the dimension reductions?
Embeddings(object = asd_d35, reduction = "pca")
Embeddings(object = asd_d35, reduction = "umap")

# How can I visualize my clustering?
DimPlot(object = asd_d35,
        reduction = "umap")

# what are the identities?
Idents(object = asd_d35)

# how can I change the identities to the cell type?
Idents(object = asd_d35) <- asd_d35$CellType

# how can I see if this changed the identities?
DimPlot(object = asd_d35,
        reduction = "umap")

# meta data in umap
FeaturePlot(object = asd_d35,
        reduction = "umap",
        features =c("percent.mito","percent.ribo"))

# genes in umap
FeaturePlot(object = asd_d35,
            reduction = "umap",
            features ="TOP2A")
```
