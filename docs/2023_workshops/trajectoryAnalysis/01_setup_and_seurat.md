## Setup and Seurat Objects

## Loading Libraries and Data

```R
# --- Load Libraries -----------------------------------------------------------
#LIB='/cluster/tufts/bio/tools/R_libs/4.0.0'
#LIB='/cluster/tufts/hpc/tools/R/4.0.0/'
.libPaths(c("","/cluster/home/jlaird01/R/x86_64-pc-linux-gnu-library/4.0/"))
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