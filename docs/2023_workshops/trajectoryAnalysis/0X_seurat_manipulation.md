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

# --- Explain the Seurat Object ------------------------------------------------

# what is in this Seurat Object?
seur

# gene names
rownames(seur)

# cell names
colnames(seur)

# what assays do I have?
Seurat::Assays(seur)

# how do I access these assays?
GetAssayData(object = seur, 
             assay = "RNA",
             slot = "counts")

# how do I switch the default assay to be used?
DefaultAssay(seur) <- "RNA"

# how do I access the meta data?
head(seur@meta.data)

# how can I access the dimension reductions?
Embeddings(object = seur, reduction = "pca")
Embeddings(object = seur, reduction = "umap")

# How can I visualize my clustering?
DimPlot(object = seur,
        reduction = "umap")

# what are the identities?
Idents(object = seur)

# how can I change the identities to the cell type?
Idents(object = seur) <- seur$CellType

# how can I see if this changed the identities?
DimPlot(object = seur,
        reduction = "umap")

# meta data in umap
FeaturePlot(object = seur,
            reduction = "umap",
            features =c("percent.mito","percent.ribo"))

# genes in umap
FeaturePlot(object = seur,
            reduction = "umap",
            features ="TOP2A")

```
