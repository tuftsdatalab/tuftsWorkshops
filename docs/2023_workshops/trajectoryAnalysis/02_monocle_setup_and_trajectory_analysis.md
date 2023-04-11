## Running Monocle3

There are typically four main commands when running Monocle3:

!!! info "[Monocle3 Workflow](https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/)"

    ![](images/monocle_workflow.png)

- `preprocess_cds()`: normalizes the data by log and size factor to address depth differences and calculates a lower dimensional space that will be used as the input for further dimensionality reduction like tSNE and UMAP.
- `reduce_dimension()`: reduces the dimensionality of the data using algorithms like UMAP or tSNE. Trajectories will be calculated through this space.
- `cluster_cells()`: clusters the cells using Louvain/Leiden community detection, and returns a cell_data_set with internally stored cluster assignments. These cluster assignments can then be assigned to cell types given that cells in a cluster are likely to be the same cell type as cells of the same type have similar gene expression patterns.
- `learn_graph()`: constructs the trajectory through clusters in a lower dimensional space to ["learn the sequence of gene expression changes each cell must go through as part of a dynamic biological process"](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/)

Here we will only run through the first three steps:

```R
# --- Running Monocle3 ---------------------------------------------------------

# let's run the monocle workflow:
# use the same # of top PCs as used for clustering
cds <- preprocess_cds(cds, num_dim = 27) 

# reduce the dimensions 
cds <- reduce_dimension(cds,umap.fast_sgd=TRUE)

# force few partitions with partition q-value set higher
cds <- cluster_cells(cds, partition_qval = 0.5) 

# let's take another look at our cell data set object
cds
```

!!! info "output"

    ```R
    class: cell_data_set 
    dim: 3233 33313 
    metadata(2): cds_version citations
    assays(1): counts
    rownames(3233): RP11-54O7.1 AURKAIP1 ... MT-ND5 MT-CYB
    rowData names(1): gene_short_name
    colnames(33313): 1_AAACCCAAGCACGGAT 1_AAACCCAAGTCATACC ... 6_TTTGTTGCATTACGGT 6_TTTGTTGTCAAACCTG
    colData names(11): orig.ident nCount_RNA ... org Size_Factor
    {=={reducedDimNames(2): PCA UMAP==}
    altExpNames(0):
    ```
Here we would like to highlight that after running our pipeline, we now have two `reducedDimNames` slots in our object! We can access them with the `reducedDims` function:

```R
# let's examine our new `reducedDimNames` slot!
head(reducedDims(cds)$UMAP)
```

!!! info "output"

    ```R
                             [,1]       [,2]
    1_AAACCCAAGCACGGAT  2.5821856 -5.6855853
    1_AAACCCAAGTCATACC -9.1280608 -6.5626800
    1_AAACCCAGTAGGAGTC  0.4445358  6.8874695
    1_AAACCCATCACTCTTA -9.6992268 -3.9201496
    1_AAACCCATCTTGAACG  5.5063716 -0.8642332
    1_AAACGAAGTGGCAACA -9.8187041  2.8391122
    ```
To understand why we only ran through the first three steps, we should examine how our cells are distributed in our dimension reduced UMAP plot:

```R
# let's inspect our cell types, clusters and partitions!
p1 <- monocle3::plot_cells(cds,
                           color_cells_by = "CellType", 
                           show_trajectory_graph = FALSE,
                           cell_size = 0.5)
p2 <- monocle3::plot_cells(cds, 
                           color_cells_by = "cluster",
                           show_trajectory_graph = FALSE,
                           cell_size = 0.5)
p3 <- monocle3::plot_cells(cds,
                           color_cells_by = "partition", 
                           show_trajectory_graph = FALSE,
                           cell_size = 0.5)

patchwork::wrap_plots(p1, p2, p3)

```

INSERT IMAGE HERE


## Subsetting Our Data

Above we can see that we have multiple clusters (usually representing our cells) and multiple partitions (usually representing groups of different cells). When Monocle3 calculates it's trajectory it will typically do so through one of these partitions. So we will subset our data to just grab the partition that contains the Cycling Progenitors, Newborn PNs, and Newborn DL PNs. By doing this we can assess how gene expression changes during cell differentiation from a Cycling Progenitors to a Newborn DL PN:

```R
# --- Subsetting our Data ------------------------------------------------------
# when we examine a trajectory in monocle3 it is useful to look at one
# partition as you are examining how gene expression changes between clusters
# in some group
cds_2 = choose_cells(cds)

# re-run the monocle3 workflow on our subset data:
# use the same # of top PCs as used for clustering
cds_2 <- preprocess_cds(cds_2, num_dim = 27) 
# reduce the dimensions 
cds_2 <- reduce_dimension(cds_2,umap.fast_sgd=TRUE)
# force few partitions with partition q-value set higher
cds_2 <- cluster_cells(cds_2, partition_qval = 0.5) 
# additionally we will run learn_graph to calculate trajectories on our subset!
cds_2 <- learn_graph(cds_2)

# now that we hav subset out data, re-run our monocle workflow, and calculated 
# trajectories, let's visualize them!
p5 <- plot_cells(cds_2, 
                 color_cells_by = "CellType",
                 show_trajectory_graph = T,
                 cell_size = 0.5)

p5 <- plot_cells(cds_2, 
                 color_cells_by = "cluster",
                 show_trajectory_graph = T,
                 cell_size = 0.5)

p6 <- plot_cells(cds_2, 
                 color_cells_by = "partition",
                 show_trajectory_graph = T,
                 cell_size = 0.5)

patchwork::wrap_plots(p4,p5,p6)
```

INSERT IMAGE HERE
