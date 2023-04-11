## Running Monocle3

There are typically four main commands when running Monocle3:

![](images/monocle_workflow.png)

- `preprocess_cds()`
- `reduce_dimension()`
- `cluster_cells()`
- `learn_graph()`

```R
# --- Running Monocle3 ---------------------------------------------------------

# let's run the monocle workflow:
# use the same # of top PCs as used for clustering
cds <- preprocess_cds(cds, num_dim = 27) 

# reduce the dimensions 
cds <- reduce_dimension(cds,umap.fast_sgd=TRUE)

# force few partitions with partition q-value set higher
cds <- cluster_cells(cds, partition_qval = 0.5) 
```
