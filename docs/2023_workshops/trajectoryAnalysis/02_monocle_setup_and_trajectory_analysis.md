## Monocle Setup and Trajectory Analysis


## Monocle Setup

```R
# --- Monocle3 Setup -----------------------------------------------------------
# We can use the SeuratWrappers package to convert our Seurat object to a 
# Monocle object. We will then be using Monocle's clustering technique to create
# clusters and partitions (groups of clusters that can be analyzed)
cds <- as.cell_data_set(asd_d35)

# note that we can use other resolutions if you aren't finding great clusters:
# 1e-5, 1e-3, 1e-1, and 0
cds <- cluster_cells(cds, resolution=1e-5)

# let's inspect our clusters and our partitions!
p1 <- plot_cells(cds, 
                 color_cells_by = "cluster",
                 show_trajectory_graph = FALSE,
                 cell_size = 0.5)
p2 <- plot_cells(cds,
                 color_cells_by = "partition", 
                 show_trajectory_graph = FALSE,
                 cell_size = 0.5)
wrap_plots(p1, p2)

# when we examine a trajectory in monocle3 it is useful to look at one
# partition as you are examining how gene expression changes between clusters
# in some group
int_sub <- subset(as.Seurat(cds, assay = NULL), monocle3_clusters == 2)
cds <- as.cell_data_set(int_sub)
```

## Trajectory Analysis

```R
# --- Trajectory Analysis ------------------------------------------------------

# Trajectory analysis in Monocle is a way to assess the relationship between 
# groups of cells based on gene expression changes. First we will need to use 
# the learn graph function to sketch out trajectories
cds <- learn_graph(cds,
                   use_partition = TRUE, 
                   verbose = FALSE)

# we can now visualize this trajectory over cell types
plot_cells(cds,
           color_cells_by = "CellType",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 4,
           cell_size = 0.5)

# let's visualize the distribution of wild type and mutant cells

# isolate only the wild type cells
wtCells = rownames(colData(cds)[colData(cds)$treat=="wt",])

# isolate only the mutant type cells
mutCells = rownames(colData(cds)[colData(cds)$treat=="mut",])

# plot the distribution of cells in each condition
wt <- plot_cells(cds[,wtCells],
                 color_cells_by = "CellType",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 group_label_size = 4,
                 cell_size = 0.5)+
  labs(title="Wild Type Cells")

mut <- plot_cells(cds[,mutCells],
                 color_cells_by = "CellType",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 group_label_size = 4,
                 cell_size = 0.5)+
  labs(title="SUV420H1 Mutant Cells")

wrap_plots(wt,mut)
```
