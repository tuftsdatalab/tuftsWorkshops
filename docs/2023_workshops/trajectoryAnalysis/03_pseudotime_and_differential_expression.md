## Pseudotime and Differential Expression

## Pseudotime

monocle introduced the concept of pseudotime which they define as: ["Pseudotime is a measure of how much progress an individual cell has made through a process such as cell differentiation."](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/) We will assess "progress" by a cell's differentiation status. So we will manually chose the starting point to be the cycling progenitors: 

```R
# --- Pseudotime ---------------------------------------------------------------
# manually chose the starting point to be the cycling progenitors
cds_2 = order_cells(cds_2) 
```

INSERT IMAGE HERE

Now let's visualize pseudotime in our UMAP plot to understand what is early pseudotime (low values - dark purple) and late pseudotime (high values - yellow):

```R
# plot pseudo time after choosing root nodes 
# (defined as the bottom in the cycling progenitors cell type)
plot_cells(cds_2,
           show_trajectory_graph = T,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 0.5)

```

INSERT IMAGE HERE

## Differential Expression

```R
# --- Differential Expression --------------------------------------------------

# run the differential expression testing 
gene_fits <- fit_models(sub_cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs %>% filter(term != "(Intercept)") %>%
  dplyr::select(
    c(gene_short_name,
      term,
      q_value,
      estimate)) %>%
  filter(q_value<0.05)
```


