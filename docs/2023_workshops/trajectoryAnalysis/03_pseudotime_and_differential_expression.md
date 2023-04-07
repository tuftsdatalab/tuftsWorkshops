## Pseudotime and Differential Expression

## Pseudotime

```R
# --- Pseudotime ---------------------------------------------------------------

# monocle introduced the concept of pseudotime which they define as:
# "Pseudotime is a measure of how much progress an individual cell has made 
# through a process such as cell differentiation."
# to d
# manually chosen to be the bottom in cycling progenitors
cds = order_cells(cds) 

# plot pseudo time after choosing root nodes 
# (defined as the bottom in the cycling progenitors cell type)
plot_cells(cds,
           show_trajectory_graph = T,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 0.5)
```

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


