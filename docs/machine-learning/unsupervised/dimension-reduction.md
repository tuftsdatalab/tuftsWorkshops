# Introduction To Dimension Reduction

!!! attention
    Please be sure to have followed the instructions on the [setup page](../setup.md)
    
- Process of reducing the number of variables to a set of principal values where variation in your data becomes apparent. Here is an example with three dimensions:

<figure markdown>
  ![Dimension Reduction Example](images/example-dim-red.png){ width="500" }
  <figcaption>Dimension Reduction Example</figcaption>
</figure>

- Here we see that most of the variation is visible along the x-y axes
- So what are the advantages:

  - simplification
  - denoising
  - variable selection
  - visualization

## Principal Component Analysis (PCA)

- PCA works by summarizing a set of continuous (quantitative) multivariate (multiple variable) data into a set of **linearly uncorrelated** variables called principal components.

### Pros

- can be used on large data
- can be used with sparse data
- preserves the structure (reproducible)

### Cons

- if one variable is on a different scale (like kg instead of g) it can bias the results. So ensure data is on one scale!
- points can be crowded with large numbers of observations and reveal no pattern
- susceptible to outliers

Let's try this in code! First we will need to do some preprocessing:

=== "R"

    ``` R
    ## load our libraries via our library path
    .libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0"))
    library(tidyverse)
    library(FactoMineR)
    library(factoextra)
    library(ggplot2)
    library(missMDA)
    library(patchwork)
    
    ## load counts/meta data
    counts <- read.csv(
      file="data/gbm_cptac_2021/data_mrna_seq_fpkm.txt",
      header = T,
      sep = "\t")

    meta <- read.csv(
      file = "data/gbm_cptac_2021/data_clinical_patient.txt",
      skip=4,
      header = T,
      sep = "\t"
    )
    
    ## ensure patient IDs match 
    ## patient IDs in counts data
    meta <- meta %>%
      mutate(PATIENT_ID = gsub("-",".",meta$PATIENT_ID)) %>%
      column_to_rownames("PATIENT_ID")
    ```

=== "Python"

    ``` py
    # still in development - sorry!
    ```

## References

- [RPubs](https://rpubs.com/Saskia/520216)
- [STHDA](http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization)
