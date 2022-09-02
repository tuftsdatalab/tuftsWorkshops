# Introduction To Linear Regression

!!! danger
    Please be sure to have followed the instructions on the [setup page](../setup.md)
    
## Linear Regression

In a regression model assesses the relationship between two quantitative variables by fitting a line to the data. 
Using our glioblastoma data, we will assess the relationship between IDH1 gene expression (a gene commonly mutated in this type of cancer) 
and TMB score (a measure of mutational burden). You can use a regression model to determine:

- how strong the relationship between these two variables is
- the value of the dependent variable given the independent variable

A linear model follows the following formula:

$$ 
y = \beta_0 + \beta_1 X + \epsilon
$$

- $y$: dependent variable
- $\beta_0$: intercept (where $y$ = 0)
- $\beta_1$: regression coefficient or slope
- $X$: independent variable
- $\epsilon$: error or our estimate (what is the variation in our regression coefficient)

This formula describes the best fit line for our data that tries to minimizes our error $\epsilon$:

![](images/linear-regression-demo.png)

Let's start by loading our data:

=== "R"

    - Sed sagittis eleifend rutrum
    - Donec vitae suscipit est
    - Nulla tempor lobortis orci

=== "Python"

    
    ```py
    ## Import our libraries
    ## Import our data set
    import pandas as pd
    import numpy as np
    from matplotlib import pyplot as plt
    import statsmodels.api as sm

    counts = pd.read_csv(
        'data/gbm_cptac_2021/data_mrna_seq_fpkm.txt' ,
        sep = '\t')
    meta = pd.read_csv(
        'data/gbm_cptac_2021/data_clinical_sample.txt' , 
        sep = '\t',
        skiprows=4)
    ## You might notice our data is on two 
    ## drastically different scales
    ## We will normalize our data

    def NormalizeData(data):
        normalized = (data - np.min(data)) / (np.max(data) - np.min(data))
        df = pd.DataFrame(normalized)
        return normalized

    norm = NormalizeData(merged)
    ## Let's visualize our data
    plt.figure(figsize=(12, 6))
    plt.plot(norm['IDH1'], norm['TMB_NONSYNONYMOUS'], 'o') 
    plt.xlabel('IDH1')
    plt.ylabel('TMB_NONSYNONYMOUS')

    plt.show()
    ## fit our linear regression model
    tmb = norm['TMB_NONSYNONYMOUS']
    idh1 = norm['IDH1'].astype(float)
    model = sm.OLS(tmb,idh1).fit()
    model.summary()
    ```

## Model Results

## Assumptions

## References

- [scribbr](https://www.scribbr.com/statistics/simple-linear-regression/)
- [STHDA](http://www.sthda.com/english/articles/40-regression-analysis/165-linear-regression-essentials-in-r/)
