# Introduction To Linear Regression

!!! danger
    Please be sure to have followed the instructions on the [setup page](setup.md)
    
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

$y$: dependent variable
$\beta_0$: intercept (where $y$ = 0)
$\beta_1$: regression coefficient or slope
$X$: independent variable
$\epsilon$: error or our estimate (what is the variation in our regression coefficient)

This formula describes the best fit line for our data that tries to minimizes our error $\epsilon$:

![](images/linear-regression-demo.png)

Let's try creating a model!

=== "R"

    - Sed sagittis eleifend rutrum
    - Donec vitae suscipit est
    - Nulla tempor lobortis orci

=== "Python"

    1. Sed sagittis eleifend rutrum
    2. Donec vitae suscipit est
    3. Nulla tempor lobortis orci

## Model Results

## Assumptions

## References

- [scribbr](https://www.scribbr.com/statistics/simple-linear-regression/)
