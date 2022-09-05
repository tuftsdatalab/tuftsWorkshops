# Introduction To Logistic Regression

!!! attention
    Please be sure to have followed the instructions on the [setup page](../setup.md)
    
## Logistic Regression

A logistic regression model attempts to classify, so for example can IDH1 gene expression to predict smoking status? So how do we do this? We could fit a linear model to our data but we would end up with something like this:

![](images/linear-v-logistic.png)

Here we see that if we used a linear model, we'd end up predicting values that are not our two catetories - smoker or non-smoker (and on the graph 0 or 1). So how do we fit our data when it is binary? We use a sigmoid function:

$$ p(X) = \frac{ e^{\beta_{0} + \beta_{1}X} }{1 + e^{\beta_{0} + \beta_{1}X} } $$

- $p(X)$ : probability of smoking status given IDH1 gene expression
- $X$ : IDH1 gene expression
- $beta_{0}$ : y intercept
- $beta_{1}$ : slope of our line

This sigmoid function creates our S-shaped curve! However we'd like our probability to be linear relationship with X. So we can manipulate this equation to get:

$$ \frac{p(X)}{1 - p(X)} = e^{\beta_{0} + \beta_{1}X}$$

Where $ \frac{p(X)}{1 - p(X)} $ is known as the **odds ratio** and this can range from $0$ to $\infty$. However, again this doesn't match our 0 to 1 scale, so we take the log of both sides to get:

$$ log(\frac{p(X)}{1 - p(X)}) = \beta_{0} + \beta_{1}X $$

Here we get $log(\frac{p(X)}{1 - p(X)})$ or the **logit function** - where a one unity increase in $X$ increases $p(X)$ by $\beta_{0}$. 

## References
