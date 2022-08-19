# Introduction To R OnDemand

## Setup

Before getting started you will need:

- Account on [Tufts HPC](https://access.tufts.edu/research-cluster-account)
- [VPN](https://access.tufts.edu/vpn) if accessing the HPC from off campus

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/) and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `3`
- `Number of cores` : `1`
- `Amount of memory` : `32GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Default`
- `Load Supporting Modules`: `curl/7.47.1 gcc/7.3.0 hdf5/1.10.4 boost/1.63.0-python3 libpng/1.6.37 java/1.8.0_60 libxml2/2.9.10 libiconv/1.16 fftw/3.3.2 gsl/2.6`

Click `Lauch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

## Introduction To RStudio

RStudio is what is known as an Integrated Development Environment or IDE. Here you can write scripts, run R code, use R packages, view plots, and manage projects. This pane is broken up into three panels:

- **The Interactive R console/Terminal (left)**
- **Environment/History/Connections (upper right)**
- **Files/Plots/Packages/Help/Viewer (lower right)**

![](images/rstudio1.png)

Now we will create an R script. R commands can be entered into the console, but saving these commands in a script will allow us to rerun these commands at a later date. To create an R script we will need to either:

- Go to `File > New File > R script`
- Click the `New File` icon and select R script

![](images/newFile.png)

## Running R Code

When running R code you have a few options:

  Running One Line/Chunk:
  
  - Put your cursor at the beginning of the line of code and hit `Ctrl + Enter` on Windows or  &#8984; + `Enter` on MacOSX.
    
  - Highlight the line/chunk of code and hit `Ctrl + Enter` or &#8984; + `Enter`.
    
  Running The Entire Script:
  
  - Clicking `Source` at the top of the script window.
    
## Calculations

Let's try running some R code! R can be used to run all sorts of calculations just like a calculator:

![](images/line-of-code.png)

You will notice that we ran this code in the script window but you can see the output in the console. 
When we work with calculations it is useful to remember the order of operations - here are their equivalents in R:

- **Parentheses:** `(`, `)`
- **Exponents:** `^` or `**`
- **Multiply:** `*`
- **Divide:** `/`
- **Add:** `+`
- **Subtract:** `-`

Let's look at some examples:

```
10 * 3^3
```
> [1] 270

```
(400 / 10) * (4e2) # 4e2 is the same as 4^2
```
> [1] 16000

You'll notice that in the last equation we added words after a `#` and the equation still ran. This is what is known as a comment, where everything after the `#` is not registered as R code. Commenting is immensely valuable for giving your code context so that you and whoever else reads it knows the purpose of a given chunk of code.

Additionally there are functions built in R to perform mathematical calculations:
```
abs(10) # absolute value
```
> [1] 10

```
sqrt(25) # square root
```
>[1] 5

```
log(10) # natural logarithm
```

>[1] 2.302585

```
log10(10) # log base 10
```
> [1] 1

## Comparisons

R can also be used to make comparisons. Here we note the operators used to do so:

- **Equals:** `==`
- **Does Not Equal:** `!=`
- **Less Than Or Equal** `<=`
- **Greater Than Or Equal** `>=`
- **Greater Than** `>`
- **Less Than* `<`

```
2 == 2
```
> [1] TRUE

```
2 != 2
```
> [1] FALSE

```
3 <= 10
```
> [1] TRUE


!!! note

    Unless the number is an integer, do not use `==` to compare. This is due to the fact that the decimal value may appear the same 
in R but from a machine level the two values can be very different. 

## Variables & Vectors

Dealing with values can be cumbersome. In R, values can be assigned to words using the `<-` operator:

```
x <- 35 # assigning a value of 35
x
x <- 40 # changing value to 40
x
```
> [1] 35

> [1] 40

You'll notice that we initially assigned `x` to a value of `35` and then updated value to `40`. This is important to keep in mind because the last value assigned to `x` will be kept. Variables can I have a combination lowercase letters, uppercase letters, underscores and periods:

```
value <- 40
biggerValue <- 45
even_bigger_value <- 50
biggest.value <- 55
```
```
value
biggerValue
even_bigger_value
biggest.value
```
> [1] 40

> [1] 45

> [1] 50

> [1] 55

!!! note

    Take note that the spelling needs to be consistent to call the variable correctly.

We can also assign a series of values in a specific order to a variable to create what is called a **vector**:

```
someVector <- 5:10
someVector
```
> [1]  5  6  7  8  9 10

## Environment

As you may have noticed we have been assigning variables and they have been added to your `Environment` window:

![](images/environment.png)

If you would like to declutter your environment, you have a few options:

- You can use the `rm()` function to remove which ever variables you'd like. To remove more than one just put a comma between variable names.
- You can clear all variables by clicking the broom icon:

![](images/remove-all.png)


## R Packages

Aside from the base functions there are thousands of custom fuctions which are bundled in R packages. We can access these functions installing them and loading them. On the Tufts HPC, libraries of packages are available. To access them you will need to specify where these packages are with the `.libPaths()` function:

```
.libPaths("","/cluster/tufts/hpc/R/4.0.0")
```
You'll can see what packagews are available in the`Packages` window:

![](images/packages.png)

To load a package you can use the `library()` function:

```
library(ggplot2)
```
