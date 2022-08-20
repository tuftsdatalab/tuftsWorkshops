## Inspecting/ Manipulating Data

## Importing Data

When importing data we use a few common functions:

* ```read.csv()``` - to read in .csv files or files separated by commas
* ```read.table()``` - to read files separated by delimiters other than commas - like spaces, tabs, semicolons, etc.
* ```openxlsx::read.xlsx()``` - to read excel files

You'll note that ```read.xlsx()``` has the prefix ```openxlsx::```. This is because the ```read. xlsx()``` function is not avaiable with base R. To get this function you will need either need to:

- specify the package that the function comes from:
  
```
openxlsx::read.xlsx()
```
  
- load the library with the package:
  
```
library(openxlsx)
read.xlsx()
```

We will now practice inspecting data frames that we will copy over from a shared location. In the `Terminal` tab enter the following command:

```
cp /cluster/tufts/bio/tools/training/intro-to-r/data/* data/
```

### read.csv()

When importing `.csv` files you'll need to specify the path to where you're file is located. So if your `.csv` file is in `data/test.csv`, you can download it like so:

```
read.csv("data/test.csv")
```

We can also extend this to URL's as well:

```
read.csv(url("https://zenodo.org/api/files/739025d8-5111-476a-9bb9-7f28a200ce8e/linked-ee-dataset-v20220524-QT_2022-07-13-sdev.csv"))
```

### read.table()

Like ```read.csv()```, ```read.table()``` can also import data. The latter function is very useful in that it can download files not delimted (a.k.a separated) by commas. So to open a ".tsv" file (a.k.a a file delimeted by a tab ```"\t"```):

```
meta <- read.table("data/metadata.tsv",sep="\t",stringsAsFactors=FALSE)
```

You'll notice in the code above that we include the option, ```stringsAsFactors=FALSE```. If this was set to ```TRUE``` it would coerce your character columns into factor columns and this isn't always desired. So here we explicitly say ```stringsAsFactors=FALSE``` to be safe.

### read.xlsx()

While files like the ones mentioned above are popular, so are excel spreadsheets. So it is worth mentioning how to read in excel data as well:

```
library(openxlsx)      
read.xlsx("data/test.xlsx")
```

Now in excel spreadsheets you may only want to pull out one page or start from a row that isn't the first. To do so you can use:

```
library(openxlsx)
read.xlsx("data/test.xlsx",sheet=1,startRow = 1,colNames = TRUE,rowNames = FALSE)
```

So here we are pulling: the document "/Documents/test.xlsx", the second sheet, starting from the fifth row, specifying we do have column names, specifying we do not have row names. 

## Inspecting Data

You will have noticed that the only data frame we saved to a variable was the `metadata.tsv` file. We are going to now examine this file:

To get a summary of each column:

```
summary(meta)

SampleID         AntibioticUsage    DaySinceExperimentStart
 Length:9           Length:9           Length:9               
 Class :character   Class :character   Class :character       
 Mode  :character   Mode  :character   Mode  :character       
                                                              
                                                              
                                                              
   Genotype         Description           OtuCount     
 Length:9           Length:9           Min.   : 175.0  
 Class :character   Class :character   1st Qu.: 279.0  
 Mode  :character   Mode  :character   Median : 452.0  
                                       Mean   : 777.8  
                                       3rd Qu.:1451.0  
                                       Max.   :1492.0  
```
    

To get the data's class:

```
class(meta)
```

> [1] data.frame

To get a display of the data's contents:

```
str(meta)

'data.frame':	9 obs. of  6 variables:
 $ SampleID               : chr  "WT.unt.1" "WT.unt.2" "WT.unt.3" "WT.unt.7" ...
 $ AntibioticUsage        : chr  "None" "None" "None" "None" ...
 $ DaySinceExperimentStart: chr  "DAY0" "DAY0" "DAY0" "DAY0" ...
 $ Genotype               : chr  "WT" "WT" "WT" "WT" ...
 $ Description            : chr  "16S_WT_unt_1_SRR2627457_1" "16S_WT_unt_2_SRR2627461_1" "16S_WT_unt_3_SRR2627463_1" "16S_WT_unt_7_SRR2627465_1" ...
 $ OtuCount               : int  1174 1474 1492 1451 314 189 279 175 452
```

 
To get the first 6 rows:

    head(iris)
 
    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    1          5.1         3.5          1.4         0.2  setosa
    2          4.9         3.0          1.4         0.2  setosa
    3          4.7         3.2          1.3         0.2  setosa
    4          4.6         3.1          1.5         0.2  setosa
    5          5.0         3.6          1.4         0.2  setosa
    6          5.4         3.9          1.7         0.4  setosa
    

To get the last 6 rows:

    tail(iris)

    Sepal.Length Sepal.Width Petal.Length Petal.Width   Species
    145          6.7         3.3          5.7         2.5 virginica
    146          6.7         3.0          5.2         2.3 virginica
    147          6.3         2.5          5.0         1.9 virginica
    148          6.5         3.0          5.2         2.0 virginica
    149          6.2         3.4          5.4         2.3 virginica
    150          5.9         3.0          5.1         1.8 virginica
    

To get the length of a vector:

    length(iris$Sepal.Length)
    
    150

To get the dimensions of a matrix/data frame:

    dim(iris)


    150 5 #(so this would be 150 rows and 5 columns)

To get the number of columns/rows:

    ncol(iris)
    
    5

    nrow(iris)

    150

To get your column names:

    colnames(iris)

    Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"

To get your row names:

    rownames(iris)

    "1"   "2"   "3"   "4"   "5" ...

Now that we know how to import our data and inspect it, we can go ahead and manipulate it!
