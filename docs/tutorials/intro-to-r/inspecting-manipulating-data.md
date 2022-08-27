## Importing Data

When importing data we use a few common functions:

* ```read.csv()``` - to read in .csv files or files separated by commas
* ```read.table()``` - to read files separated by delimiters other than commas - like spaces, tabs, semicolons, etc.
* ```openxlsx::read.xlsx()``` - to read excel files

You'll note that ```read.xlsx()``` has the prefix ```openxlsx::```. This is because the ```read. xlsx()``` function is not avaiable with base R. To get this function you will need either need to:

- specify the package that the function comes from:
  
```R
openxlsx::read.xlsx()
```
  
- load the library with the package:
  
```R
library(openxlsx)
read.xlsx()
```

We will now practice inspecting data frames that we will copy over from a shared location. In the `Terminal` tab enter the following command:

```bash
cp /cluster/tufts/bio/tools/training/intro-to-r/data/* data/
```

### read.csv()

When importing `.csv` files you'll need to specify the path to where you're file is located. So if your `.csv` file is in `data/test.csv`, you can download it like so:

```R
read.csv("data/test.csv")
```

We can also extend this to URL's as well:

```R
read.csv(url("https://zenodo.org/api/files/739025d8-5111-476a-9bb9-7f28a200ce8e/linked-ee-dataset-v20220524-QT_2022-07-13-sdev.csv"))
```

### read.table()

Like ```read.csv()```, ```read.table()``` can also import data. The latter function is very useful in that it can download files not delimted (a.k.a separated) by commas. So to open a ".tsv" file (a.k.a a file delimeted by a tab ```"\t"```):

```R
meta <- read.table("data/metadata.tsv",sep="\t",stringsAsFactors=FALSE)
```

You'll notice in the code above that we include the option, ```stringsAsFactors=FALSE```. If this was set to ```TRUE``` it would coerce your character columns into factor columns and this isn't always desired. So here we explicitly say ```stringsAsFactors=FALSE``` to be safe.

### read.xlsx()

While files like the ones mentioned above are popular, so are excel spreadsheets. So it is worth mentioning how to read in excel data as well:

```R
library(openxlsx)      
read.xlsx("data/test.xlsx")
```

Now in excel spreadsheets you may only want to pull out one page or start from a row that isn't the first. To do so you can use:

```R
library(openxlsx)
read.xlsx("data/test.xlsx",sheet=1,startRow = 1,colNames = TRUE,rowNames = FALSE)
```

So here we are pulling: the document "/Documents/test.xlsx", the second sheet, starting from the fifth row, specifying we do have column names, specifying we do not have row names. 

## Inspecting Data

You might have noticed that the only data frame we saved to a variable was the `metadata.tsv` file. We are going to now examine this file:

To get a summary of each column:

```R
summary(meta)
```

```
   SampleID         AntibioticUsage    DaySinceExperimentStart   Genotype         Description           OtuCount     
 Length:9           Length:9           Length:9                Length:9           Length:9           Min.   : 175.0  
 Class :character   Class :character   Class :character        Class :character   Class :character   1st Qu.: 279.0  
 Mode  :character   Mode  :character   Mode  :character        Mode  :character   Mode  :character   Median : 452.0  
                                                                                                     Mean   : 777.8  
                                                                                                     3rd Qu.:1451.0  
                                                                                                     Max.   :1492.0 
```
    

To get the data's class:

```R
class(meta)
```

```
[1] data.frame
```

To get a display of the data's contents:

```R
str(meta)
```

```
'data.frame':	9 obs. of  6 variables:
 $ SampleID               : chr  "WT.unt.1" "WT.unt.2" "WT.unt.3" "WT.unt.7" ...
 $ AntibioticUsage        : chr  "None" "None" "None" "None" ...
 $ DaySinceExperimentStart: chr  "DAY0" "DAY0" "DAY0" "DAY0" ...
 $ Genotype               : chr  "WT" "WT" "WT" "WT" ...
 $ Description            : chr  "16S_WT_unt_1_SRR2627457_1" "16S_WT_unt_2_SRR2627461_1" "16S_WT_unt_3_SRR2627463_1" "16S_WT_unt_7_SRR2627465_1" ...
 $ OtuCount               : int  1174 1474 1492 1451 314 189 279 175 452
```

 
To get the first 6 rows:

```R
head(meta)
```

```
  SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount
1   WT.unt.1            None                    DAY0       WT   16S_WT_unt_1_SRR2627457_1     1174
2   WT.unt.2            None                    DAY0       WT   16S_WT_unt_2_SRR2627461_1     1474
3   WT.unt.3            None                    DAY0       WT   16S_WT_unt_3_SRR2627463_1     1492
4   WT.unt.7            None                    DAY0       WT   16S_WT_unt_7_SRR2627465_1     1451
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189
```
    

To get the last 6 rows:

```R
tail(meta)
```

```
SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount
4   WT.unt.7            None                    DAY0       WT   16S_WT_unt_7_SRR2627465_1     1451
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189
7 WT.day3.15    Streptomycin                    DAY3       WT 16S_WT_day3_15_SRR2628507_1      279
8 WT.day3.14    Streptomycin                    DAY3       WT 16S_WT_day3_14_SRR2627471_1      175
9  WT.day3.9    Streptomycin                    DAY3       WT  16S_WT_day3_9_SRR2628504_1      452
```
    

To get the length of a vector:

```R
length(meta$Genotype)
```

```
[1] 9
```

To get the dimensions of a matrix/data frame:

```R
dim(meta) # answer is given in number of rows, then number of columns
```

```
[1] 9 6
```

To get the number of columns/rows:

```R
ncol(meta)
```

```
[1] 6
```

```R
nrow(meta)
```

```
[1] 9
```

To get your column names:

```R
colnames(meta)
```

```
[1] "SampleID"                "AntibioticUsage"        
[3] "DaySinceExperimentStart" "Genotype"               
[5] "Description"             "OtuCount" 
```

To get your row names:

```R
rownames(meta)
```

```
[1] "1" "2" "3" "4" "5" "6" "7" "8" "9"
```

Now that we know how to import our data and inspect it, we can go ahead and manipulate it!

## Manipulating Data

So now that we have downloaded and inspected our data we can get to manipulating it! So to start, let's talk about accessing parts of your data. To grab the first column in a data frame/matrix you can do so like:

```R
meta[,1]
```

```
[1] "WT.unt.1"   "WT.unt.2"   "WT.unt.3"   "WT.unt.7"   "WT.day3.11"
[6] "WT.day3.13" "WT.day3.15" "WT.day3.14" "WT.day3.9" 
```

To grab the first row:

```R
meta[1,]
```

```
  SampleID AntibioticUsage DaySinceExperimentStart Genotype
1 WT.unt.1            None                    DAY0       WT
                Description OtuCount
1 16S_WT_unt_1_SRR2627457_1     1174
```
    

Now if your data is a data frame you have a special way of accessing coluns with the ```$``` operator:

```R
meta$AntibioticUsage
```

```
[1] "None"         "None"         "None"         "None"         "Streptomycin"
[6] "Streptomycin" "Streptomycin" "Streptomycin" "Streptomycin"
```

This comes in handy for readability. While you can grab your data by column number, it is much easier to read that you are grabbing Sepal Length. To grab mulitple columns/rows, you can do the following for both data frames and matrices:

```R
meta[,c(2,4,6)] # grabbing the 2nd, 4th, and 6th columns
```

```
  AntibioticUsage Genotype OtuCount
1            None       WT     1174
2            None       WT     1474
3            None       WT     1492
4            None       WT     1451
5    Streptomycin       WT      314
6    Streptomycin       WT      189
7    Streptomycin       WT      279
8    Streptomycin       WT      175
9    Streptomycin       WT      452
```

In a data frame, to access columns you can be more specific and specify by column name:

```R
meta[,c("SampleID","Genotype","OtuCount")]
```

```
    SampleID Genotype OtuCount
1   WT.unt.1       WT     1174
2   WT.unt.2       WT     1474
3   WT.unt.3       WT     1492
4   WT.unt.7       WT     1451
5 WT.day3.11       WT      314
6 WT.day3.13       WT      189
7 WT.day3.15       WT      279
8 WT.day3.14       WT      175
9  WT.day3.9       WT      452
```

Now if we wanted to add a new column we could add one like so:

```R
meta$Day <- c(0,0,0,0,3,3,3,3,3) # name of new column comes after the $ sign
meta
```

```
    SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount Day
1   WT.unt.1            None                    DAY0       WT   16S_WT_unt_1_SRR2627457_1     1174   0
2   WT.unt.2            None                    DAY0       WT   16S_WT_unt_2_SRR2627461_1     1474   0
3   WT.unt.3            None                    DAY0       WT   16S_WT_unt_3_SRR2627463_1     1492   0
4   WT.unt.7            None                    DAY0       WT   16S_WT_unt_7_SRR2627465_1     1451   0
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314   3
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189   3
7 WT.day3.15    Streptomycin                    DAY3       WT 16S_WT_day3_15_SRR2628507_1      279   3
8 WT.day3.14    Streptomycin                    DAY3       WT 16S_WT_day3_14_SRR2627471_1      175   3
9  WT.day3.9    Streptomycin                    DAY3       WT  16S_WT_day3_9_SRR2628504_1      452   3
```

## Subsetting Data

To subset our data we need to know a little bit about the different logical operators:

| Operator | Description |
:-------|:-----|
| > | greater than | 
| >= | greater than or equal |
| < | less than |
| <= | less than or equal |
| == | equals | 
| != | not equal |
| & | and |
| \| | or|

Let's go through a few of these!

Subsetting so that we only have rows where the `OtuCount` is greater than 1000:

```R
meta[meta$OtuCount > 1000,]
```

```
  SampleID AntibioticUsage DaySinceExperimentStart Genotype               Description OtuCount Day
1 WT.unt.1            None                    DAY0       WT 16S_WT_unt_1_SRR2627457_1     1174   0
2 WT.unt.2            None                    DAY0       WT 16S_WT_unt_2_SRR2627461_1     1474   0
3 WT.unt.3            None                    DAY0       WT 16S_WT_unt_3_SRR2627463_1     1492   0
4 WT.unt.7            None                    DAY0       WT 16S_WT_unt_7_SRR2627465_1     1451   0
```

Subsetting so that we only have rows where `OtuCount` is less than 400:

```R
meta[meta$OtuCount < 400,]
```

```
    SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount Day
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314   3
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189   3
7 WT.day3.15    Streptomycin                    DAY3       WT 16S_WT_day3_15_SRR2628507_1      279   3
8 WT.day3.14    Streptomycin                    DAY3       WT 16S_WT_day3_14_SRR2627471_1      175   3
```

Subsetting so that we only have rows where the `AntibioticUsage` is equal to `Stretomycin`:

```R
 meta[meta$AntibioticUsage == "Streptomycin",]
```
 
```
    SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount Day
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314   3
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189   3
7 WT.day3.15    Streptomycin                    DAY3       WT 16S_WT_day3_15_SRR2628507_1      279   3
8 WT.day3.14    Streptomycin                    DAY3       WT 16S_WT_day3_14_SRR2627471_1      175   3
9  WT.day3.9    Streptomycin                    DAY3       WT  16S_WT_day3_9_SRR2628504_1      452   3
```

Subsetting so that we only have rows where the `AntibioticUsage` is not equal to `Stretomycin`:

```R
meta[meta$AntibioticUsage != "Streptomycin",]
```
 
```
  SampleID AntibioticUsage DaySinceExperimentStart Genotype               Description OtuCount Day
1 WT.unt.1            None                    DAY0       WT 16S_WT_unt_1_SRR2627457_1     1174   0
2 WT.unt.2            None                    DAY0       WT 16S_WT_unt_2_SRR2627461_1     1474   0
3 WT.unt.3            None                    DAY0       WT 16S_WT_unt_3_SRR2627463_1     1492   0
4 WT.unt.7            None                    DAY0       WT 16S_WT_unt_7_SRR2627465_1     1451   0
```

Subsetting so that we only have rows where the `AntibioticUsage` equals `Steptomycin` or the `OtuCount` is less than `300`:

```R
meta[meta$AntibioticUsage == "Streptomycin" | meta$OtuCount < 300,]
```

```
    SampleID AntibioticUsage DaySinceExperimentStart Genotype                 Description OtuCount Day
5 WT.day3.11    Streptomycin                    DAY3       WT 16S_WT_day3_11_SRR2628505_1      314   3
6 WT.day3.13    Streptomycin                    DAY3       WT 16S_WT_day3_13_SRR2628506_1      189   3
7 WT.day3.15    Streptomycin                    DAY3       WT 16S_WT_day3_15_SRR2628507_1      279   3
8 WT.day3.14    Streptomycin                    DAY3       WT 16S_WT_day3_14_SRR2627471_1      175   3
9  WT.day3.9    Streptomycin                    DAY3       WT  16S_WT_day3_9_SRR2628504_1      452   3
```

## Using Dplyr

When subsetting data we should also mention the R package `dplyr`. This package has functionality to neatly modify data frames using the `%>%` operator to separate your subsetting operations. Let's go through a quick example:

```R
library(dplyr)

meta %>%
  filter(OtuCount < 1400) %>%    # filter rows with OtuCount less than 1400
  select(c(SampleID,AntibioticUsage,Genotype,OtuCount)) %>%   # Select certain rows
  group_by(AntibioticUsage) %>%   # group data by some column
  mutate(HighOtuCount = OtuCount > 1000)   # add a new column 
```

```
# A tibble: 6 × 5
# Groups:   AntibioticUsage [2]
  SampleID   AntibioticUsage Genotype OtuCount HighOtuCount
  <chr>      <chr>           <chr>       <int> <lgl>       
1 WT.unt.1   None            WT           1174 TRUE        
2 WT.day3.11 Streptomycin    WT            314 FALSE       
3 WT.day3.13 Streptomycin    WT            189 FALSE       
4 WT.day3.15 Streptomycin    WT            279 FALSE       
5 WT.day3.14 Streptomycin    WT            175 FALSE       
6 WT.day3.9  Streptomycin    WT            452 FALSE  
```

!!! tip
    For more dplyr data wrangling tips check out the [Data Wrangling with dplyr and tidyr Cheat Sheet](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
 

