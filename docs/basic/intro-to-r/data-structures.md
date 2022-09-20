# Data Structures

## Data Types

So far we have only dealt with numeric values. However, we are not limited to just numbers we can store:

* ```numeric``` - numeric values that can contain whole numbers and decimals
* ```character``` - text value that is made a text value by adding quotes. So for example ```1 2 3``` is a numeric data, but ```"1" "2" "3"``` is character data
* ```integer``` - limited to just whole numbers, but will take up less memory than numeric data. We can specify an integer by adding `L` to a number (e.g. `1L` or `3L`)
* ```logical``` - These are boolean values so ```TRUE```/```T``` or ```FALSE```/```F```.
* ```complex``` - complex number such as ```1+6i```

??? question "What data type do you think count data would be?"
    - numeric
    - character
    - integer
    - logical
    - complex

## Data Structures

So we have all this lovely data to play with and in R we typically organize in a few ways:

### Vectors

Vectors are collections of data, like a:

  collection of numbers - ```c(1,2,3)```
  collection of characters -  ```c("1","2","3")```
  collection of logical values - ```c(TRUE,FALSE,TRUE)```

!!! tip
    It should be noted that a vector needs to be a collection of the **same type** of data. You will also note that each list is separated by commas and surrounded by ```c()```. This is necessary to create vectors so make sure to remember the ```c()```!

### Factors

Factors can be used to store categorical data and can be created like this:

  ```R
  size <- c("small", "medium", "small", "large", "medium")
  size
  ```
  
  ```
  [1] "small"  "medium" "small"  "large"  "medium"
  ```
  
  ```R
  size <- factor(size,levels = c("small","medium","large"))
  size
  ```
  
  ```
  [1] small  medium small  large  medium
  
  Levels: small medium large
  ```
  
Now we have turned this character vector into a factor vector! These will come in handy when we start breaking down data by category.

### Matrices

A matrix can be created by combining vectors of the **same length and same data type**. They are used frequently when performing operations on numeric data but can include other data types. In R we can create a matrix with the `matrix()` function:

```R
matrix(data=1:9,nrow = 3,ncol=3)
```
```
     [,1] [,2] [,3]
[1,]    1    4    7
[2,]    2    5    8
[3,]    3    6    9
```

Here we take a vector and specify how many columns and how many rows we'd like. You'll also note that we have multiple arguments in our function:

- `data` : for our input data
- `nrow` : number of rows
- `ncol` : number of columns

!!! note
    While it is true that vectors need to be the same length, they can have different numbers of data points. So for instance:
    
    ```R
    vec1 <- c(1,2,3)
    vec2 <- c(7,NA,4)
    ```
    These vectors are the same length, BUT `vec1` has 3 data points and `vec2` has 2 data points with an `NA` is a filler.

### Data Frames

Data frames are also collections of vectors of the **same length**. However, they do not need to be the same data type. Here we create a data.frame with the `data.frame()` function:

```R
data.frame(
characters=c("past","present","future"),
numbers=c(1,2,3),
logical=c(TRUE,FALSE,TRUE),
integer=c(1L,2L,3L)
)
```
```
  characters numbers logical integer
1       past       1    TRUE       1
2    present       2   FALSE       2
3     future       3    TRUE       3
```

!!! note 
    The inputs here are not arguments specific to the `data.frame()` function. Instead, we are adding columns and the name before the `=` sign is our column name

### Lists

Lists are collections of data that **do not** need to be the same type or length. We can create lists with the `list()` function:

```R
list(
data.frame=data.frame(numbers=1:3,characters=c("past","present","future")),
numbers=1:5,
characters=c("past","present","future")
)
```
```
$data.frame
  numbers characters
1       1       past
2       2    present
3       3     future

$numbers
[1] 1 2 3 4 5

$characters
[1] "past"    "present" "future" 
```

!!! done "Time for a 10 minute break!"
