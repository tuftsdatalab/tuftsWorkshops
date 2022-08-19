# Functions and Flow

## Functions

Functions are operations we can perform on our various data structures to get some result. We typically like to make functions modular so they perform one specific task and not whole pipelines. Here is the general format for a function:

```
functionName <- function(x){
  result <- operation(x)
  return(result)
}
```

So here we see that we assign some operation to a name, here it we just call it ```functionName```. Then the function takes an input, ```x```. Inside the function our result is obtained by doing some operation on our input. Finally we then use ```return()``` to return that result. Let's try making a function that will square the input:

```
squareInput <- function(x){
  result <- x * x
  return(result)
}

squareInput(5)
> 25
```

!!! note
   Without `return()` in the function, R will return the last variable in the function. So you can leave out `return()`, however it is best to specify what you are returning for clarity.
   
## Additional Arguements

Functions can have as little or as many arguements as needed. So in the example above we used one arguement, `x`. Let's try using more than one arguement:

```
item_in_vector_func <- function(vector, item){
  item_in_vector <- item %in% vector
  return(item_in_vector)
}
item_in_vector_func(vector = c(1,2,3), item=2)
```

> [1] TRUE

Here we used a new operator, `%in%`, which does exactly what it sounds like - it checks whether some value is in another set of values. We should also note, you can specify values in functions:

```
item_in_vector_func <- function(vector=1:10, item=NULL){
  item_in_vector <- item %in% vector
  return(item_in_vector)
}
```

Here we specify default values for the `item_in_vector()` function. Unless we change them when we call the function, these values will remain. So when you call a function from a package it's a good idea to check what the default values are.

## Control Flow

Now what if we don't want to perform an operation until a condition is met? For this we need an if/else statement:

```
x <- 3

if (x == 10){
  print("x equals 10")
} else{
  print("x does not equal 10")
}
```

> [1] "x does not equal 10"

Now what if we wanted to include multiple conditions?

```
x <- 3

if (x == 10){
  print("x equals 10")
} else if (x > 2){
  print("x is greater than 2")
} else if (x > 1){
  print("x is greater than 1")
} else{
  print("x does not equal  10")
}
```

> [1] "x is greater than 2"

Here notice something, x meets 2 of the conditions:

- `x > 2`
- `x > 1`

However, we note that the conditional statement is broken when `x` meets the first condition in the order above. 

## Loops

Operations can be repeated with a `for` loop:

```
for (i in 1:5){
  print(i*i)
}
```

```
[1] 1
[1] 4
[1] 9
[1] 16
[1] 25
```
Here we see that `i` is a substitute for some value in the sequence provided - in this case `1,2,3,4,5`. We can also nest a loop inside a loop like so:

```
for (i in 1:3){
  for(j in 3:5){
    print(i*j)
  }
}
```
```
[1] 3
[1] 4
[1] 5
[1] 6
[1] 8
[1] 10
[1] 9
[1] 12
[1] 15
```

You'll notice that for each value i was multiplied by each value j. So:

```
1*3
1*4
1*5
2*3
2*5
2*6
3*3
3*4
3*5
```
