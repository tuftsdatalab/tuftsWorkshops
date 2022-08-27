## Functions

So far we have used functions in base python or python modules. But what if we want to create our own? Well here is the general formula to do so:

```py
# load the module we need
import numpy as np 

# define a function to return geometric mean
def geometric_mean(values):
    return np.exp(np.mean(np.log(values)))
   
# call our function
geometric_mean(coverage)                 
    
```

```
209.88855396892262
```

Here we see that we use `def` to define the function and `return` to specify what value you'd like to return. We then call our function and use our `coverage` variable as our set of values. The geometric mean of the set of values in `coverage` are then returned.

### Function Documentation

To clarify the purpose of your function you can add a **multiline string** to your function using three quotes `'''`:

```py
# define a function to return geometric mean
def geometric_mean(values):
    ''' This function takes a list of
    values and returns the geometric mean 
    of those values'''
    return np.exp(np.mean(np.log(values)))
```

This multiline string is also accessible when we run our function through the `help` function:

```py
help(geometric_mean)
```

```py
Help on function geometric_mean in module __main__:

geometric_mean(values)
    This function takes a list of
    values and returns the geometric mean 
    of those values
```

## Variable Scope

Variable naming can be difficult and sometimes variable names might need to be reused. Normally, when we use a variable name over again, we change the value of that variable. However, if we assign the same variable in and outside a function the values do not get overwritten:

```py
x = 45

def print_x():
    x = 30
    return x
```

```
x
```

```
45
```

```py
print_x()
```

```
30
```
