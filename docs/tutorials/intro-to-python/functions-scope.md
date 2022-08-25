## Functions

So far we have used functions in base python or python modules. But what if we want to create our own? Well here is the general formula to do so:

```
# load the module we need
import numpy as np 

# define a function to return geometric mean
def geometric_mean(values):
    return np.exp(np.mean(np.log(values)))
   
# apply our function
geometric_mean(coverage)                 
    
```

```
209.88855396892262
```

Here we see that we use `def` to define the function and `return` to specify what value you'd like to return. 

### Function Documentation

To clarify the purpose of your function you can add a **multiline string** to your function using three quotes `'''`:

```
# define a function to return geometric mean
def geometric_mean(values):
    ''' This function takes a list of
    values and returns the geometric mean 
    of those values'''
    return np.exp(np.mean(np.log(values)))
```

This multiline string is also accessible when we run our function through the `help` function:

```
help(geometric_mean)
```

```
Help on function geometric_mean in module __main__:

geometric_mean(values)
    This function takes a list of
    values and returns the geometric mean 
    of those values
```
