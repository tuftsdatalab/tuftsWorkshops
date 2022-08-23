# Libraries

Libraries are collections of functions called modules that can be imported and used in your script. Let's use the `math` library to grab constants:

```
import math
math.e
```

```
2.718281828459045
```

Now how about functions:

```
math.log2(25)
```

```
4.643856189774724
```

!!! tip
    If you ever need assistance with a library, try using the `help()` function to grab more information (e.g. `help(math)`).
    
## Importing Parts of Libraries & Using Aliases

Sometimes you'll only need a few things from a library. To grab just those few things use the following approach:

```
from math import log2, e
math.log2(25)
```

```
4.643856189774724
```

Now sometimes the name of a library is just too long to continuously type out. For this we can use an **alias**

```
from math import log2 as l2
math.l2(25)
```

```
4.643856189774724
```

Here we abbreviate `log2` from the `math` package to `l2`.

## Data Frames 

In data analysis we often work with tabular data, or two dimensional data with columns and rows. Columns will typically contain the same type of data and rows will be one sample with different observations. We commonly read in tabular data using the `pandas` module:

```
dataframe <- 
```
