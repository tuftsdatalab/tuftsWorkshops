## Variables

In Python, we store values using names called variables. We can assign a variable with an `=` sign:

```
max_coverage = 6000
minCoverage = 35
antibiotic = 'Streptomycin'
antibiotic2 = 'Penicillin'
```
You will note a few things about variables:

- can incorporate letters, digits and underscores
- cannot start with a digit
- these are case sensitive

Variables, once we assign them to some value, can be passed into functions to accomplish certain tasks. Let's use the simplist use case, the `print()` function:

```
print('The maximum coverage is ', max_coverage)
```
```
The maximum coverage is  6000
```

!!! note
    variables are available to use between blocks. However, the order in which you run blocks matters so make sure to run your code blocks in order!
    
## Data Types

- `integer`: a positive/negative whole number (34, -675)
- `float`: a floating point number (4.67, -2034.67)
- `string`: a character string written with either single or double quotes ('Streptomycin', "antibiotic")
- `bool`: a TRUE/FALSE value

So you have a variable, how do you determine the type? Well we can use the `type()` function:

```
type(max_coverage)
```

```
int
```

If you want to convert between data types you can specify with the following functions:

- `int()`: to convert to an integer
- `float()`: to convert to a floating point number
- `str()`: to convert to a string


## Calculations

You can use Python like a calculator using the following symbols:

|Operator	|Name	|Example	|
|-|-|-|
|+	|Addition	|x + y|
|-	|Subtraction	|x - y|
|*	|Multiplication	|x * y|
|/	|Division	|x / y|
|%	|Modulus	|x % y|
|**	|Exponentiation	|x ** y|
|//	|Floor division	|x // y|

Let's try an few example:

```
35 / 7 - 5 + 4 * 4 + 2**2
```

```
20.0
```

We note that Python calculations follow the order of operations when performing a calculation. We should also bring up two non-standard operations that you may or may not be familiar with: Modulus and Floor division. Modulus is the remainder after division so:

```
7 % 2
```

```
1
```

Floor division is a division operation for which you round the result down to a whole number:

```
7 // 2
```

```
3
```

## Strings & Operators

You can use `+` and `*` with string data as well to add and multiply, take for instance:

```
antibiotic + antibiotic
```

```
'StreptomycinStreptomycin'
```

```
antibiotic * 4
```

```
'StreptomycinStreptomycinStreptomycinStreptomycin'
```

## Indexing

Unlike the other data types, strings have lengths. We can use the `len()` function to  check how long  a string is:

```
print(antibiotic)
len(antibiotic)
```

```
'Streptomycin'
12
```

We can slice strings if needed to! However, the letters you are grabbing are **0-indexed** meaning that the first letter is letter 0, the second letter is letter 1, and so on:

```
antibiotic[0]
```

```
'S'
```

```
antibiotic[1]
```

```
't'
```

We can grab more letters using the format `[start:stop]`:

```
antibiotic[1:5]
```

```
'trep'
```
