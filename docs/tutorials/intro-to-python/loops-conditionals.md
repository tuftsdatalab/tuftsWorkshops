## Loops

Loops perform some operation on a value in a set of values. Let's go through an example using our `coverage` list from the previous note:

```
for i in coverage:
    print(i)
```

```
200
34
300
98
789
542
```

Here we see that `i` is a substitute for some value in the sequence provided - in this case 200, `34, 300, 98, 789, 542`. 

### Nested Loops

Loops can also **nested** where a loop is placed inside a loop:

```
for i in [1,2]:
    for j in coverage:
        print(j*2)
```

```
200
34
300
98
789
542
400
68
600
196
1578
1084
```

Here we move through the loop and for every value in the first list (`[1,2]`), Then for each pass of the first loop we move through values the second list (`[200, 34, 300, 98, 789, 542]`). Finally for each value `i` we then multiply by each value `j`. 

### Pass Statement

If you want a placeholder for your loop, meaning no operation is performed, use the `pass` statement:

```
for i in coverage:
    pass
```

## Conditionals

If we were interested in performing some operation on a value only if a condition is met, we can use an `if` statement:

```
for i in coverage:
    if i > 500:
        print(i)
    else:
        pass
```

```
789
542
```

Here we use the comparison operators we mentioned in the [Libraries & Data Frames Topic Note](libraries-data-frames.md) to only print values in `coverage` if they are larger than `500`.

### Multiple Conditionals

To perform operations based on multiple conditions you can add in `elif` statements:

```
for i in coverage:
    if i > 500:
        print(i)
    elif i < 500:
        print('This value is less than 500')
        
```

```
This value is less than 500
This value is less than 500
This value is less than 500
This value is less than 500
789
542
```
