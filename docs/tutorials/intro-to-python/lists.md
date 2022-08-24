# Lists

A data frame is not the only way to store data, we can also create lists of values which can be the same data type or different data types. Here is an example:

```
coverage = [200, 34, 900, 423, 98, 789]
```

### Grabbing List Values

We can also grab these values by their index, which again are **zero-indexed** (meaning they start at zero). Here is an example of grabbing the 3rd item in the list:

```
coverage[2]
```

```
900
```
## Adding/Deleting Values

We can change items in the list as well by assigning them to different values:

```
coverage.append(542)
coverage
```

```
[200, 34, 300, 423, 98, 789, 542]
```
