## Plotting with Plotly

While there are other plotting libraries, we will focus on `plotly` for the following reasons:

- has the ability to zoom 
- images can be downloaded as `png` files
- select features can highlight features of the plot

Let's make a scatterplot:

```
import plotly.express as px
fig = px.scatter(df,                      # the data we are using
                 x="Day",                 # x axis data
                 y="OtuCount",            # y axis data
                 color='Day',             # how to color our data
                 template="simple_white") # what theme we would like
fig.show()
```

![](images/scatterplot.png)

We can add a trend line as well:

```
import plotly.express as px
fig = px.scatter(df,
                 x="Day",
                 y="OtuCount",
                 color='Day',
                 template="simple_white",
                 trendline="ols")         # add in a trend line
fig.show()
```

![](images/trend-line.png)

