# Plot the results of the search

This function plots the results of the search. It uses ggplot2 to create
a line plot of the number of publications over time.

## Usage

``` r
plot_time_interval(x, interval = "Weeks")
```

## Arguments

- x:

  The results group_by_time_interval

- interval:

  The time interval to plot. Default is "Weeks". Choose "Weeks",
  "Months", "Years".

## Value

A ggplot object.
