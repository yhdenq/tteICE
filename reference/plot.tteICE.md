# Plot method for 'tteICE' objects

This function plots the estimated potential cumulative incidence
functions or treatment effect curve with pointwise confidence intervals.

## Usage

``` r
# S3 method for class 'tteICE'
plot(
  x,
  type = c("ate", "inc")[1],
  decrease = FALSE,
  conf.int = 0.95,
  xlab = "Time",
  xlim = NULL,
  ylim = NULL,
  plot.configs = list(),
  ...
)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- type:

  Which plot to create: `type="ate"` indicates to plot the estimated
  treatment effects; `type="inc"` indicates to plot the estimated
  cumulative incidence functions (CIFs).

- decrease:

  Corresponds to the argument in
  [`plot_ate`](https://mephas.github.io/tteICE/reference/plot_ate.md)
  and
  [`plot_inc`](https://mephas.github.io/tteICE/reference/plot_inc.md).

- conf.int:

  \#' Confidence level for the pointwise confidence intervals If
  `conf.int = NULL`, no confidence intervals are provided.

- xlab:

  Label for the x-axis.

- xlim:

  A numeric vector of length 2 specifying the limits of the x-axis. If
  `xlim=NULL` (default), the range is determined automatically from the
  data.

- ylim:

  A numeric vector of length 2 giving the limits of the y-axis. If
  `ylim=NULL` (default), the range is determined automatically by the
  type of plot, corresponding to the argument in
  [`plot_ate`](https://mephas.github.io/tteICE/reference/plot_ate.md)
  and
  [`plot_inc`](https://mephas.github.io/tteICE/reference/plot_inc.md).

- plot.configs:

  A named `list` of additional plot configurations. See details in
  [`plot_ate`](https://mephas.github.io/tteICE/reference/plot_ate.md)
  and
  [`plot_inc`](https://mephas.github.io/tteICE/reference/plot_inc.md)

- ...:

  Other arguments in function
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html) or
  function [`curve`](https://rdrr.io/r/graphics/curve.html)

## Value

Plot the results from a tteICE object

## See also

[`plot_ate`](https://mephas.github.io/tteICE/reference/plot_ate.md),
[`plot_inc`](https://mephas.github.io/tteICE/reference/plot_inc.md),
[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.md),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.md),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE.md)

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
bmt$A = A

## simple model fitting and plotting
library(survival)
fit1 = tteICE(Surv(t2,d4,type = "mstate")~A, data=bmt)
plot(fit1, type="ate")

plot(fit1, type="inc")



## plot cumulative incidence functions with p-values
fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
plot(fit2, type="inc", decrease=TRUE, ylim=c(0,1),
     plot.configs=list(show.p.value=TRUE))


## plot treatment effects for semicompeting risk data
fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
plot(fit3, type="ate", ylim=c(-1,1), xlab="time",
     plot.configs=list(col="red"))


```
