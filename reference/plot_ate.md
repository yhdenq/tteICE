# Plot estimated treatment effects

This function plots the estimated treatment effect, defined as the
difference in potential cumulative incidences under treated and control
groups, along with pointwise confidence intervals.

## Usage

``` r
plot_ate(
  fit,
  decrease = FALSE,
  conf.int = 0.95,
  xlab = "Time",
  ylim = c(-1, 1),
  xlim = NULL,
  plot.configs = list(ylab = NULL, main = NULL, lty = 1, lwd = 2, col = "black",
    add.null.line = TRUE, null.line.lty = 2, ci.lty = 5, ci.lwd = 1.5, ci.col =
    "darkgrey"),
  ...
)
```

## Arguments

- fit:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- decrease:

  A logical value indicating the type of curve difference to display. If
  `decrease = FALSE` (default), the difference in cumulative incidence
  functions (CIFs) is plotted. If `decrease = TRUE`, the difference in
  survival functions is plotted instead.

- conf.int:

  Confidence level for the pointwise confidence intervals If
  `conf.int = NULL`, no confidence intervals are provided.

- xlab:

  Label for the x-axis.

- ylim:

  A numeric vector of length 2 specifying the limits of the y-axis.
  Defaults to `ylim = c(-1, 1)`.

- xlim:

  A numeric vector of length 2 specifying the limits of the x-axis. If
  `xlim = NULL` (default), the limits are determined automatically from
  the data.

- plot.configs:

  A named `list` of additional plot configurations. Common entries
  include:

  - `ylab`: character, label for the y-axis (default: `ylab=NULL`, use
    the default label).

  - `main`: character, title for the plot (default: `main=NULL`, use the
    default label).

  - `lty`: line type for effect curve (default: `lty=1`).

  - `lwd`: line width for effect curve (default: `lwd=2`).

  - `col`: line color for effect curve (default: `col="black"`).

  - `add.null.line`: logical, whether to draw a horizontal line at 0
    (default: `add.null.line=TRUE`, add the null line).

  - `null.line.lty`: line type for horizontal line at 0 (default:
    `null.line.lty=2`.

  - `ci.lty`: line type for confidence interval curves (default:
    `ci.lty=5`).

  - `ci.lwd`: line width for confidence interval curves (default:
    `ci.lwd=1.5`).

  - `ci.col`: line color for confidence interval curves (default:
    `ci.col="darkgrey"`).

- ...:

  Additional graphical arguments passed to function
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html) or
  function [`curve`](https://rdrr.io/r/graphics/curve.html)

## Value

Plot the average treatment effect (ATE) results from a tteICE object

## See also

[`plot.default`](https://rdrr.io/r/graphics/plot.default.html),
[`points`](https://rdrr.io/r/graphics/points.html),
[`curve`](https://rdrr.io/r/graphics/curve.html),
[`plot.tteICE`](https://mephas.github.io/tteICE/reference/plot.tteICE.md)

## Examples

``` r
## Load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
bmt$A = A

## simple model fitting and plotting
library(survival)
fit = tteICE(Surv(t2,d4,type = "mstate")~A, data=bmt)
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38
plot_ate(fit)
#> Error: object 'fit' not found

## model fitting using competing risk data
fit1 = surv.tteICE(A, bmt$t2, bmt$d4, 'composite')
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38

## Plot asymptotic confidence intervals based on explicit formulas
plot_ate(fit1, ylim=c(-0.4,0.4))
#> Error: object 'fit1' not found

## Plot bootstrap confidence intervals
fit2 = surv.tteICE(A, bmt$t2, bmt$d4, 'natural', nboot=50) ## SE=0??
plot_ate(fit2, ylim=c(-0.4,0.4))


## Model with semicompeting risk data
fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38

## Plot asymptotic confidence intervals based on explicit formulas
plot_ate(fit3, ylim=c(-0.4,0.4),
         plot.configs=list(add.null.line=FALSE))
#> Error: object 'fit3' not found

## Plot bootstrap confidence intervals
fit4 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2,
                  "composite", nboot=50)           ## SE=0??
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38

plot_ate(fit4, ylim=c(-0.4,0.4),
         plot.configs=list(add.null.line=FALSE, lty=2, main=""))
#> Error: object 'fit4' not found
```
