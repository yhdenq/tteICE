# Plot estimated cumulative incidence functions (CIFs)

This function plots the estimated potential cumulative incidence
function, along with pointwise confidence intervals.

## Usage

``` r
plot_inc(
  fit,
  decrease = FALSE,
  conf.int = 0.95,
  xlab = "Time",
  xlim = NULL,
  ylim = c(0, 1),
  plot.configs = list(ylab = NULL, main = NULL, lty = 1, lwd = 2, ci.lty = 5, ci.lwd =
    1.5, legend = c("Treated", "Control"), col = c("brown", "darkcyan"), legend.cex =
    0.9, show.p.value = TRUE),
  ...
)
```

## Arguments

- fit:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- decrease:

  A logical variable indicating the type of curve to display. If
  `decrease = FALSE` (default), cumulative incidence functions (CIFs)
  are plotted. If `decrease = TRUE`, survival functions are plotted
  instead.

- conf.int:

  Confidence level for the pointwise confidence intervals If
  `conf.int = NULL`, no confidence intervals are provided.

- xlab:

  Label for the x-axis.

- xlim:

  A numeric vector of length 2 specifying the limits of the x-axis. If
  `xlim = NULL` (default), the limits are determined automatically from
  the data.

- ylim:

  A numeric vector of length 2 specifying the limits of the y-axis.
  Defaults to `ylim = c(0, 1)`.

- plot.configs:

  A named `list` of additional plot configurations. Common entries
  include:

  - `ylab`: character, label for the y-axis (default: `ylab=NULL`, use
    the default label).

  - `main`: character, title for the plot (default: `main=NULL`, use the
    default label).

  - `lty`: line type for the curve (default: `lty=1`).

  - `lwd`: line width for the curve (default: `lwd=2`).

  - `ci.lty`: line type for confidence interval curves (default:
    `ci.lty=5, ci.lwd=1.5`).

  - `ci.lwd`: line width for confidence interval curves (default:
    `ci.lwd=1.5`).

  - `legend`: legend for the two group (default:
    `legend=c('Treated','Control')`).

  - `col`: color of the curve for the two group (default:
    `col=c('brown','darkcyan')`).

  - `legend.cex`: font size for the legend (default: `legend.cex=0.9`).

  - `show.p.value`: whether to show the p-value between two groups
    (default: `show.p.value=TRUE`, show the p-value)

- ...:

  Additional graphical arguments passed to function
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html) or
  function [`curve`](https://rdrr.io/r/graphics/curve.html)

## Value

Plot the cumulative incidence function results from a tteICE object

## See also

[`plot.default`](https://rdrr.io/r/graphics/plot.default.html),
[`points`](https://rdrr.io/r/graphics/points.html),
[`curve`](https://rdrr.io/r/graphics/curve.html),
[`plot.tteICE`](https://mephas.github.io/tteICE/reference/plot.tteICE.md)

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
bmt$A = A

## simple model fitting and plotting
library(survival)
fit = tteICE(Surv(t2,d4,type = "mstate")~A, data=bmt)
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38
plot_inc(fit)
#> Error: object 'fit' not found

## model fitting using competing risk data
fit1 = surv.tteICE(A, bmt$t2, bmt$d4, 'treatment')

## plot asymptotic confidence intervals based on explicit formulas
plot_inc(fit1, ylim=c(0,1),
         plot.configs=list(legend=c('AML','ALL'), show.p.value=FALSE) )


## plot bootstrap confidence intervals
fit2 = surv.tteICE(A, bmt$t2, bmt$d4, 'treatment', nboot=50)
plot_inc(fit2, ylim=c(0,1),
         plot.configs=list(legend=c('AML','ALL')))


## model fitting using semicompeting risk data
fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38

## plot asymptotic confidence intervals based on explicit formulas
plot_inc(fit3, ylim=c(0,1), plot.configs=list(add.null.line=FALSE))
#> Error: object 'fit3' not found

## plot bootstrap confidence intervals
fit4 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2,
                  "composite", nboot=50) ##??
#> Error in data.frame(time = tt, cumhaz1 = cumhaz1, cumhaz0 = cumhaz0): arguments imply differing number of rows: 131, 95, 38
plot_inc(fit4, ylim=c(0,1),
         plot.configs=list(lty=2, lwd=3, main="My title"))
#> Error: object 'fit4' not found
```
