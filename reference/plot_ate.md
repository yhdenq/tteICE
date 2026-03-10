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
