# Predict method for 'tteICE' objects at specific time points

This function predicts the potential cumulative incidence function and
treatment effect at specific time points.

## Usage

``` r
# S3 method for class 'tteICE'
predict(object, timeset = NULL, ...)
```

## Arguments

- object:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- timeset:

  Time at which to predict the risk. If `timeset=NULL`, risks will be
  predict at the quartiles of the maximum follow-up time.

- ...:

  Other arguments in function
  [`predict`](https://rdrr.io/r/stats/predict.html)

## Value

A matrix with each row being time points, potential cumulative
incidences (under treated and under control), treatment effects,
standard errors, and P-values.

predict a tteICE object. The meanings of each row are: time points,
potential cumulative incidences (under treated and under control),
treatment effects, standard errors, and P-values.

## See also

[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.md),
[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.md),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE.md)
[`surv.boot`](https://mephas.github.io/tteICE/reference/surv.boot.md)

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
X = as.matrix(bmt[,c('z1','z3','z5')])

## predict results at specified time points
## model fitting using semicompeting risk data
fit1 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
predict(fit1, timeset=c(670,2000))
#>               670        2000
#> CIF1   0.53226259  0.58641246
#> se1    0.05012612  0.04988569
#> CIF0   0.63767315  0.63767315
#> se0    0.07929563  0.07929563
#> ATE   -0.10541056 -0.05126070
#> se     0.09381057  0.09368233
#> p.val  0.26116016  0.58425802

## predict results without specifying any time points
## model fitting using competing risk data
fit2 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
predict(fit2)
#>               660        1320        1980         2640
#> CIF1   0.53226259  0.58641246  0.58641246  0.629905604
#> se1    0.05012612  0.04988569  0.04988569  0.061748891
#> CIF0   0.60870186  0.63767315  0.63767315  0.637673151
#> se0    0.08026005  0.07929563  0.07929563  0.079295626
#> ATE   -0.07643926 -0.05126070 -0.05126070 -0.007767547
#> se     0.09462718  0.09368233  0.09368233  0.100502347
#> p.val  0.41920919  0.58425802  0.58425802  0.938395060

## a simpler way
library(survival)
fit3 = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
              data=bmt, strategy="composite", method='eff')
#> Warning: type= 'mstate' is deprecated, use a factor variable as status
#> Error in `[.cox.zph`(cox.zph(fit1, terms = FALSE), , 3): argument "..1" is missing, with no default
predict(fit3, timeset=c(670,2000))
#> Error: object 'fit3' not found
predict(fit3)
#> Error: object 'fit3' not found

```
