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

[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html),
[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.html),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE-package.html)
[`surv.boot`](https://mephas.github.io/tteICE/reference/surv.boot.md)

[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.html),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE-package.html)

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
predict(fit3, timeset=c(670,2000))
#>               670        2000
#> CIF1   0.52459066  0.58362684
#> se1    0.05140949  0.05108944
#> CIF0   0.70103992  0.70103992
#> se0    0.06510731  0.06510731
#> ATE   -0.17644926 -0.11741308
#> se     0.08295721  0.08275925
#> p.val  0.03342080  0.15597758
predict(fit3)
#>               660        1320        1980        2640
#> CIF1   0.52459066  0.58362684  0.58362684  0.63470334
#> se1    0.05140949  0.05108944  0.05108944  0.05874485
#> CIF0   0.67836383  0.70103992  0.70103992  0.70103992
#> se0    0.06623568  0.06510731  0.06510731  0.06510731
#> ATE   -0.15377316 -0.11741308 -0.11741308 -0.06633659
#> se     0.08384570  0.08275925  0.08275925  0.08769219
#> p.val  0.06665373  0.15597758  0.15597758  0.44936692

```
