# Summary method for 'tteICE' objects

This function summarizes the results

## Usage

``` r
# S3 method for class 'tteICE'
summary(object, ...)
```

## Arguments

- object:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- ...:

  Other arguments in function
  [`summary`](https://rdrr.io/r/base/summary.html)

## Value

A list that consists of summaries of a tteICE object: data type,
strategy, estimation method, maximum follow-up time, sample size,
treated sample size, controlled sample size, p-value, and predicted
risks at quartiles

## See also

[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.md),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.md),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE.md),
[`print.tteICE`](https://mephas.github.io/tteICE/reference/print.tteICE.md)

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
bmt$A = A
X = as.matrix(bmt[,c('z1','z3','z5')])

## Composite variable strategy,
## nonparametric estimation without covariates
fit1 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
summary(fit1)
#> $dtype
#> [1] "smcmprsk"
#> 
#> $strategy
#> [1] "composite"
#> 
#> $method
#> [1] "np"
#> 
#> $maxt
#> [1] 2640
#> 
#> $n
#> [1] 137
#> 
#> $n1
#> [1] 99
#> 
#> $n0
#> [1] 38
#> 
#> $p.val
#> [1] 0.5906978
#> 
#> $est
#>               660        1320        1980         2640
#> CIF1   0.53226259  0.58641246  0.58641246  0.629905604
#> se1    0.05012612  0.04988569  0.04988569  0.061748891
#> CIF0   0.60870186  0.63767315  0.63767315  0.637673151
#> se0    0.08026005  0.07929563  0.07929563  0.079295626
#> ATE   -0.07643926 -0.05126070 -0.05126070 -0.007767547
#> se     0.09462718  0.09368233  0.09368233  0.100502347
#> p.val  0.41920919  0.58425802  0.58425802  0.938395060
#> 
#> attr(,"class")
#> [1] "summary.tteICE"

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

library(survival)
fit3 = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5, 
              data=bmt, strategy="composite", method='eff')
summary(fit3)
#> $dtype
#> [1] "cmprsk"
#> 
#> $strategy
#> [1] "composite"
#> 
#> $method
#> [1] "eff"
#> 
#> $maxt
#> [1] 2640
#> 
#> $n
#> [1] 137
#> 
#> $n1
#> [1] 99
#> 
#> $n0
#> [1] 38
#> 
#> $p.val
#> [1] 0.1365913
#> 
#> $est
#>               660        1320        1980        2640
#> CIF1   0.52459066  0.58362684  0.58362684  0.63470334
#> se1    0.05140949  0.05108944  0.05108944  0.05874485
#> CIF0   0.67836383  0.70103992  0.70103992  0.70103992
#> se0    0.06623568  0.06510731  0.06510731  0.06510731
#> ATE   -0.15377316 -0.11741308 -0.11741308 -0.06633659
#> se     0.08384570  0.08275925  0.08275925  0.08769219
#> p.val  0.06665373  0.15597758  0.15597758  0.44936692
#> 
#> attr(,"class")
#> [1] "summary.tteICE"
```
