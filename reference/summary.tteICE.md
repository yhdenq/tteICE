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
#> Input:
#> scr.tteICE(A = A, Time = bmt$t1, status = bmt$d1, Time_int = bmt$t2, 
#>     status_int = bmt$d2, strategy = "composite")
#> -----------------------------------------------------------------------
#> Data type: smcmprsk 
#> Strategy: composite 
#> Estimation method: np 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.591 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>          660   1320   1980   2640
#> CIF1   0.532  0.586  0.586  0.630
#> se1    0.050  0.050  0.050  0.062
#> CIF0   0.609  0.638  0.638  0.638
#> se0    0.080  0.079  0.079  0.079
#> ATE   -0.076 -0.051 -0.051 -0.008
#> se     0.095  0.094  0.094  0.101
#> p.val  0.419  0.584  0.584  0.938
#> 

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
#> Input:
#> tteICE(formula = Surv(t2, d4, type = "mstate") ~ A | z1 + z3 + 
#>     z5, data = bmt, strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: cmprsk 
#> Strategy: composite 
#> Estimation method: eff 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.137 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>          660   1320   1980   2640
#> CIF1   0.525  0.584  0.584  0.635
#> se1    0.051  0.051  0.051  0.059
#> CIF0   0.678  0.701  0.701  0.701
#> se0    0.066  0.065  0.065  0.065
#> ATE   -0.154 -0.117 -0.117 -0.066
#> se     0.084  0.083  0.083  0.088
#> p.val  0.067  0.156  0.156  0.449
#> 
```
