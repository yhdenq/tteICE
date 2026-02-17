# Print method for 'tteICE' objects

This function summarizes the results

## Usage

``` r
# S3 method for class 'tteICE'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- digits:

  The digits of the results

- ...:

  Other arguments in function
  [`print.default`](https://rdrr.io/r/base/print.default.html)

## Value

Print the summary of a tteICE object

## See also

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

## print the results
fit1 = surv.tteICE(A, bmt$t2, bmt$d4, "composite")
print(fit1)
#> Input:
#> surv.tteICE(A = A, Time = bmt$t2, cstatus = bmt$d4, strategy = "composite")
#> -----------------------------------------------------------------------
#> Data type: competing risks 
#> Strategy: composite variable strategy 
#> Estimation method: nonparametric estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.5907 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5323  0.5864  0.5864  0.6299
#> se1    0.0501  0.0499  0.0499  0.0617
#> CIF0   0.6087  0.6377  0.6377  0.6377
#> se0    0.0803  0.0793  0.0793  0.0793
#> ATE   -0.0764 -0.0513 -0.0513 -0.0078
#> se     0.0946  0.0937  0.0937  0.1005
#> p.val  0.4192  0.5843  0.5843  0.9384
#> 

fit2 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
print(fit2, digits=2)
#> Input:
#> scr.tteICE(A = A, Time = bmt$t1, status = bmt$d1, Time_int = bmt$t2, 
#>     status_int = bmt$d2, strategy = "composite")
#> -----------------------------------------------------------------------
#> Data type: semicompeting risks 
#> Strategy: composite variable strategy 
#> Estimation method: nonparametric estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.59 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>         660  1320  1980  2640
#> CIF1   0.53  0.59  0.59  0.63
#> se1    0.05  0.05  0.05  0.06
#> CIF0   0.61  0.64  0.64  0.64
#> se0    0.08  0.08  0.08  0.08
#> ATE   -0.08 -0.05 -0.05 -0.01
#> se     0.09  0.09  0.09  0.10
#> p.val  0.42  0.58  0.58  0.94
#> 

library(survival)
fit3 = tteICE(Surv(t2, d4, type = "mstate")~A,
              data=bmt, strategy="composite", method='eff')
print(fit3, digits=3)
#> Input:
#> tteICE(formula = Surv(t2, d4, type = "mstate") ~ A, data = bmt, 
#>     strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: competing risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
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
```
