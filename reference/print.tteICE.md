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
#> P-value of the average treatment effect: 0.591 

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
```
