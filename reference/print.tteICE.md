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
#> Error in data.frame(time = tt, cumhaz1 = haz1, cumhaz0 = haz0): arguments imply differing number of rows: 131, 95, 38
print(fit1)
#> Error: object 'fit1' not found

fit2 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
#> Error in data.frame(time = tt, cumhaz1 = haz1, cumhaz0 = haz0): arguments imply differing number of rows: 131, 95, 38
print(fit2, digits=2)
#> Error: object 'fit2' not found

library(survival)
fit3 = tteICE(Surv(t2, d4, type = "mstate")~A,
              data=bmt, strategy="composite", method='eff')
#> Error in data.frame(time = tt, cumhaz1 = haz1, cumhaz0 = haz0): arguments imply differing number of rows: 131, 95, 38
print(fit3, digits=3)
#> Error: object 'fit3' not found
```
