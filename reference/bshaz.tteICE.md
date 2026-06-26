# Baseline hazards of 'tteICE' objects

This function extracts the baseline cumulative hazards in the survival
models

## Usage

``` r
# S3 method for class 'tteICE'
bshaz(x)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

## Value

A data frame of baseline cumulative hazards in the working Kaplan-Meier
or Cox models, stratified by treatment groups. The first column is time,
the following columns are baseline cumulative hazards.

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
X = as.matrix(bmt[,c('z1','z3','z5')])
bmt$A = A

fit = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
 data=bmt, strategy="whileon", method='eff')
#> Warning: type= 'mstate' is deprecated, use a factor variable as status
#> Error in `[.cox.zph`(cox.zph(fit11, terms = FALSE), , 3): argument "..1" is missing, with no default
bshaz(fit)
#> Error: object 'fit' not found
```
