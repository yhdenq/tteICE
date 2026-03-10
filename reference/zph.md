# S3 method to checking proportional hazards

This function checks the proportional hazards assumption in the Cox
model using Schoenfeld residuals. This function only return results for
strategies based on efficient influence functions.

## Usage

``` r
zph(x)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

## Value

A list of P-values of testing the proportional hazards (PH) assumption
in the working Cox models, for each covariate and a global test,
stratified by treatment groups.
