# S3 method of baseline hazards

This function extracts the baseline cumulative hazards in the survival
models

## Usage

``` r
bshaz(x)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

## Value

A data frame of baseline cumulative hazards in the working Kaplan-Meier
or Cox models, stratified by treatment groups from the model object.
