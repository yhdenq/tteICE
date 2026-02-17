# Calculate standard errors for estimated CIFs and treatment effects

This function calculates the standard error for the estimated potential
cumulative incidence function and treatment effect. Two methods to
calculate the standard error are considered: the asymptotic standard
error based on the explicit formula and bootstrapping.

## Usage

``` r
surv.boot(fit, nboot = 0, seed = NULL)
```

## Arguments

- fit:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- nboot:

  Number of resamplings in the boostrapping method. If `nboot` is 0 or
  1, then asymptotic standard error based on the explicit form is
  calculated instead of bootstrapping.

- seed:

  Seed for bootstrapping.

## Value

A list including

- time:

  Time points in both groups.

- cif1:

  Estimated cumulative incidence function in the treated group.

- cif0:

  Estimated cumulative incidence function in the control group.

- se1:

  Standard error of the estimated cumulative incidence function in the
  treated group.

- se0:

  Standard error of the estimated cumulative incidence function in the
  control group.

- ate:

  Estimated treatment effect (difference in cumulative incidence
  functions).

- se:

  Standard error of the estimated treatment effect.

- strategy:

  Strategy used.

- method:

  Estimation method used.

## See also

[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.html),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html),
[`tteICE`](https://mephas.github.io/tteICE/reference/tteICE-package.html)
