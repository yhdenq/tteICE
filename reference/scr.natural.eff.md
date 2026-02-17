# Fit CIFs using hypothetical strategy (I) for semicompeting risks data, based on efficient influence functions

This function estimates the potential cumulative incidence function
based on efficient influence functions using hypothetical strategy
(semicompeting risks data structure). Cox models are employed for
survival models. The intercurrent event is only permitted under treated
if it would occur under control.

## Usage

``` r
scr.natural.eff(A, Time, status, Time_int, status_int, X = NULL)
```

## Arguments

- A:

  Treatment indicator, 1 for treatment and 0 for control.

- Time:

  Time to the primary (terminal) event.

- status:

  Indicator of the primary (terminal) event, 1 for event and 0 for
  censoring.

- Time_int:

  Time to the intercurrent event.

- status_int:

  Indicator of the intercurrent event, 1 for event and 0 for censoring.

- X:

  Baseline covariates.

## Value

A list including

- time1:

  Time points in the treated group.

- time0:

  Time points in the control group.

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

- time:

  Time points in both groups.

- ate:

  Estimated treatment effect (difference in cumulative incidence
  functions).

- se:

  Standard error of the estimated treatment effect.

- p.val:

  P value of testing the treatment effect based on the efficient
  influence function of the restricted mean survival time lost by the
  end of study.

## Details

## See also

[`scr.natural`](https://mephas.github.io/tteICE/reference/scr.natural.md),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.md)
