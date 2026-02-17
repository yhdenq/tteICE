# Fit CIFs using hypothetical strategy (II) for competing risks data

This function nonparametrically estimates the potential cumulative
incidence function using hypothetical strategy (competing risks data
structure). The intercurrent event is assumed to be absent in the
hypothetical scenario.

## Usage

``` r
surv.removed(A, Time, cstatus, weights = rep(1, length(A)))
```

## Arguments

- A:

  Treatment indicator, 1 for treatment and 0 for control.

- Time:

  Time to event.

- cstatus:

  Indicator of event, 1 for the primary event, 2 for the intercurrent
  event, 0 for censoring.

- weights:

  Weight for each subject.

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

  P value of testing the treatment effect based on logrank test.

## Details

## See also

[`surv.removed.eff`](https://mephas.github.io/tteICE/reference/surv.removed.eff.md),
[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.md)
