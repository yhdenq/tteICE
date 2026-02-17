# Fit CIFs using while on treatment strategy for semicompeting risks data

This function nonparametrically estimates the potential cumulative
incidence function using while on treatment strategy (semicompeting
risks data structure). This strategy can be understood as the competing
risks model, which gives the subdistribution of the primary event.

## Usage

``` r
scr.whileon(A, Time, status, Time_int, status_int, weights = rep(1, length(A)))
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

  P value of testing the treatment effect based on Gray test.

## Details

## See also

[`scr.whileon.eff`](https://mephas.github.io/tteICE/reference/scr.whileon.eff.md),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html)
