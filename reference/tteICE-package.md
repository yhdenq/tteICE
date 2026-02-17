# tteICE: Treatment Effect Estimation for Time-to-Event Data with Intercurrent Events

This package aims to analyze treatment effects in clinical trials with
time-to-event outcomes is complicated by intercurrent events. This
package implements methods for estimating and inferring the cumulative
incidence functions for time-to-event (TTE) outcomes with intercurrent
events (ICE) under the five strategies outlined in the ICH E9 (R1)
addendum, see Deng (2025) <doi:10.1002/sim.70091>. This package can be
used for analyzing data from both randomized controlled trials and
observational studies. In general, the data involve a primary outcome
event and, potentially, an intercurrent event. Two data structures are
allowed: competing risks, where only the time to the first event is
recorded, and semicompeting risks, where the times to both the primary
outcome event and intercurrent event (or censoring) are recorded. For
estimation methods, nonparametric estimation (which does not use
covariates) and semiparametrically efficient estimation are presented.

## Details

Main functions:

- [`tteICE`](https://mephas.github.io/tteICE/reference/tteICE-package.html)
  Using formula to fit cumulative incidence functions (CIFs) for
  competing/semicompeting risk time-to-event data with intercurrent
  events.

- [`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html)
  Fit CIFs for semicompeting risk time-to-event data with intercurrent
  events.

- [`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.html)
  Fit CIFs for competing risk time-to-event with intercurrent events.

- [`plot.tteICE`](https://mephas.github.io/tteICE/reference/plot.tteICE.md)
  Plot results from 'tteICE' objects.

- [`print.tteICE`](https://mephas.github.io/tteICE/reference/print.tteICE.md)
  Print a short summary of results from 'tteICE' objects

- [`summary.tteICE`](https://mephas.github.io/tteICE/reference/summary.tteICE.md)
  Summarize results from 'tteICE' objects

- [`predict.tteICE`](https://mephas.github.io/tteICE/reference/predict.tteICE.md)
  Predict risks for 'tteICE' objects at specific time points

- [`tteICEShiny`](https://mephas.github.io/tteICE/reference/tteICE-package.html)
  Interactive Shiny app for the 'tteICE' package

Example data:

- [`bmt`](https://mephas.github.io/tteICE/reference/bmt.md) Data from
  Section 1.3 of Klein and Moeschberger (1997)

## See also

Useful links:

- <https://github.com/mephas/tteICE>

- <https://mephas.github.io/tteICE/>

- Report bugs at <https://github.com/mephas/tteICE/issues>

## Author

**Maintainer**: Yi Zhou <yzhou@pku.edu.cn>

Authors:

- Yuhao Deng <dengyuhao@pku.edu.cn>
