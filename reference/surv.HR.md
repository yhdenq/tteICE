# Estimate hazard ratios

This function estimates the hazard ratio for time-to event data under
ICH E9 (R1) to address intercurrent events. Multiple strategies except
the principal stratum strategy are allowed.

## Usage

``` r
surv.HR(
  A,
  Time,
  cstatus,
  strategy = "composite",
  cov1 = NULL,
  conf.int = 0.95,
  weights = NULL,
  subset = NULL
)
```

## Arguments

- A:

  Treatment indicator, 1 for treatment and 0 for control.

- Time:

  Time to event.

- cstatus:

  Indicator of event, 1 for the primary event, 2 for the intercurrent
  event, 0 for censoring.

- strategy:

  Strategy to address intercurrent events, `"treatment"` indicating
  treatment policy strategy, `"composite"` indicating composite variable
  strategy, `"natural"` indicating hypothetical strategy (Scenario I,
  controlling the hazard of intercurrent events), `"removed"` indicating
  hypothetical strategy (Scenario II, removing intercurrent events), and
  `"whileon"` indicating while on treatment strategy.

- cov1:

  Baseline covariates.

- conf.int:

  Level of the confidence interval.

- weights:

  Weight for each subject (not applied to the while on treatment
  strategy).

- subset:

  Subset, either numerical or logical.

## Value

A list including

- logHR:

  Estimated log hazard ratio (logHR) of the treatment effect on the
  primary event.

- se:

  Standard error of the estimated log hazard ratio (logHR).

- CI:

  Confidence interval of the hazard ratio (HR).

- p.val:

  P value of the hazard ratio.

## Details

## Examples

``` r
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)

## composite variable strategy
fit = surv.HR(A, bmt$t2, bmt$d4, "composite")

## while on treatment strategy
X = bmt[,c('z1','z3','z5')]
fit = surv.HR(A, bmt$t2, bmt$d4, "whileon", cov1=X)

```
