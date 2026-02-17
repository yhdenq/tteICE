# Fit CIFs for semicompeting risks time-to-event data with intercurrent events.

This function estimates the potential cumulative incidence function for
time-to event data under ICH E9 (R1) to address intercurrent events. The
input data should be of a semicompeting risks structure.

## Usage

``` r
scr.tteICE(
  A,
  Time,
  status,
  Time_int,
  status_int,
  strategy = "composite",
  cov1 = NULL,
  method = "np",
  weights = NULL,
  subset = NULL,
  na.rm = FALSE,
  nboot = 0,
  seed = 0
)
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

- strategy:

  Strategy to address intercurrent events, `"treatment"` indicating
  treatment policy strategy, `"composite"` indicating composite variable
  strategy, `"natural"` indicating hypothetical strategy (Scenario I,
  controlling the hazard of intercurrent events), `"removed"` indicating
  hypothetical strategy (Scenario II, removing intercurrent events),
  `"whileon"` indicating while on treatment strategy, and `"principal"`
  indicating principal stratum strategy.

- cov1:

  Baseline covariates.

- method:

  Estimation method, `"np"` indicating nonparametric estimation, `"ipw"`
  indicating invserse treatment probability weighting, `"eff"`
  indicating semiparametrically efficient estimation based on efficient
  influence functions.

- weights:

  Weight for each subject.

- subset:

  Subset, either numerical or logical.

- na.rm:

  Whether to remove missing values.

- nboot:

  Number of resamplings in the boostrapping method. If `nboot` is 0 or
  1, then asymptotic standard error based on the explicit form is
  calculated instead of bootstrapping.

- seed:

  Seed for bootstrapping.

## Value

A list including the fitted object and input variables.

## Details

- Background:

  Intercurrent events refer to the events occurring after treatment
  initiation of clinical trials that affect either the interpretation of
  or the existence of the measurements associated with the clinical
  question of interest. The International Conference on Harmonization
  (ICH) E9 (R1) addendum proposed five strategies to address
  intercurrent events, namely, treatment policy strategy, composite
  variable strategy, while on treatment strategy, hypothetical strategy,
  and principal stratum strategy. To answer a specific scientific
  question, a strategy with a particular estimand is chosen before the
  study design.

- Model:

  We adopt the potential outcomes framework that defines a causal
  estimand as the contrast between functionals of potential outcomes.
  Consider a randomized controlled trial with \\n\\ individuals randomly
  assigned to one of two treatment conditions, denoted by \\w\\, where
  \\w = 1\\ represents the active treatment (a test drug) and \\w = 0\\
  represents the control (placebo). Assume that all patients adhere to
  their treatment assignments and do not discontinue treatment.
  Associated with individual \\i = 1, ..., n\\ are two potential
  time-to-event primary outcomes \\T_i(1)\\ and \\T_i(0)\\, if any,
  which represent the time durations from treatment initiation to the
  primary outcome event under two treatment assignments respectively.
  Let \\R_i(1)\\ and \\R_i(0)\\ denote the occurrence time of potential
  intercurrent events, if any, under the two treatment assignments,
  respectively. Intercurrent events are considered as absent if no
  post-treatment intercurrent events occur until the end of study.

- Estimand:

  We adopt the potential cumulative incidences under both treatment
  assignments as the target estimands. Potential cumulative incidences
  describe the probability of time-to-event outcomes occurring at each
  time point. We define the treatment effect as the contrast of two
  potential cumulative incidences. Cumulative incidences are model-free
  and collapsible, enjoying causal interpretations.

## References

Deng, Y., Han, S., & Zhou, X. H. (2025). Inference for Cumulative
Incidences and Treatment Effects in Randomized Controlled Trials With
Time-to-Event Outcomes Under ICH E9 (R1). *Statistics in Medicine*.
[doi:10.1002/sim.70091](https://doi.org/10.1002/sim.70091)

## See also

[`surv.boot`](https://mephas.github.io/tteICE/reference/surv.boot.md),
[`surv.tteICE`](https://mephas.github.io/tteICE/reference/surv.tteICE.html)

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
X = as.matrix(bmt[,c('z1','z3','z5')])

## Composite variable strategy,
## nonparametric estimation without covariates
fit1 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "composite")
# \donttest{
fit10 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "aa") ## warning message
#> Warning: Please choose a strategy from the following:
#>  treatment, composite, natural, removed, whileon, principal
#> 
#>             composite variable strategy is used by default
# }

## Hypothetical strategy (natural effects),
## nonparametric estimation with inverse probability weighting
fit2 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "natural", X, method='ipw')

## nonparametric estimation with weights as non-standardized inverse probability score
ps = predict(glm(A ~ X, family='binomial'), type='response')
w = A/ps + (1-A)/(1-ps)
fit2 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "natural", weights=w)

## Hypothetical strategy (removing intercurrent events),
## semiparametrically efficient estimation with covariates
fit3 = scr.tteICE(A, bmt$t1, bmt$d1, bmt$t2, bmt$d2, "removed", X, method='eff')
```
