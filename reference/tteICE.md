# Using formula to fit CIFs for time-to-event data with intercurrent events

This function estimates the potential cumulative incidence function for
time-to event data under ICH E9 (R1) to address intercurrent events. The
input data can be competing or semicompeting risks data structure.

## Usage

``` r
tteICE(
  formula,
  add.scr = NULL,
  data,
  strategy = "composite",
  method = "np",
  weights = NULL,
  subset = NULL,
  na.rm = FALSE,
  nboot = 0,
  seed = 0
)
```

## Arguments

- formula:

  An object of class "formula" (or one that can be coerced to that
  class). A symbolic description of the model to be fitted. For example,
  `formula=Surv(time, status, type="mstate")~treatment | baseline.covariate`.
  The details of model specification are given under ‘Details’.

- add.scr:

  Required for semicompeting data. An object of class "Surv" (or one
  that can be coerced to that class). For example,
  `add.scr=~Surv(time.intercurrent, status.intercurrent)`. The details
  of model specification are given under ‘Details’.

- data:

  Data or object coercible by as.data.frame to a data frame, containing
  the variables in the model.

- strategy:

  Strategy to address intercurrent events, `"treatment"` indicating
  treatment policy strategy, `"composite"` indicating composite variable
  strategy, `"natural"` indicating hypothetical strategy (Scenario I,
  controlling the hazard of intercurrent events), `"removed"` indicating
  hypothetical strategy (Scenario II, removing intercurrent events),
  `"whileon"` indicating while on treatment strategy, and `"principal"`
  indicating principal stratum strategy.

- method:

  Estimation method, `"np"` indicating nonparametric estimation, `"np"`
  indicating inverse treatment probability weighting, `"eff"` indicating
  semiparametrically efficient estimation based on efficient influence
  functions.

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

- Formula specifications:

  The formula should be set as the following two ways.

  When data take format of competing risk data, set the first argument
  `formula = Surv(time, status, type="mstate") ~ treatment | covariate1+covariate2`
  or `formula = Surv(time, status)~ A` without any baseline covariates,
  where `status`=0,1,2 (1 for the primary event, 2 for the intercurrent
  event, and 0 for censoring).

  When data take the format of semicompeting risk data, set the first
  argument
  `formula = Surv(time, status) ~ treatment | covariate1+covariate2` or
  `formula = Surv(time, status)~ A` without any baseline covariates,
  where `status`=0,1 (1 for the primary event and 0 for censoring). In
  addition, the second argument
  `add.scr = ~ Surv(time.intercurrent, status.intercurrent)` is
  required.

## References

Deng, Y., Han, S., & Zhou, X. H. (2025). Inference for Cumulative
Incidences and Treatment Effects in Randomized Controlled Trials With
Time-to-Event Outcomes Under ICH E9 (R1). *Statistics in Medicine*.
[doi:10.1002/sim.70091](https://doi.org/10.1002/sim.70091)

## See also

[`surv.boot`](https://mephas.github.io/tteICE/reference/surv.boot.md),
[`scr.tteICE`](https://mephas.github.io/tteICE/reference/scr.tteICE.html)

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
X = as.matrix(bmt[,c('z1','z3','z5')])
bmt$A = A

library(survival)
## Composite variable strategy,
## nonparametric estimation without covariates
## Composite variable strategy,
## nonparametric estimation without covariates

## model fitting for competing risk data without covariates
fit1 = tteICE(Surv(t2, d4, type = "mstate")~A,
 data=bmt, strategy="composite", method='eff')
print(fit1)
#> Input:
#> tteICE(formula = Surv(t2, d4, type = "mstate") ~ A, data = bmt, 
#>     strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: competing risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.5907 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5323  0.5864  0.5864  0.6299
#> se1    0.0501  0.0499  0.0499  0.0617
#> CIF0   0.6087  0.6377  0.6377  0.6377
#> se0    0.0803  0.0793  0.0793  0.0793
#> ATE   -0.0764 -0.0513 -0.0513 -0.0078
#> se     0.0946  0.0937  0.0937  0.1005
#> p.val  0.4192  0.5843  0.5843  0.9384
#> 

## model fitting for competing risk data without covariates 
## with bootstrap confidence intervals
fit.bt1 = tteICE(Surv(t2, d4, type = "mstate")~A,
 data=bmt, strategy="composite", method='eff', nboot=20, seed=2)
print(fit.bt1)
#> Input:
#> tteICE(formula = Surv(t2, d4, type = "mstate") ~ A, data = bmt, 
#>     strategy = "composite", method = "eff", nboot = 20, seed = 2)
#> -----------------------------------------------------------------------
#> Data type: competing risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.5907 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5323  0.5864  0.5864  0.6299
#> se1    0.0511  0.0574  0.0574  0.0586
#> CIF0   0.6087  0.6377  0.6377  0.6377
#> se0    0.0673  0.0716  0.0716  0.0716
#> ATE   -0.0764 -0.0513 -0.0513 -0.0078
#> se     0.0857  0.0916  0.0916  0.0912
#> p.val  0.3723  0.5759  0.5759  0.9321
#> 

## model fitting for competing risk data with covariates
fit2 = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
 data=bmt, strategy="composite", method='eff')
print(fit2)
#> Input:
#> tteICE(formula = Surv(t2, d4, type = "mstate") ~ A | z1 + z3 + 
#>     z5, data = bmt, strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: competing risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.1366 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5246  0.5836  0.5836  0.6347
#> se1    0.0514  0.0511  0.0511  0.0587
#> CIF0   0.6784  0.7010  0.7010  0.7010
#> se0    0.0662  0.0651  0.0651  0.0651
#> ATE   -0.1538 -0.1174 -0.1174 -0.0663
#> se     0.0838  0.0828  0.0828  0.0877
#> p.val  0.0667  0.1560  0.1560  0.4494
#> 

## model fitting for semicompeting risk data without covariates
fitscr1 = tteICE(Surv(t1, d1)~A, ~Surv(t2, d2),
 data=bmt, strategy="composite", method='eff')
print(fitscr1)
#> Input:
#> tteICE(formula = Surv(t1, d1) ~ A, add.scr = ~Surv(t2, d2), data = bmt, 
#>     strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: semicompeting risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.5907 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5323  0.5864  0.5864  0.6299
#> se1    0.0501  0.0499  0.0499  0.0617
#> CIF0   0.6087  0.6377  0.6377  0.6377
#> se0    0.0803  0.0793  0.0793  0.0793
#> ATE   -0.0764 -0.0513 -0.0513 -0.0078
#> se     0.0946  0.0937  0.0937  0.1005
#> p.val  0.4192  0.5843  0.5843  0.9384
#> 

## model fitting for semicompeting risk data without covariates
fitscr2 = tteICE(Surv(t1, d1)~A|z1+z3+z5, ~Surv(t2, d2),
 data=bmt, strategy="composite", method='eff')
print(fitscr2)
#> Input:
#> tteICE(formula = Surv(t1, d1) ~ A | z1 + z3 + z5, add.scr = ~Surv(t2, 
#>     d2), data = bmt, strategy = "composite", method = "eff")
#> -----------------------------------------------------------------------
#> Data type: semicompeting risks 
#> Strategy: composite variable strategy 
#> Estimation method: semiparametrically efficient estimation 
#> Observations: 137 (including 99 treated and 38 control)
#> Maximum follow-up time: 2640 
#> P-value of the average treatment effect: 0.1366 
#> -----------------------------------------------------------------------
#> The estimated cumulative incidences and treatment effects at quartiles:
#>           660    1320    1980    2640
#> CIF1   0.5246  0.5836  0.5836  0.6347
#> se1    0.0514  0.0511  0.0511  0.0587
#> CIF0   0.6784  0.7010  0.7010  0.7010
#> se0    0.0662  0.0651  0.0651  0.0651
#> ATE   -0.1538 -0.1174 -0.1174 -0.0663
#> se     0.0838  0.0828  0.0828  0.0877
#> p.val  0.0667  0.1560  0.1560  0.4494
#> 

```
