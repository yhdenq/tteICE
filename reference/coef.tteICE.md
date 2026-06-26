# Coefficients of 'tteICE' objects

This function extracts the coefficients in the Cox models

## Usage

``` r
# S3 method for class 'tteICE'
coef(object, ...)
```

## Arguments

- object:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

- ...:

  Other arguments in function
  [`coef.default`](https://rdrr.io/r/stats/coef.html)

## Value

A list of coefficients of covariates in the working Cox models,
stratified by treatment groups. For the treatment policy strategy and
composite variable strategy, only one Cox model is fit (for the primary
outcome event or the composite event). In these two strategies, `coef1`
is the coefficients in the treated group, `coef0` is the coefficients in
the control group. For other strategies, Cox models are fitted for each
event (primary outcome event and intercurrent event). In these
strategies, `coef11` is the coefficients for the primary outcome event
in the treatment group, `coef10` is the coefficients for the primary
outcome event in the control group, `coef21` is the coefficients for the
intercurrent event in the treated group, `coef20` is the coefficients
for the intercurrent in the control group.

## Examples

``` r
## load data
data(bmt)
bmt = transform(bmt, d4=d2+d3)
A = as.numeric(bmt$group>1)
X = as.matrix(bmt[,c('z1','z3','z5')])
bmt$A = A

fit = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
 data=bmt, strategy="whileon", method='eff')
#> Warning: type= 'mstate' is deprecated, use a factor variable as status
coef(fit)
#>    Primary, A=1         SE Primary, A=0         SE    ICE, A=1         SE
#> z1   0.01807550 0.01966552   0.05860219 0.04900616 -0.01067268 0.02039951
#> z3  -0.11852413 0.37564863   0.28192137 0.73926362 -0.34003123 0.36772740
#> z5   0.08175745 0.38396510  -1.33574510 0.81648428  0.37432047 0.38170065
#>       ICE, A=0         SE
#> z1  0.06910952 0.04865528
#> z3 -0.49625608 0.66123273
#> z5  0.72289467 0.63150591
```
