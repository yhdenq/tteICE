# Checking proportional hazards of 'tteICE' objects

This function checks the proportional hazards assumption in the Cox
model using Schoenfeld residuals. This function only return results for
strategies based on efficient influence functions.

## Usage

``` r
# S3 method for class 'tteICE'
zph(x)
```

## Arguments

- x:

  A fitted object returned by the function `tteICE`, `surv.tteICE`, or
  `scr.tteICE`.

## Value

A list of P-values of testing the proportional hazards (PH) assumption
in the working Cox models, for each covariate and a global test,
stratified by treatment groups. For the treatment policy strategy and
composite variable strategy, only one Cox model is fit (for the primary
outcome event or the composite event). In these two strategies, `ph1` is
the P-values in the treated group, `ph0` is the P-values in the control
group. For other strategies, Cox models are fitted for each event
(primary outcome event and intercurrent event). In these strategies,
`ph11` is the P-values for the primary outcome event in the treatment
group, `ph10` is the P-values for the primary outcome event in the
control group, `ph21` is the P-values for the intercurrent event in the
treated group, `ph20` is the P-values for the intercurrent in the
control group.

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
#> Error in as.data.frame.default(x[[i]], optional = TRUE, stringsAsFactors = stringsAsFactors): cannot coerce class ‘"cox.zph"’ to a data.frame
print(fit$ph)
#> Error: object 'fit' not found
zph(fit)
#> Error: object 'fit' not found

plot(fit$ph$ph11)
#> Error: object 'fit' not found
plot(fit$ph$ph10)
#> Error: object 'fit' not found


## No results when method is nonparametric
fit.np = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
 data=bmt, strategy="whileon", method='np')
#> Warning: type= 'mstate' is deprecated, use a factor variable as status
print(fit.np$ph)
#> NULL

```
