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
print(fit$ph)
#> $ph11
#>        chisq df    p
#> Xz1    0.166  1 0.68
#> Xz3    0.767  1 0.38
#> Xz5    0.310  1 0.58
#> GLOBAL 1.131  3 0.77
#> 
#> $ph10
#>        chisq df     p
#> Xz1     5.38  1 0.020
#> Xz3     0.11  1 0.740
#> Xz5     3.19  1 0.074
#> GLOBAL  6.77  3 0.079
#> 
#> $ph21
#>         chisq df    p
#> Xz1    0.0178  1 0.89
#> Xz3    0.0956  1 0.76
#> Xz5    0.0905  1 0.76
#> GLOBAL 0.2487  3 0.97
#> 
#> $ph20
#>         chisq df    p
#> Xz1    0.9115  1 0.34
#> Xz3    0.0551  1 0.81
#> Xz5    2.5960  1 0.11
#> GLOBAL 5.5810  3 0.13
#> 
zph(fit)
#> $ph11
#>        chisq df    p
#> Xz1    0.166  1 0.68
#> Xz3    0.767  1 0.38
#> Xz5    0.310  1 0.58
#> GLOBAL 1.131  3 0.77
#> 
#> $ph10
#>        chisq df     p
#> Xz1     5.38  1 0.020
#> Xz3     0.11  1 0.740
#> Xz5     3.19  1 0.074
#> GLOBAL  6.77  3 0.079
#> 
#> $ph21
#>         chisq df    p
#> Xz1    0.0178  1 0.89
#> Xz3    0.0956  1 0.76
#> Xz5    0.0905  1 0.76
#> GLOBAL 0.2487  3 0.97
#> 
#> $ph20
#>         chisq df    p
#> Xz1    0.9115  1 0.34
#> Xz3    0.0551  1 0.81
#> Xz5    2.5960  1 0.11
#> GLOBAL 5.5810  3 0.13
#> 

plot(fit$ph$ph11)



plot(fit$ph$ph10)





## No results when method is nonparametric
fit.np = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
 data=bmt, strategy="whileon", method='np')
print(fit.np$ph)
#> NULL

```
