#' @title Baseline hazards of 'tteICE' objects
#'
#' @description This function extracts the baseline cumulative hazards in the survival models
#'
#' @param x
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @return
#' A data frame of baseline cumulative hazards in the working Kaplan-Meier or Cox models, stratified by treatment groups.
#' The first column is time, the following columns are baseline cumulative hazards.
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#' bmt$A = A
#'
#' fit = tteICE(Surv(t2, factor(d4))~A|z1+z3+z5,
#'  data=bmt, strategy="whileon", method='eff')
#' bshaz(fit)
#'
#' @method bshaz tteICE
#' @export
#'
bshaz.tteICE <- function(x) {
  x$cumhaz
}
