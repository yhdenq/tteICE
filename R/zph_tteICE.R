#' @title Checking proportional hazards of 'tteICE' objects
#'
#' @description This function checks the proportional hazards assumption in the Cox models using Schoenfeld residuals
#'
#' @param object 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @examples
#' ## load data
#' data(bmt)
#' bmt = transform(bmt, d4=d2+d3)
#' A = as.numeric(bmt$group>1)
#' X = as.matrix(bmt[,c('z1','z3','z5')])
#' bmt$A = A
#' library(survival)
#' fit = tteICE(Surv(t2, d4, type = "mstate")~A|z1+z3+z5,
#'  data=bmt, strategy="whileon", method='eff')
#' print(fit$ph)
#' plot(fit$ph$ph11)
#' plot(fit$ph$ph10)
#'
#' @return
#' A list of P-values of testing the proportional hazards (PH) assumption in the working Cox models, for each 
#' covariate and a global test, stratified by treatment groups. 
#' For the treatment policy strategy and composite variable strategy, only one Cox model is fit (for the primary 
#' outcome event or the composite event). In these two strategies, \code{ph1} is the P-values in the treated 
#' group, \code{ph0} is the P-values in the control group. For other strategies, Cox models are fitted for 
#' each event (primary outcome event and intercurrent event). In these strategies, \code{ph11} is the P-values 
#' for the primary outcome event in the treatment group, \code{ph10} is the P-values for the primary outcome 
#' event in the control group, \code{ph21} is the P-values for the intercurrent event in the treated group, 
#' \code{ph20} is the P-values for the intercurrent in the control group.
#'
#' @export

coef.tteICE <- function(object) {
  object$ph
}
