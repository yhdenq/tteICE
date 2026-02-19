#' @title Fit CIFs using hypothetical strategy (II) for semicompeting risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function based on
#' efficient influence functions using hypothetical strategy (semicompeting risks data structure).
#' Cox models are employed for survival models. The intercurrent event is assumed to be absent
#' in the hypothetical scenario.
#'
#' @param A Treatment indicator, 1 for treatment and 0 for control.
#'
#' @param Time Time to the primary (terminal) event.
#'
#' @param status Indicator of the primary (terminal) event, 1 for event and 0 for censoring.
#'
#' @param Time_int Time to the intercurrent event.
#'
#' @param status_int Indicator of the intercurrent event, 1 for event and 0 for censoring.
#'
#' @param X Baseline covariates.
#'
#'
#' @return A list including
#' \describe{
#' \item{time1}{Time points in the treated group.}
#' \item{time0}{Time points in the control group.}
#' \item{cif1}{Estimated cumulative incidence function in the treated group.}
#' \item{cif0}{Estimated cumulative incidence function in the control group.}
#' \item{se1}{Standard error of the estimated cumulative incidence function in the treated group.}
#' \item{se0}{Standard error of the estimated cumulative incidence function in the control group.}
#' \item{time}{Time points in both groups.}
#' \item{ate}{Estimated treatment effect (difference in cumulative incidence functions).}
#' \item{se}{Standard error of the estimated treatment effect.}
#' \item{p.val}{P value of testing the treatment effect based on the efficient influence function of
#' the restricted mean survival time lost by the end of study.}
#' \item{coef}{Coefficients of covariates in the working Cox models for each event.}
#' \item{ph}{P values of the proportional hazards assumption in the working Cox models for each event.}
#' }
#'
#' @details
#' \describe{
#' The hypothetical strategy envisions a hypothetical clinical trial condition where the occurrence
#' of intercurrent events is restricted in certain ways. By doing so, the distribution of potential
#' outcomes under the hypothetical scenario can capture the impact of intercurrent events explicitly
#' through a pre-specified criterion. We use \eqn{T'(w)}, \eqn{w = 1, 0} to denote the time to the
#' primary outcome event in the hypothetical scenario. The time-dependent treatment effect specific
#' to this hypothetical scenario is written as
#' \eqn{\tau(t) = P(T'(1) < t) - P(T'(0) < t),}
#' representing the difference in probabilities of experiencing primary outcome events during \eqn{(0,t)}
#' in the pre-specified hypothetical scenario under active treatment and placebo. \cr
#' The key question is how to envision \eqn{T'(w)}. We manipulate the hazard specific to intercurrent
#' event \eqn{\lambda_2(t; w)} while assuming the cause-specific hazard specific to the primary outcome event
#' under no intercurrent events \eqn{\lambda_1(t; w)} remains unchanged. Specifically, we envision that
#' intercurrent events are absent in the hypothetical scenario for all individuals, so
#' \eqn{\lambda_2'(t;0) = \Lambda_2'(t;1) = 0}. This hypothetical scenario leads to an estimand called the
#' marginal cumulative incidence. The treatment effect corresponds to the controlled direct effect with the
#' intercurrent events removed.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.removed}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.removed.eff <- function(A,Time,status,Time_int,status_int,X=NULL){
  Time = (Time + Time_int - abs(Time-Time_int))/2
  cstatus = status + 2*status_int
  cstatus[cstatus>2] = 2
  fit = surv.removed.eff(A,Time,cstatus,X)
  return(fit)
}
