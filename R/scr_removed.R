#' @title Fit CIFs using hypothetical strategy (II) for semicompeting risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using hypothetical strategy (semicompeting risks data structure). The intercurrent event is assumed to
#' be absent in the hypothetical scenario.
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
#' @param weights Weight for each subject.
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
#' \item{p.val}{P value of testing the treatment effect based on logrank test.}
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
#' event \eqn{\lambda_2(t; w)} while assuming the hazard specific to the primary outcome event
#' \eqn{\lambda_1(t; w)} remains unchanged. Specifically, we envision that intercurrent events are
#' absent in the hypothetical scenario for all individuals, so \eqn{\lambda_2'(t;0) = \Lambda_2'(t;1) = 0}.
#' This hypothetical scenario leads to an estimand called the marginal cumulative incidence. The treatment
#' effect corresponds to the controlled direct effect with the intercurrent events removed.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.removed.eff}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.removed <- function(A,Time,status,Time_int,status_int,weights=rep(1,length(A))){
  Time = (Time + Time_int - abs(Time-Time_int))/2
  cstatus = status + 2*status_int
  cstatus[cstatus>2] = 2
  fit = surv.removed(A,Time,cstatus,weights)
  return(fit)
}
