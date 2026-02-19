#' @title Fit CIFs using composite variable strategy for semicompeting risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using composite variable strategy (semicompeting risks data structure). This strategy adopts the
#' first occurrence of either the intermediate or primary event as the event of interest.
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
#' The composite variable strategy addresses the problem of intercurrent events by expanding the
#' outcome variables. It aggregates the intercurrent event and the primary outcome event into a single
#' composite outcome variable. The idea is not new in the context of progression-free survival,
#' where the composite outcome variable is defined as the occurrence of either a non-terminal event
#' (e.g., cancer progression) or a terminal event (e.g., death). One widely used composite outcome
#' variable has the form \eqn{Q(w) = \min\{T(w), R(w)\}} for \eqn{w = 1, 0}. When this simple form
#' is adopted, the difference in counterfactual cumulative incidences is
#' \eqn{\tau(t) = P( Q(1) < t ) - P( Q(0) < t ),}
#' representing the difference in probabilities of experiencing either intercurrent events or primary
#' outcome events during \eqn{(0,t)} under active treatment and placebo.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.composite.eff}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.composite <- function(A,Time,status,Time_int,status_int,weights=rep(1,length(A))){
  Time = (Time + Time_int - abs(Time-Time_int))/2
  cstatus = status + 2*status_int
  cstatus[cstatus>2] = 2
  fit = surv.composite(A,Time,cstatus,weights)
  return(fit)
}
