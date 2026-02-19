#' @title Fit CIFs using while on treatment strategy for semicompeting risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using while on treatment strategy (semicompeting risks data structure). This strategy can be understood
#' as the competing risks model, which gives the subdistribution of the primary event.
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
#' \item{p.val}{P value of testing the treatment effect based on Gray test.}
#' }
#'
#' @details
#' \describe{
#' The while on treatment strategy considers the measure of outcome variables taken only up to
#' the occurrence of intercurrent events. The failures of primary outcome events should not be
#' counted in the cumulative incidences if intercurrent events occurred. The difference in
#' counterfactual cumulative incidences under this strategy is
#' \eqn{\tau(t) = P(T(1) < t, R(1) \geq t) - P(T(0) < t, R(0) \geq t),}
#' representing the difference in probabilities of experiencing primary outcome events without
#' intercurrent events during \eqn{(0,t)} under active treatment and placebo. The cumulative
#' incidence function is also known as the cause-specific cumulative incidence or subdistribution
#' function. \cr
#' The while on treatment strategy is closely related to the competing risks model. However,
#' for causal interpretations, it is worth emphasizing that the hazard of \eqn{R(1)} may differ
#' from that of \eqn{R(0)}, leading to vast difference in the underlying features of individuals
#' who have not experienced the primary outcome event between treatment conditions at any time
#' \eqn{t \in (0,t^*)}, where \eqn{t^*} is the end of study. When the scientific question of
#' interest is the impact of treatment on the primary outcome event, the estimand \eqn{\tau(t)}
#' is hard to interpret if a systematic difference in the risks of intercurrent events between two
#' treatment conditions under comparison is anticipated.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.whileon.eff}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.whileon <- function(A,Time,status,Time_int,status_int,weights=rep(1,length(A))){
  Time = (Time + Time_int - abs(Time-Time_int))/2
  cstatus = status + 2*status_int
  cstatus[cstatus>2] = 2
  fit = surv.whileon(A,Time,cstatus,weights)
  return(fit)
}
