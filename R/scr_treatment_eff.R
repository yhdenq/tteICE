#' @title Fit CIFs using treatment policy strategy for semicompeting risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function based on efficient
#' influence functions using treatment policy strategy (semicompeting risks data structure). Cox models are
#' employed for the survival model. This strategy ignores the intercurrent event and uses the time to
#' the primary event as it was recorded.
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
#' \item{p.val}{P value of testing the treatment effect based on the efficient influence function
#' of the restricted mean survival time lost by the end of study.}
#' \item{coef}{Coefficients of covariates in the working Cox models for each event.}
#' \item{ph}{P values of the proportional hazards assumption in the working Cox models for each event.}
#' }
#'
#' @details
#' \describe{
#' The treatment policy strategy addresses the problem of intercurrent events by expanding
#' the initial treatment conditions to a treatment policy. This strategy is applicable
#' only if intercurrent events do not hinder primary outcome events. The treatments under
#' comparison are now two treatment policies: \eqn{(w, R(w))}, where \eqn{w = 1, 0}. One policy
#' \eqn{(1,R(1))} involves administering the test drug, along with any naturally occurring
#' intercurrents, whereas the other policy \eqn{(0,R(0))} involves administering a placebo,
#' along with any naturally occurring intercurrents. Thus, the potential outcomes are
#' \eqn{T(1,R(1))} and \eqn{T(0,R(0))}. Instead of comparing the test drug and placebo themselves,
#' the contrast of interest is made between the two treatment policies. The difference in
#' cumulative incidences under the two treatment policies is then
#' \eqn{\tau(t) = P(T(1, R(1)) < t) - P(T(0, R(0)) < t),}{ATE_tp}
#' representing the difference in probabilities of experiencing primary outcome events during
#' \eqn{(0,t)} under active treatment and placebo. \cr
#' The average treatment effect \eqn{\tau^{\text{tp}}(t)} has a meaningful causal interpretation
#' only when \eqn{T(1, R(1))} and \eqn{T(0, R(0))} are well defined. Because the treatment policy
#' includes the occurrence of the intercurrent event as natural, the entire treatment policy is
#' determined by manipulating the initial treatment condition $w$ only. Therefore, we can simplify
#' the notations \eqn{T(w, R(w)) = T(w)} in defining estimands. As such,
#' \eqn{\tau(t) = P(T(1)) < t) - P(T(0) < t)} as the intention-to-treat analysis.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.treatment}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.treatment.eff <- function(A,Time,status,Time_int,status_int,X=NULL){
  fit = surv.treatment.eff(A,Time,status,X)
  return(fit)
}
