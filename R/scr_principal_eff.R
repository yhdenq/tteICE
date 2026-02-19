#' @title Fit CIFs using principal stratum strategy for semicompeting risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function based on 
#' efficient influence functions using principal stratum strategy (semicompeting risks data 
#' structure). Cox models are employed for survival models. The estimand is defined in a subpopulation 
#' where intercurrent events would never occur regardless of treatment conditions.
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
#' The principal stratum strategy aims to stratify the population into subpopulations based on the joint
#' potential occurrences of intercurrent events under the two treatment assignments \eqn{(R(1), R(0))}.
#' Suppose we are interested in a principal stratum comprised of individuals who would never experience
#' intercurrent events, regardless of which treatment they receive. This principal stratum can be indicated
#' by \eqn{\{R(1)=R(0)=\infty\}}. The treatment effect is now defined within this subpopulation,
#' \eqn{\tau(t) = P(T(1) < t \mid R(1)=R(0)=\infty) - P(T(0) < t \mid R(1)=R(0)=\infty),}
#' representing the difference in probabilities of experiencing primary outcome events during \eqn{(0,t)}
#' under active treatment and placebo in the subpopulation that will not experience intercurrent events
#' regardless of treatment during \eqn{(0,t)}. A principal ignorability assumption is made for identification.
#' If the size of the target principal stratum is small, the results could be highly variable.
#' }
#'
#' @seealso \code{\link[tteICE]{scr.principal}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.principal.eff <- function(A,Time,status,Time_int,status_int,X=NULL){
  Time = (Time + Time_int - abs(Time-Time_int))/2
  cstatus = status + 2*status_int
  cstatus[cstatus>2] = 2
  fit = surv.principal.eff(A,Time,cstatus,X)
  return(fit)
}



