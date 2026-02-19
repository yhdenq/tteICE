#' @title Fit CIFs using composite variable strategy for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using composite variable strategy (competing risks data structure).
#'  This strategy adopts the first occurrence of either the intermediate
#' or primary event as the event of interest.
#'
#' @param A Treatment indicator, 1 for treatment and 0 for control.
#'
#' @param Time Time to event.
#'
#' @param cstatus Indicator of event, 1 for the primary event, 2 for the intercurrent event, 0 for censoring.
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
#' @seealso \code{\link[tteICE]{surv.composite.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.composite <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit1 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus>0)[s1], weights=weights[s1])
  fit0 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus>0)[s0], weights=weights[s0])
  cif1 = c(0, 1 - exp(-fit1$cumhaz))
  cif0 = c(0, 1 - exp(-fit0$cumhaz))
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  surv_diff = survdiff(Surv(Time,cstatus>0)~A)
  p = pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail=FALSE)
  tt1 = c(0,fit1$time)
  tt0 = c(0,fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  ate = .matchy(cif1,tt1,tt)-.matchy(cif0,tt0,tt)
  se = sqrt(.matchy(se1,tt1,tt)^2+.matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}
