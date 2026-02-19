#' @title Fit CIFs using hypothetical strategy (II) for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using hypothetical strategy (competing risks data structure). The intercurrent event is assumed to
#' be absent in the hypothetical scenario.
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
#' \item{cumhaz}{Baseline cumulative hazards in the survival models.}
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
#' @seealso \code{\link[tteICE]{surv.removed.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.removed <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit1 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus==1)[s1], weights=weights[s1])
  fit0 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus==1)[s0], weights=weights[s0])
  tt1 = c(0, fit1$time)
  tt0 = c(0, fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  cumhaz1 = c(0, fit1$cumhaz)
  cumhaz0 = c(0, fit0$cumhaz)
  cumhaz = data.frame(time=tt,cumhaz1=cumhaz1,cumhaz0=cumhaz0)
  cif1 = 1 - exp(-cumhaz1)
  cif0 = 1 - exp(-cumhaz0)
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  surv_diff = survdiff(Surv(Time,cstatus==1)~A)
  p = pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail=FALSE)
  ate = .matchy(cif1,tt1,tt)-.matchy(cif0,tt0,tt)
  se = sqrt(.matchy(se1,tt1,tt)^2+.matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,cumhaz=cumhaz))
}
