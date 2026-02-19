#' @title Fit CIFs using treatment policy strategy for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using treatment policy strategy (competing risks data structure). This strategy ignores the intercurrent
#' event and uses the time to the primary event as it was recorded.
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
#' @seealso \code{\link[tteICE]{surv.treatment.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#' @keywords internal

surv.treatment <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit1 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus==1)[s1], weights=weights[s1])
  fit0 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus==1)[s0], weights=weights[s0])
  tt1 = c(0,fit1$time)
  tt0 = c(0,fit0$time)
  tt = sort(unique(c(tt1,tt0)))
  cumhaz1 = c(0, fit1$cumhaz)
  cumhaz0 = c(0, fit0$cumhaz)
  cif1 = 1 - exp(-cumhaz1)
  cif0 = 1 - exp(-cumhaz0)
  se1 = c(0, fit1$std.err * fit1$surv)
  se0 = c(0, fit0$std.err * fit0$surv)
  se1[is.na(se1)] = rev(na.omit(se1))[1]
  se0[is.na(se0)] = rev(na.omit(se0))[1]
  cumhaz = data.frame(time=tt,cumhaz1=.matchy(cumhaz1,tt1,tt),cumhaz0=.matchy(cumhaz0,tt0,tt))
  surv_diff = survdiff(Surv(Time,cstatus==1)~A)
  p = pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail=FALSE)
  ate = .matchy(cif1,tt1,tt)-.matchy(cif0,tt0,tt)
  se = sqrt(.matchy(se1,tt1,tt)^2+.matchy(se0,tt0,tt)^2)
  return(list(time1=tt1,time0=tt0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,cumhaz=cumhaz))
}
