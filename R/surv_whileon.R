#' @title Fit CIFs using while on treatment strategy for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using while on treatment strategy (competing risks data structure). This strategy can be understood
#' as the competing risks model, which gives the subdistribution of the primary event.
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
#' @seealso \code{\link[tteICE]{surv.whileon.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.whileon <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit11 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus==1)[s1], weights=weights[s1])
  fit10 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus==1)[s0], weights=weights[s0])
  fit21 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus>1)[s1], weights=weights[s1])
  fit20 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus>1)[s0], weights=weights[s0])
  time1 = c(0, fit11$time)
  time0 = c(0, fit10$time)
  fit11 = rbind(0,cbind(fit11$cumhaz,fit11$std.err))
  fit10 = rbind(0,cbind(fit10$cumhaz,fit10$std.err))
  fit21 = rbind(0,cbind(fit21$cumhaz,fit21$std.err))
  fit20 = rbind(0,cbind(fit20$cumhaz,fit20$std.err))
  dcif1 = exp(-fit11[,1]-fit21[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  V1 = fit11[,2]^2
  V0 = fit21[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit21[,1])+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G1+G0)
  V1 = fit10[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  M1[is.infinite(M1)] = 0
  M0[is.infinite(M0)] = 0
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G1+G0)
  tt = sort(unique(c(time1,time0)))
  ate = .matchy(cif1,time1,tt)-.matchy(cif0,time0,tt)
  se = sqrt(.matchy(se1,time1,tt)^2+.matchy(se0,time0,tt)^2)
  p = cuminc(Time,cstatus,group=A)$Tests[1,2]
  return(list(time1=time1,time0=time0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p))
}
