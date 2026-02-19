#' @title Fit CIFs using principal stratum strategy for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using principal stratum strategy (competing risks data structure). The estimand is defined in a
#' subpopulation where intercurrent events would never occur regardless of treatment conditions.
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
#' \item{p.val}{P value of testing the treatment effect, which is not available under this strategy.}
#' \item{cumhaz}{Baseline cumulative hazards in the survival models.}
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
#' @seealso \code{\link[tteICE]{surv.principal.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.principal <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit11 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus==1)[s1], weights=weights[s1])
  fit10 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus==1)[s0], weights=weights[s0])
  fit21 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus>1)[s1], weights=weights[s1])
  fit20 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus>1)[s0], weights=weights[s0])
  time1 = c(0, fit11$time)
  time0 = c(0, fit10$time)
  tt = sort(unique(c(time1,time0)))
  fit11 = rbind(0,cbind(fit11$cumhaz,fit11$std.err))
  fit10 = rbind(0,cbind(fit10$cumhaz,fit10$std.err))
  fit21 = rbind(0,cbind(fit21$cumhaz,fit21$std.err))
  fit20 = rbind(0,cbind(fit20$cumhaz,fit20$std.err))
  cumhaz = data.frame(time=tt,cumhaz11=.matchy(fit11[,1],time1,tt),cumhaz10=.matchy(fit10[,1],time0,tt),
                      cumhaz21=.matchy(fit21[,1],time1,tt),cumhaz20=.matchy(fit20[,1],time0,tt))
  S1 = exp(-fit11[,1]-fit21[,1])
  S0 = exp(-fit10[,1]-fit20[,1])
  dcif1 = S1*diff(c(0,fit11[,1]))
  dcif0 = S0*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  PR1 = min(S1 + cif1)
  PR0 = min(S0 + cif0)
  V1 = fit11[,2]^2
  V0 = fit21[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif1^2 + cumsum(M1*(S1+cif1)^2) -
    2*cif1*cumsum(M1*(S1+cif1))
  G0 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  G3 = cif1^2/PR1^2*sum(M1*(S1+cif1-PR1)^2)
  G2 = cif1^2/PR1^2*sum(M0*(cif1-PR1)^2)
  G5 = 2*cif1/PR1*(cumsum(M1*(S1+cif1)^2) + cumsum(M1)*cif1*PR1 -
                     cumsum(M1*(S1+cif1))*(PR1+cif1))
  G4 = 2*cif1/PR1*(cumsum(M0*cif1^2) + cumsum(M0)*cif1*PR1 -
                     cumsum(M0*cif1)*(PR1+cif1))
  se1 = sqrt(G1+G0+G3+G2-G5-G4)/PR1
  se1[is.nan(se1)] = rev(na.omit(se1))[1]

  V1 = fit10[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G1 = cumsum(M1)*cif0^2 + cumsum(M1*(S0+cif0)^2) -
    2*cif0*cumsum(M1*(S0+cif0))
  G0 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  G3 = cif0^2/PR0^2*sum(M1*(S0+cif0-PR0)^2)
  G2 = cif0^2/PR0^2*sum(M0*(cif0-PR0)^2)
  G5 = 2*cif0/PR0*(cumsum(M1*(S0+cif0)^2) + cumsum(M1)*cif0*PR0 -
                     cumsum(M1*(S0+cif0))*(PR0+cif0))
  G4 = 2*cif0/PR0*(cumsum(M0*cif0^2) + cumsum(M0)*cif0*PR0 -
                     cumsum(M0*cif0)*(PR0+cif0))
  se0 = sqrt(G1+G0+G3+G2-G5-G4)/PR0
  se0[is.nan(se0)] = rev(na.omit(se0))[1]
  cif1 = cif1/PR1
  cif0 = cif0/PR0
  ate = .matchy(cif1,time1,tt)-.matchy(cif0,time0,tt)
  se = sqrt(.matchy(se1,time1,tt)^2+.matchy(se0,time0,tt)^2)
  return(list(time1=time1,time0=time0,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=NULL,cumhaz=cumhaz))
}
