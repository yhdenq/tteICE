#' @title Fit CIFs using hypothetical strategy (I) for competing risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using hypothetical strategy (competing risks data structure). The intercurrent event is only
#' permitted under treated if it would occur under control.
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
#' \eqn{\lambda_1(t; w)} remains unchanged. Specifically, we envision that the intercurrent events that
#' occurred when individuals were assigned to test drugs were only permitted if these intercurrent events
#' would have also occurred if these individuals had been assigned to the placebo. In this hypothetical
#' scenario, when assigned to placebo, individuals would be equally likely to experience intercurrent
#' events as they are assigned to placebo in the real-world trial in terms of the hazards; when assigned
#' to test drug, the hazard of intercurrent events would be identical to that if assigned to placebo in
#' the real-world trial. That is, \eqn{\lambda_2'(t;0) = \lambda_2'(t;1) = \lambda_2(t;0)}. The treatment
#' effect corresponds to the natural direct effect with the hazard of intercurrent events set at
#' the level under control.
#' }
#'
#' @seealso \code{\link[tteICE]{surv.natural.eff}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.natural <- function(A,Time,cstatus,weights=rep(1,length(A))){
  n = length(A)
  s1 = (A==1); n1 = sum(s1)
  s0 = (A==0); n0 = sum(s0)
  fit11 = survfitKM(factor(rep(1,n1)), Surv(Time,cstatus==1)[s1], weights=weights[s1])
  fit10 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus==1)[s0], weights=weights[s0])
  fit20 = survfitKM(factor(rep(0,n0)), Surv(Time,cstatus>1)[s0], weights=weights[s0])
  time = c(0, unique(sort(c(fit11$time,fit10$time,fit20$time))))
  fit11 = .matchy(rbind(0,cbind(fit11$cumhaz,fit11$std.err)),c(0,fit11$time),time)
  fit10 = .matchy(rbind(0,cbind(fit10$cumhaz,fit10$std.err)),c(0,fit10$time),time)
  fit20 = .matchy(rbind(0,cbind(fit20$cumhaz,fit20$std.err)),c(0,fit20$time),time)
  dcif1 = exp(-fit11[,1]-fit20[,1])*diff(c(0,fit11[,1]))
  dcif0 = exp(-fit10[,1]-fit20[,1])*diff(c(0,fit10[,1]))
  cif1 = cumsum(dcif1)
  cif0 = cumsum(dcif0)
  V1 = fit11[,2]^2
  V0 = fit20[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  V0[is.infinite(V0)] = max(V0[!is.infinite(V0)])
  M1 = diff(c(0,V1))
  M0 = diff(c(0,V0))
  G11 = cumsum(M1)*cif1^2 + cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1)^2) -
    2*cif1*cumsum(M1*(exp(-fit11[,1]-fit20[,1])+cif1))
  G01 = cumsum(M0)*cif1^2 + cumsum(M0*cif1^2) - 2*cif1*cumsum(M0*cif1)
  se1 = sqrt(G11+G01)
  V1 = fit10[,2]^2
  V1[is.infinite(V1)] = max(V1[!is.infinite(V1)])
  M1 = diff(c(0,V1))
  G10 = cumsum(M1)*cif0^2 + cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0)^2) -
    2*cif0*cumsum(M1*(exp(-fit10[,1]-fit20[,1])+cif0))
  G00 = cumsum(M0)*cif0^2 + cumsum(M0*cif0^2) - 2*cif0*cumsum(M0*cif0)
  se0 = sqrt(G10+G00)
  dcif = cif1-cif0
  G0 = cumsum(M0)*dcif^2 + cumsum(M0*dcif^2) - 2*dcif*cumsum(M0*dcif)
  se = sqrt(G11+G10+G0)
  surv_diff = survdiff(Surv(Time,cstatus==1)~A)
  p = pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail=FALSE)
  ate = cif1-cif0
  return(list(time1=time,time0=time,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=time,ate=ate,se=se,p.val=p))
}
