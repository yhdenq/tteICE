#' @title Fit CIFs using hypothetical strategy (I) for competing risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function
#' based on efficient influence functions using hypothetical strategy (competing risks
#' data structure). Cox models are employed for survival models. The intercurrent event
#' is only permitted under treated if it would occur under control.
#'
#' @param A Treatment indicator, 1 for treatment and 0 for control.
#'
#' @param Time Time to event.
#'
#' @param cstatus Indicator of event, 1 for the primary event, 2 for the intercurrent event, 0 for censoring.
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
#' \item{cumhaz}{Baseline cumulative hazards in the working Cox models.}
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
#' @seealso \code{\link[tteICE]{surv.natural}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @keywords internal

surv.natural.eff <- function(A,Time,cstatus,X=NULL){
  n = length(A)
  if (is.null(X)){
    return(surv.natural(A,Time,cstatus))
  }
  X = as.matrix(scale(X))
  ips = .ipscore(A,X)
  fit11 = coxph(Surv(Time,cstatus==1)~X, subset=(A==1))
  fit10 = coxph(Surv(Time,cstatus==1)~X, subset=(A==0))
  fit21 = coxph(Surv(Time,cstatus>1)~X, subset=(A==1))
  fit20 = coxph(Surv(Time,cstatus>1)~X, subset=(A==0))
  fit1c = coxph(Surv(Time,cstatus==0)~X, subset=(A==1))
  fit0c = coxph(Surv(Time,cstatus==0)~X, subset=(A==0))
  tt11 = c(0,basehaz(fit11)$time)
  tt10 = c(0,basehaz(fit10)$time)
  tt21 = c(0,basehaz(fit21)$time)
  tt20 = c(0,basehaz(fit20)$time)
  tt = sort(unique(c(tt11,tt10,tt21,tt20)))
  K = length(tt)
  Xb11 = X%*%fit11$coefficients
  Xb10 = X%*%fit10$coefficients
  Xb21 = X%*%fit21$coefficients
  Xb20 = X%*%fit20$coefficients
  Xb1c = X%*%fit1c$coefficients
  Xb0c = X%*%fit0c$coefficients
  cumhaz11 = .matchy(c(0,basehaz(fit11,centered=FALSE)$hazard),tt11,tt)
  cumhaz10 = .matchy(c(0,basehaz(fit10,centered=FALSE)$hazard),tt10,tt)
  cumhaz21 = .matchy(c(0,basehaz(fit21,centered=FALSE)$hazard),tt21,tt)
  cumhaz20 = .matchy(c(0,basehaz(fit20,centered=FALSE)$hazard),tt20,tt)
  cumhaz = data.frame(time=tt,cumhaz11=cumhaz11,cumhaz10=cumhaz10,cumhaz21=cumhaz21,cumhaz20=cumhaz20)
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = .matchy(c(0,basehaz(fit1c,centered=FALSE)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = .matchy(c(0,basehaz(fit0c,centered=FALSE)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  cumhaz1 = cbind(0,cumhaz11+cumhaz21)[,1:K]
  cumhaz0 = cbind(0,cumhaz10+cumhaz20)[,1:K]
  cumhaz2 = cbind(0,cumhaz11+cumhaz20)[,1:K]
  dN1 = sapply(tt, function(l) (Time==l)*(cstatus==1))
  dN2 = sapply(tt, function(l) (Time==l)*(cstatus>1))
  Y = sapply(tt, function(l) as.numeric(Time>=l))
  lam11 = t(apply(cbind(0,cumhaz11),1,diff))
  lam10 = t(apply(cbind(0,cumhaz10),1,diff))
  lam21 = t(apply(cbind(0,cumhaz21),1,diff))
  lam20 = t(apply(cbind(0,cumhaz20),1,diff))
  S1 = exp(-cumhaz1-cbind(0,cumhaz1c)[,1:K])
  S0 = exp(-cumhaz0-cbind(0,cumhaz0c)[,1:K])
  dMP11 = (dN1-Y*lam11)/S1
  dMP21 = (dN2-Y*lam21)/S1
  dMP10 = (dN1-Y*lam10)/S0
  dMP20 = (dN2-Y*lam20)/S0
  cif1 = t(apply(exp(-cumhaz2)*lam11,1,cumsum))
  cif0 = t(apply(exp(-cumhaz0)*lam10,1,cumsum))
  cif1x = A*ips*t(apply((exp(-cumhaz2)+cif1)*dMP11,1,cumsum))-
    A*ips*cif1*t(apply(dMP11,1,cumsum))+
    (1-A)*ips*t(apply(cif1*dMP20,1,cumsum))-
    (1-A)*ips*cif1*t(apply(dMP20,1,cumsum))+cif1
  cif0x = (1-A)*ips*t(apply(exp(-cumhaz0)*dMP10,1,cumsum))-
    (1-A)*ips*cif0*t(apply(dMP10+dMP20,1,cumsum))+
    (1-A)*ips*t(apply(cif0*(dMP10+dMP20),1,cumsum))+cif0
  cif1 = colMeans(cif1x,na.rm=TRUE)
  cif0 = colMeans(cif0x,na.rm=TRUE)
  se1 = apply(cif1x,2,sd,na.rm=TRUE)/sqrt(n)
  se0 = apply(cif0x,2,sd,na.rm=TRUE)/sqrt(n)
  ate = cif1-cif0
  se = apply(cif1x-cif0x,2,sd,na.rm=TRUE)/sqrt(n)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Ti = (tt<0.99*max(tt))
  Tt = sum((cif1-cif0)*diff(c(0,tt))*Ti)
  IFt = colSums(t(eif1-eif0)*diff(c(0,tt))*Ti,na.rm=TRUE)
  Vt = sd(IFt,na.rm=TRUE)/sqrt(n)
  p = 2*pnorm(-abs(Tt/Vt))
  coef11 = fit11$coefficients / attr(X,"scaled:scale")
  coef10 = fit10$coefficients / attr(X,"scaled:scale")
  coef21 = fit21$coefficients / attr(X,"scaled:scale")
  coef20 = fit20$coefficients / attr(X,"scaled:scale")
  coef = list(coef11=coef11,coef10=coef10,coef21=coef21,coef20=coef20)
  ph11 = cox.zph(fit11, terms=FALSE)
  ph10 = cox.zph(fit10, terms=FALSE)
  ph21 = cox.zph(fit21, terms=FALSE)
  ph20 = cox.zph(fit20, terms=FALSE)
  ph = list(ph11=ph11,ph10=ph10,ph21=ph21,ph20=ph20)
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,
              coef=coef,ph=ph,cumhaz=cumhaz))
}
