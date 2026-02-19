#' @title Fit CIFs using principal stratum strategy for competing risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function
#' based on efficient influence functions using principal stratum strategy (competing
#' risks data structure). Cox models are employed for survival models. The estimand is defined in a subpopulation where
#' intercurrent events would never occur regardless of treatment conditions.
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
#' @seealso \code{\link[tteICE]{surv.principal}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @export

surv.principal.eff <- function(A,Time,cstatus,X=NULL){
  n = length(A)
  if (is.null(X)){
    return(surv.principal(A,Time,cstatus))
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
  cumhaz11 = exp(Xb11)%*%t(cumhaz11)
  cumhaz10 = .matchy(c(0,basehaz(fit10,centered=FALSE)$hazard),tt10,tt)
  cumhaz10 = exp(Xb10)%*%t(cumhaz10)
  cumhaz21 = .matchy(c(0,basehaz(fit21,centered=FALSE)$hazard),tt21,tt)
  cumhaz21 = exp(Xb21)%*%t(cumhaz21)
  cumhaz20 = .matchy(c(0,basehaz(fit20,centered=FALSE)$hazard),tt20,tt)
  cumhaz20 = exp(Xb20)%*%t(cumhaz20)
  cumhaz1c = .matchy(c(0,basehaz(fit1c,centered=FALSE)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = .matchy(c(0,basehaz(fit0c,centered=FALSE)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  cumhaz1 = cbind(0,cumhaz11+cumhaz21)[,1:K]
  cumhaz0 = cbind(0,cumhaz10+cumhaz20)[,1:K]
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
  cif11 = t(apply(exp(-cumhaz1)*lam11,1,cumsum))
  cif10 = t(apply(exp(-cumhaz0)*lam10,1,cumsum))
  cif1x = A*ips*t(apply(exp(-cumhaz1)*dMP11,1,cumsum))-
    A*ips*cif11*t(apply(dMP11+dMP21,1,cumsum))+
    A*ips*t(apply(cif11*(dMP11+dMP21),1,cumsum))+cif11
  cif0x = (1-A)*ips*t(apply(exp(-cumhaz0)*dMP10,1,cumsum))-
    (1-A)*ips*cif10*t(apply(dMP10+dMP20,1,cumsum))+
    (1-A)*ips*t(apply(cif10*(dMP10+dMP20),1,cumsum))+cif10
  cif.wo1 = colMeans(cif1x,na.rm=TRUE)
  cif.wo0 = colMeans(cif0x,na.rm=TRUE)
  if.wo1 = t(t(cif1x)-cif.wo1)
  if.wo0 = t(t(cif0x)-cif.wo0)
  cif1x = A*ips*exp(-cumhaz1)*t(apply(dMP11+dMP21,1,cumsum))+1-exp(-cumhaz1)
  cif0x = (1-A)*ips*exp(-cumhaz0)*t(apply(dMP10+dMP20,1,cumsum))+1-exp(-cumhaz0)
  cif.cv1 = colMeans(cif1x,na.rm=TRUE)
  cif.cv0 = colMeans(cif0x,na.rm=TRUE)
  if.cv1 = t(t(cif1x)-cif.cv1)
  if.cv0 = t(t(cif0x)-cif.cv0)
  cif1x = t((t(if.wo1)+cif.wo1)/min(1-cif.cv1+cif.wo1))+
    ((if.cv1-if.wo1)[,ncol(if.cv1)]/min(1-cif.cv1+cif.wo1)^2)%*%t(cif.wo1)
  cif0x = t((t(if.wo0)+cif.wo0)/min(1-cif.cv0+cif.wo0))+
    ((if.cv0-if.wo0)[,ncol(if.cv1)]/min(1-cif.cv0+cif.wo0)^2)%*%t(cif.wo0)
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
  coef11 = fit11$coefficients
  coef10 = fit10$coefficients
  coef21 = fit21$coefficients
  coef20 = fit20$coefficients
  coef = list(coef11=coef11,coef10=coef10,coef21=coef21,coef20=coef20)
  ph11 = cox.zph(fit11)
  ph10 = cox.zph(fit10)
  ph21 = cox.zph(fit21)
  ph20 = cox.zph(fit20)
  ph = list(ph11=ph11,ph10=ph10,ph21=ph21,ph20=ph20)
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,
              coef=coef,ph=ph))
}



