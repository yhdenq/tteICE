#' @title Fit CIFs using composite variable strategy for competing risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function based on
#' efficient influence functions using composite variable strategy (competing risks data structure).
#' Cox models are employed for survival models. This strategy adopts the first occurrence of either
#' the intermediate or primary event as the event of interest.
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
#' @seealso \code{\link[tteICE]{surv.composite}}, \code{\link[tteICE]{surv.tteICE}}
#'
#'
#' @export

surv.composite.eff <- function(A,Time,cstatus,X=NULL){
  n = length(A)
  if (is.null(X)){
    return(surv.composite(A,Time,cstatus))
  } 
  X = as.matrix(scale(X))
  ips = .ipscore(A,X)
  fit1 = coxph(Surv(Time,cstatus>0)~X, subset=(A==1))
  fit0 = coxph(Surv(Time,cstatus>0)~X, subset=(A==0))
  fit1c = coxph(Surv(Time,cstatus==0)~X, subset=(A==1))
  fit0c = coxph(Surv(Time,cstatus==0)~X, subset=(A==0))
  tt1 = c(0,basehaz(fit1)$time)
  tt0 = c(0,basehaz(fit0)$time)
  tt = sort(unique(c(tt1,tt0)))
  K = length(tt)
  Xb1 = X%*%fit1$coefficients
  Xb0 = X%*%fit0$coefficients
  Xb1c = X%*%fit1c$coefficients
  Xb0c = X%*%fit0c$coefficients
  cumhaz1 = .matchy(c(0,basehaz(fit1,centered=FALSE)$hazard),tt1,tt)
  cumhaz1 = exp(Xb1)%*%t(cumhaz1)
  cumhaz0 = .matchy(c(0,basehaz(fit0,centered=FALSE)$hazard),tt0,tt)
  cumhaz0 = exp(Xb0)%*%t(cumhaz0)
  cumhaz1c = .matchy(c(0,basehaz(fit1c,centered=FALSE)$hazard),c(0,basehaz(fit1c)$time),tt)
  cumhaz0c = .matchy(c(0,basehaz(fit0c,centered=FALSE)$hazard),c(0,basehaz(fit0c)$time),tt)
  cumhaz1c = exp(Xb1c)%*%t(cumhaz1c)
  cumhaz0c = exp(Xb0c)%*%t(cumhaz0c)
  dN = sapply(tt, function(l) (Time==l)*(cstatus>0))
  Y = sapply(tt, function(l) as.numeric(Time>=l))
  lam1 = t(apply(cbind(0,cumhaz1),1,diff))
  lam0 = t(apply(cbind(0,cumhaz0),1,diff))
  S1 = cbind(1,exp(-cumhaz1-cumhaz1c))[,1:K]
  S0 = cbind(1,exp(-cumhaz0-cumhaz0c))[,1:K]
  dMP1 = (dN-Y*lam1)/S1
  dMP0 = (dN-Y*lam0)/S0
  cif1x = A*ips*exp(-cumhaz1)*t(apply(dMP1,1,cumsum))+1-exp(-cumhaz1)
  cif0x = (1-A)*ips*exp(-cumhaz0)*t(apply(dMP0,1,cumsum))+1-exp(-cumhaz0)
  cif1 = colMeans(cif1x,na.rm=TRUE)
  cif0 = colMeans(cif0x,na.rm=TRUE)
  se1 = apply(cif1x,2,sd,na.rm=TRUE)/sqrt(n)
  se0 = apply(cif0x,2,sd,na.rm=TRUE)/sqrt(n)
  ate = cif1-cif0
  se = sqrt(se1^2+se0^2)
  eif1 = t(t(cif1x)-cif1)
  eif0 = t(t(cif0x)-cif0)
  Ti = (tt<0.99*max(tt))
  Tt = sum((cif1-cif0)*diff(c(0,tt))*Ti)
  IFt = colSums(t(eif1-eif0)*diff(c(0,tt))*Ti,na.rm=TRUE)
  Vt = sd(IFt,na.rm=TRUE)/sqrt(n)
  p = 2*pnorm(-abs(Tt/Vt))
  coef1 = fit1$coefficients
  coef0 = fit0$coefficients
  coef = list(coef1=coef1,coef0=coef0)
  ph1 = cox.zph(fit1)
  ph0 = cox.zph(fit0)
  ph = list(ph1=ph1,ph0=ph0)
  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,
              coef=coef,ph=ph))
}
