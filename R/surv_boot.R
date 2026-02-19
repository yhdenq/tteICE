#' @title 
#' Calculate standard errors for estimated CIFs and treatment effects
#'
#' @description 
#' This function calculates the standard error for the estimated potential cumulative incidence function
#' and treatment effect. Two methods to calculate the standard error are considered: the asymptotic standard error
#' based on the explicit formula and bootstrapping.
#'
#' @param fit 
#' A fitted object returned by the function \code{tteICE}, \code{surv.tteICE}, or \code{scr.tteICE}.
#'
#' @param nboot 
#' Number of resamplings in the boostrapping method. If \code{nboot} is 0 or 1, then
#' asymptotic standard error based on the explicit form is calculated instead of bootstrapping.
#'
#' @param seed Seed for bootstrapping.
#'
#'
#' @return A list including
#' \describe{
#' \item{time}{Time points in both groups.}
#' \item{cif1}{Estimated cumulative incidence function in the treated group.}
#' \item{cif0}{Estimated cumulative incidence function in the control group.}
#' \item{se1}{Standard error of the estimated cumulative incidence function in the treated group.}
#' \item{se0}{Standard error of the estimated cumulative incidence function in the control group.}
#' \item{ate}{Estimated treatment effect (difference in cumulative incidence functions).}
#' \item{se}{Standard error of the estimated treatment effect.}
#' \item{strategy}{Strategy used.}
#' \item{method}{Estimation method used.}
#' }
#'
#'
#' @seealso 
#' \code{\link[tteICE]{surv.tteICE}}, \code{\link[tteICE]{scr.tteICE}}, 
#' \code{\link[tteICE]{tteICE}}
#'
#' @keywords internal

surv.boot <- function(fit,nboot=0,seed=NULL){
  N = length(fit$A)
  time1 = fit$time1
  time0 = fit$time0
  if ((is.null(time1)&is.null(time1))) {
    time1 = time0 = fit$time
  }
  tt = sort(unique(c(0,time1,time0)))
  cif1 = .matchy(fit$cif1,fit$time1,tt)
  cif0 = .matchy(fit$cif0,fit$time0,tt)
  se1 = .matchy(fit$se1,fit$time1,tt)
  se0 = .matchy(fit$se0,fit$time0,tt)
  ate = cif1-cif0
  se = .matchy(fit$se,fit$time,tt)
  p.val = fit$p.val
  if (nboot>1){
    cif1l = cif0l = te = NULL
    if (is.null(seed)) seed = 0
    set.seed(seed)
    for(b in 1:nboot){
      subset = sample(1:N,replace=TRUE)
      if (fit$dtype=='cmprsk'){
        fitb = surv.tteICE(fit$A,fit$Time,fit$cstatus,fit$strategy,fit$cov1,fit$method,
                           fit$weights,subset=subset,na.rm=FALSE,nboot=-1)
      } else {
        fitb = scr.tteICE(fit$A,fit$Time,fit$status,fit$Time_int,fit$status_int,fit$strategy,fit$cov1,fit$method,
                          fit$weights,subset=subset,na.rm=FALSE,nboot=-1)
      }
      cifb1 = .matchy(fitb$cif1,fitb$time1,tt)
      cifb0 = .matchy(fitb$cif0,fitb$time0,tt)
      cif1l = rbind(cif1l, cifb1)
      cif0l = rbind(cif0l, cifb0)
      te = rbind(te, cifb1-cifb0)
    }
    se1 = apply(cif1l,2,sd,na.rm=TRUE)
    se0 = apply(cif0l,2,sd,na.rm=TRUE)
    se = apply(te,2,sd)
  }
  return(list(time=tt,cif1=cif1,cif0=cif0,ate=ate,se1=se1,se0=se0,se=se,p.val=p.val,
              strategy=fit$strategy,method=fit$method,dtype=fit$dtype))
}
