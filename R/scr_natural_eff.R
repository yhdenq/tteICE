#' @title Fit CIFs using hypothetical strategy (I) for semicompeting risks data, based on efficient influence functions
#'
#' @description This function estimates the potential cumulative incidence function
#' based on efficient influence functions using hypothetical strategy (semicompeting risks
#' data structure). Cox models are employed for survival models. The intercurrent event
#' is only permitted under treated if it would occur under control.
#'
#' @param A Treatment indicator, 1 for treatment and 0 for control.
#'
#' @param Time Time to the primary (terminal) event.
#'
#' @param status Indicator of the primary (terminal) event, 1 for event and 0 for censoring.
#'
#' @param Time_int Time to the intercurrent event.
#'
#' @param status_int Indicator of the intercurrent event, 1 for event and 0 for censoring.
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
#' effect corresponds to the natural direct effect, with the hazard of intercurrent events set at
#' the level under control. Markovness is assumed in estimation.
#' }
#' 
#' @seealso \code{\link[tteICE]{scr.natural}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.natural.eff <- function(A,Time,status,Time_int,status_int,X=NULL){
  n = length(A)
  if (is.null(X)) {
    return(scr.natural(A,Time,status,Time_int,status_int))
  }
  X = as.matrix(scale(X))
  df = data.frame(id=1:n, Td=Time, Dd=status, Tr=Time_int, Dr=status_int, X=X, A=A)
  id <- df$id
  mg = tmerge(data1=df, data2=df, id=id, event=.event(df$Td,df$Dd))
  mg = tmerge(mg, df, id=id, Revent = .tdc(df$Tr))
  xvars = grep("^X", names(mg), value=TRUE)
  cox_formula = reformulate(c("Revent",xvars), response="Surv(tstart,tstop,event)")
  tt = sort(unique(c(Time,Time_int)))
  l = length(tt)
  fit11 = coxph(cox_formula, data=mg, cluster=mg$id, subset=(A==1))
  time1 = basehaz(fit11,centered=FALSE)$time
  lamd = diff(c(0,basehaz(fit11,centered=FALSE)$hazard))
  lam_d = .matchy(lamd, time1, tt, TRUE)
  lam_od1 = sapply(1:l, function(t) lam_d[t]*exp( X%*%fit11$coefficients[-1]))
  lam_ord1 = lam_od1 * exp(fit11$coefficients[1])
  fit10 = coxph(cox_formula, data=mg, cluster=mg$id, subset=(A==0))
  time0 = basehaz(fit10,centered=FALSE)$time
  lamd = diff(c(0,basehaz(fit10,centered=FALSE)$hazard))
  lam_d = .matchy(lamd, time0, tt, TRUE)
  lam_od0 = sapply(1:l, function(t) lam_d[t]*exp(X%*%fit10$coefficients[-1]))
  lam_ord0 = lam_od0 * exp(fit10$coefficients[1])
  fit21 = coxph(Surv(Time_int,status_int)~X, subset=(A==1))
  time1 = basehaz(fit21,centered=FALSE)$time
  lamr = diff(c(0,basehaz(fit21,centered=FALSE)$hazard))
  lam_r = .matchy(lamr, time1, tt, TRUE)
  lam_or1 = sapply(1:l, function(t) lam_r[t]*exp(X%*%fit21$coefficients))
  fit20 = coxph(Surv(Time_int,status_int)~X, subset=(A==0))
  time0 = basehaz(fit20,centered=FALSE)$time
  lamr = diff(c(0,basehaz(fit20,centered=FALSE)$hazard))
  lam_r = .matchy(lamr, time0, tt, TRUE)
  lam_or0 = sapply(1:l, function(t) lam_r[t]*exp(X%*%fit20$coefficients))
  fit1c = coxph(Surv(Time,status==0)~X, subset=(A==1))
  time1 = basehaz(fit1c,centered=FALSE)$time
  lamc = diff(c(0,basehaz(fit1c,centered=FALSE)$hazard))
  lam_c = .matchy(lamc, time1, tt, TRUE)
  lam_c1 = sapply(1:l, function(t) lam_c[t]*exp(X%*%fit1c$coefficients))
  fit0c = coxph(Surv(Time,status==0)~X, subset=(A==0))
  time0 = basehaz(fit0c,centered=FALSE)$time
  lamc = diff(c(0,basehaz(fit0c,centered=FALSE)$hazard))
  lam_c = .matchy(lamc, time0, tt, TRUE)
  lam_c0 = sapply(1:l, function(t) lam_c[t]*exp(X%*%fit0c$coefficients))
  ips = .ipscore(A,X)

  # observable incidence
  lam_od_A = A*lam_od1 + (1-A)*lam_od0
  lam_or_A = A*lam_or1 + (1-A)*lam_or0
  lam_ord_A = A*lam_ord1 + (1-A)*lam_ord0
  Lam_o_A = t(apply(lam_od_A+lam_or_A, 1, cumsum))
  Lam_or_A = t(apply(lam_ord_A, 1, cumsum))
  lam_c_A = A*lam_c1 + (1-A)*lam_c0
  Lam_c_A = t(apply(lam_c_A, 1, cumsum))
  SC = exp(-Lam_c_A)

  dF_or_A = exp(-Lam_o_A)*lam_or_A
  dF_od_A = exp(-Lam_o_A)*lam_od_A
  F_or_A = t(apply(dF_or_A, 1, cumsum))
  F_od_A = t(apply(dF_od_A, 1, cumsum))
  dF_or__A = t(apply(dF_or_A*exp(Lam_or_A), 1, cumsum))*exp(-Lam_or_A)
  dF_ord_A = dF_or__A*lam_ord_A
  F_ord_A = t(apply(dF_ord_A, 1, cumsum))
  F_A = colMeans(F_od_A+F_ord_A)

  a = c(1,0,1)
  lam_od_a = a[1]*lam_od1 + (1-a[1])*lam_od0
  lam_or_a = a[2]*lam_or1 + (1-a[2])*lam_or0
  lam_ord_a = a[3]*lam_ord1 + (1-a[3])*lam_ord0
  Lam_o_a = t(apply(lam_od_a+lam_or_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a, 1, cumsum))

  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  Fd = F_od_a + F_ord_a

  # od
  Y_o = sapply(tt, function(t) (Time>=t)*(Time_int>=t))
  PY_o = (1-F_or_A-F_od_A) * SC
  dM_od = sapply(tt, function(t) (Time==t)*status*(Time_int>=t)) - Y_o*lam_od_A
  dQ_od = dM_od * (A==a[1])*ips /PY_o
  dM_or = sapply(tt, function(t) (Time_int==t)*status_int*(Time>=t)) - Y_o*lam_or_A
  dQ_or = dM_or * (A==a[2])*ips / PY_o
  dQ_od[PY_o==0] = 0
  dQ_or[PY_o==0] = 0
  Q_or = t(apply(dQ_or,1,cumsum))
  Q_od = t(apply(dQ_od,1,cumsum))
  G1_od = dQ_od - (Q_od+Q_or)*lam_od_a
  G_od = t(apply(G1_od*exp(-Lam_o_a), 1, cumsum))
  # ord
  Y_or = sapply(tt, function(t) (Time>=t)*(Time_int<=t)*status_int)
  PY_or = (F_or_A - F_ord_A) * SC
  dM_ord = sapply(tt, function(t) (Time==t)*status*status_int*(Time_int<=t)) - Y_or*lam_ord_A
  dQ_ord = dM_ord * (A==a[3])*ips / PY_or
  dQ_ord[PY_or==0] = 0
  Q_ord = t(apply(dQ_ord,1,cumsum))
  G1_or = dQ_or - (Q_od+Q_or)*lam_or_a
  G1_or = t(apply(G1_or*exp(Lam_or_a-Lam_o_a), 1, cumsum))
  G2_or = dQ_ord - Q_ord*lam_ord_a
  G2_or1 = t(apply(exp(Lam_or_a-Lam_o_a)*lam_or_a, 1, cumsum))
  G2_or2 = G2_or1 * G2_or
  G3_or = t(apply(exp(Lam_or_a-Lam_o_a)*lam_or_a*Q_ord, 1, cumsum))
  G_ord = t(apply(((G1_or+G3_or)*lam_ord_a+G2_or2)*exp(-Lam_or_a), 1, cumsum))

  EIF1 = G_od + G_ord
  Feff1 = Fd + EIF1

  a = c(0,0,0)
  lam_od_a = a[1]*lam_od1 + (1-a[1])*lam_od0
  lam_or_a = a[2]*lam_or1 + (1-a[2])*lam_or0
  lam_ord_a = a[3]*lam_ord1 + (1-a[3])*lam_ord0
  Lam_o_a = t(apply(lam_od_a+lam_or_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a, 1, cumsum))

  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  Fd = F_od_a + F_ord_a

  # od
  Y_o = sapply(tt, function(t) (Time>=t)*(Time_int>=t))
  PY_o = (1-F_or_A-F_od_A) * SC
  dM_od = sapply(tt, function(t) (Time==t)*status*(Time_int>=t)) - Y_o*lam_od_A
  dQ_od = dM_od * (A==a[1])*ips /PY_o
  dM_or = sapply(tt, function(t) (Time_int==t)*status_int*(Time>=t)) - Y_o*lam_or_A
  dQ_or = dM_or * (A==a[2])*ips / PY_o
  dQ_od[PY_o==0] = 0
  dQ_or[PY_o==0] = 0
  Q_or = t(apply(dQ_or,1,cumsum))
  Q_od = t(apply(dQ_od,1,cumsum))
  G1_od = dQ_od - (Q_od+Q_or)*lam_od_a
  G_od = t(apply(G1_od*exp(-Lam_o_a), 1, cumsum))
  # ord
  Y_or = sapply(tt, function(t) (Time>=t)*(Time_int<=t)*status_int)
  PY_or = (F_or_A - F_ord_A) * SC
  dM_ord = sapply(tt, function(t) (Time==t)*status*status_int*(Time_int<=t)) - Y_or*lam_ord_A
  dQ_ord = dM_ord * (A==a[3])*ips / PY_or
  dQ_ord[PY_or==0] = 0
  Q_ord = t(apply(dQ_ord,1,cumsum))
  G1_or = dQ_or - (Q_od+Q_or)*lam_or_a
  G1_or = t(apply(G1_or*exp(Lam_or_a-Lam_o_a), 1, cumsum))
  G2_or = dQ_ord - Q_ord*lam_ord_a
  G2_or1 = t(apply(exp(Lam_or_a-Lam_o_a)*lam_or_a, 1, cumsum))
  G2_or2 = G2_or1 * G2_or
  G3_or = t(apply(exp(Lam_or_a-Lam_o_a)*lam_or_a*Q_ord, 1, cumsum))
  G_ord = t(apply(((G1_or+G3_or)*lam_ord_a+G2_or2)*exp(-Lam_or_a), 1, cumsum))

  EIF0 = G_od + G_ord
  Feff0 = Fd + EIF0

  cif1 = apply(Feff1, 2, mean, na.rm=TRUE)
  cif0 = apply(Feff0, 2, mean, na.rm=TRUE)
  se1 = apply(Feff1, 2, sd, na.rm=TRUE) / sqrt(n)
  se0 = apply(Feff0, 2, sd, na.rm=TRUE) / sqrt(n)
  ate = cif1 - cif0
  se = apply(Feff1-Feff0, 2, sd, na.rm=TRUE) / sqrt(n)
  Ti = (tt<0.99*max(tt))
  Tt = sum((cif1-cif0)*diff(c(0,tt))*Ti)
  IFt = colSums(t(EIF1-EIF0)*diff(c(0,tt))*Ti)
  Vt = sd(IFt, na.rm=TRUE)/sqrt(n)
  p = 2*pnorm(-abs(Tt/Vt))
  if (tt[1]>0) {
    tt = c(0,tt); cif1 = c(0,cif1); cif0 = c(0,cif0)
    se1 = c(0,se1); se0 = c(0,se0); ate = c(0,ate); se = c(0,se)
  }
  coef11 = fit11$coefficients
  coef10 = fit10$coefficients
  coef21 = fit21$coefficients
  coef20 = fit20$coefficients
  coef = list(coef11=coef11,coef10=coef10,coef21=coef21,coef20=coef20)
  ph11 = cox.zph(fit11, terms=FALSE)
  ph10 = cox.zph(fit10, terms=FALSE)
  ph21 = cox.zph(fit21, terms=FALSE)
  ph20 = cox.zph(fit20, terms=FALSE)
  ph = list(ph11=ph11,ph10=ph10,ph21=ph21,ph20=ph20)

  return(list(time1=tt,time0=tt,cif1=cif1,cif0=cif0,se1=se1,se0=se0,
              time=tt,ate=ate,se=se,p.val=p,
              coef=coef,ph=ph))
}
