#' @title Fit CIFs using hypothetical strategy (I) for semicompeting risks data
#'
#' @description This function nonparametrically estimates the potential cumulative incidence function
#' using hypothetical strategy (semicompeting risks data structure). The intercurrent event is only
#' permitted under treated if it would occur under control.
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
#' \item{cumhaz}{Baselime cumulative hazards in the survival models.}
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
#' @seealso \code{\link[tteICE]{scr.natural.eff}}, \code{\link[tteICE]{scr.tteICE}}
#'
#'
#' @keywords internal

scr.natural <- function(A,Time,status,Time_int,status_int,weights=rep(1,length(A))){
  Td = Time; Dd = status
  Tr = Time_int; Dr = status_int
  tseq = sort(unique(c(Td[Dd==1],Tr[Dr==1])))
  K = length(tseq)
  haz1.1 = haz2.1 = haz3.1 = haz1t.1 = haz2t.1 = haz3t.1 = 0
  haz1.0 = haz2.0 = haz3.0 = haz1t.0 = haz2t.0 = haz3t.0 = 0
  F1.01 = F2.01 = F3.01 = F1t.01 = F2t.01 = F3t.01 = F3t.n.01 = 0
  F1.00 = F2.00 = F3.00 = F1t.00 = F2t.00 = F3t.00 = F3t.n.00 = 0
  G1.1.01 = G1.2.01 = G1.3.01 = G2.1.01 = G2.2.01 = G2.3.01 = G3.01 = 0
  G1.1.00 = G1.2.00 = G1.3.00 = G2.1.00 = G2.2.00 = G2.3.00 = G3.00 = 0
  G2.4 = 0
  for (k in 1:K){
    t0 = tseq[k]
    Y1.1 = sum(A*weights*(Td>=t0)*(Tr>=t0))
    dhaz1.1 = sum(A*weights*(Td==t0)*Dd*(1-Dr))/Y1.1
    if (Y1.1==0) dhaz1.1 = 0
    Y1.0 = sum((1-A)*weights*(Td>=t0)*(Tr>=t0))
    dhaz1.0 = sum((1-A)*weights*(Td==t0)*Dd*(1-Dr))/Y1.0
    if (Y1.0==0) dhaz1.0 = 0
    Y2.1 = sum(A*weights*(Td>=t0)*(Tr>=t0))
    dhaz2.1 = sum(A*weights*(Tr==t0)*Dr)/Y2.1
    if (Y2.1==0) dhaz2.1 = 0
    Y2.0 = sum((1-A)*weights*(Td>=t0)*(Tr>=t0))
    dhaz2.0 = sum((1-A)*weights*(Tr==t0)*Dr)/Y2.0
    if (Y2.0==0) dhaz2.0 = 0
    Y3.1 = sum(A*weights*(Td>=t0)*(Tr<=t0)*Dr)
    dhaz3.1 = sum(A*weights*(Td==t0)*Dd*Dr)/Y3.1
    if (Y3.1==0) dhaz3.1 = 0
    Y3.0 = sum((1-A)*weights*(Td>=t0)*(Tr<=t0)*Dr)
    dhaz3.0 = sum((1-A)*weights*(Td==t0)*Dd*Dr)/Y3.0
    if (Y3.0==0) dhaz3.0 = 0

    cumhaz = data.frame(time=tseq, cumhaz11=cumsum(dhaz1.1), cumhaz10=cumsum(dhaz1.0),
                       cumhaz21=cumsum(dhaz2.1), cumhaz20=cumsum(dhaz2.0),
                       cumhaz31=cumsum(dhaz3.1), cumhaz30=cumsum(dhaz3.0))
    if (tseq[1]!=0) cumhaz = rbind(0,cumhaz)

    dG1.1.01 = dhaz1.1/Y1.1
    dG1.1.00 = dhaz1.0/Y1.0
    dG1.2.01 = F3t.n.01*dhaz1.1/Y1.1
    dG1.2.00 = F3t.n.00*dhaz1.0/Y1.0
    dG1.3.01 = F3t.n.01^2*dhaz1.1/Y1.1
    dG1.3.00 = F3t.n.00^2*dhaz1.0/Y1.0
    dG2.1.01 = dG2.1.00 = dhaz2.0/Y2.0
    dG2.2.01 = (1-F1t.01-F3t.01)*exp(haz3t.1)*dhaz2.0/Y2.0
    dG2.2.00 = (1-F1t.00-F3t.00)*exp(haz3t.0)*dhaz2.0/Y2.0
    dG2.4 = (1-F1t.01-F3t.01)*(1-F1t.00-F3t.00)*exp(haz3t.1+haz3t.0)*dhaz2.0/Y2.0
    dG2.3.01 = (1-F1t.01-F3t.01)^2*exp(2*haz3t.1)*dhaz2.0/Y2.0
    dG2.3.00 = (1-F1t.00-F3t.00)^2*exp(2*haz3t.0)*dhaz2.0/Y2.0
    dG3.01 = F3t.n.01^2*dhaz3.1/Y3.1
    dG3.00 = F3t.n.00^2*dhaz3.0/Y3.0
    if (Y1.1==0) dG1.1.01=dG1.2.01=dG1.3.01=0
    if (Y1.0==0) dG1.1.00=dG1.2.00=dG1.3.00=0
    if (Y2.0==0) dG2.1.01=dG2.1.00=dG2.2.01=dG2.2.00=dG2.3.01=dG2.3.00=dG2.4=0
    if (Y3.1==0) dG3.01=0
    if (Y3.0==0) dG3.00=0
    G1.1.01 = append(G1.1.01, dG1.1.01)
    G1.1.00 = append(G1.1.00, dG1.1.00)
    G1.2.01 = append(G1.2.01, dG1.2.01)
    G1.2.00 = append(G1.2.00, dG1.2.00)
    G1.3.01 = append(G1.3.01, dG1.3.01)
    G1.3.00 = append(G1.3.00, dG1.3.00)
    G2.1.01 = append(G2.1.01, dG2.1.01)
    G2.1.00 = append(G2.1.00, dG2.1.00)
    G2.2.01 = append(G2.2.01, dG2.2.01)
    G2.2.00 = append(G2.2.00, dG2.2.00)
    G2.3.01 = append(G2.3.01, dG2.3.01)
    G2.3.00 = append(G2.3.00, dG2.3.00)
    G2.4 = append(G2.4, dG2.4)
    G3.01 = append(G3.01, dG3.01)
    G3.00 = append(G3.00, dG3.00)
    F1t.01 = F1t.01 + exp(-haz1t.1-haz2t.0)*dhaz1.1
    F1t.00 = F1t.00 + exp(-haz1t.0-haz2t.0)*dhaz1.0
    F2t.01 = F2t.01 + exp(-haz1t.1-haz2t.0)*dhaz2.0
    F2t.00 = F2t.00 + exp(-haz1t.0-haz2t.0)*dhaz2.0
    F3t.n.01 = F3t.n.01 + exp(-haz1t.1-haz2t.0+haz3t.1)*dhaz2.0
    F3t.n.00 = F3t.n.00 + exp(-haz1t.0-haz2t.0+haz3t.0)*dhaz2.0
    haz3.1 = append(haz3.1, haz3t.1)
    haz3.0 = append(haz3.0, haz3t.0)
    F3t.01 = F2t.01 - F3t.n.01*exp(-haz3t.1)
    F3t.00 = F2t.00 - F3t.n.00*exp(-haz3t.0)
    #F3t.01 = F3t.01 + (F2t.01-F3t.01)*dhaz3.1
    #F3t.00 = F3t.00 + (F2t.00-F3t.00)*dhaz3.0
    haz1t.1 = haz1t.1 + dhaz1.1
    haz1t.0 = haz1t.0 + dhaz1.0
    haz2t.1 = haz2t.1 + dhaz2.1
    haz2t.0 = haz2t.0 + dhaz2.0
    haz3t.1 = haz3t.1 + dhaz3.1
    haz3t.0 = haz3t.0 + dhaz3.0
    haz1.1 = append(haz1.1, haz1t.1)
    haz1.0 = append(haz1.0, haz1t.0)
    haz2.1 = append(haz2.1, haz2t.1)
    haz2.0 = append(haz2.0, haz2t.0)
    F1.01 = append(F1.01, F1t.01)
    F1.00 = append(F1.00, F1t.00)
    F2.01 = append(F2.01, F2t.01)
    F2.00 = append(F2.00, F2t.00)
    F3.01 = append(F3.01, F3t.01)
    F3.00 = append(F3.00, F3t.00)
  }
  Fhaz.01 = 1-F1.01-F3.01
  Fhaz.00 = 1-F1.00-F3.00
  G1.01 = Fhaz.01^2*cumsum(G1.1.01)-
    2*Fhaz.01*exp(-haz3.1)*cumsum(G1.2.01)+
    exp(-haz3.1)^2*cumsum(G1.3.01)
  G1.00 = Fhaz.00^2*cumsum(G1.1.00)-
    2*Fhaz.00*exp(-haz3.0)*cumsum(G1.2.00)+
    exp(-haz3.0)^2*cumsum(G1.3.00)
  G2.01 = Fhaz.01^2*cumsum(G2.1.01)-
    2*Fhaz.01*exp(-haz3.1)*cumsum(G2.2.01)+
    exp(-2*haz3.1)*cumsum(G2.3.01)
  G2.de = (Fhaz.00-Fhaz.01)^2*cumsum(G2.1.00)+
    2*(Fhaz.00-Fhaz.01)*(exp(-haz3.1)*cumsum(G2.2.01)-exp(-haz3.0)*cumsum(G2.2.00))+
    exp(-2*haz3.1)*cumsum(G2.3.01)+exp(-2*haz3.0)*cumsum(G2.3.00)-
    2*exp(-haz3.1-haz3.0)*cumsum(G2.4)
  G2.00 = Fhaz.00^2*cumsum(G2.1.00)-
    2*Fhaz.00*exp(-haz3.0)*cumsum(G2.2.00)+
    exp(-2*haz3.0)*cumsum(G2.3.00)
  G3.01 = exp(-2*haz3.1)*cumsum(G3.01)
  G3.00 = exp(-2*haz3.0)*cumsum(G3.00)
  G1.01[is.na(G1.01)] = 0
  G1.00[is.na(G1.00)] = 0
  G2.01[is.na(G2.01)] = 0
  G2.00[is.na(G2.00)] = 0
  G2.de[is.na(G2.de)] = 0
  cif0 = F1.00+F3.00
  cif2 = F1.01+F3.01
  ate = cif2 - cif0
  se0 = sqrt(G1.00+G2.00+G3.00)
  se2 = sqrt(G1.01+G2.01+G3.01)
  se = sqrt(G1.01+G1.00+G2.de+G3.01+G3.00)
  if (tseq[1]==0){
    tseq = tseq[-1]; cif2 = cif2[-1]; cif0 = cif0[-1]
    se2 = se2[-1]; se0 = se0[-1]; ate = ate[-1]; se = se[-1]
  }
  tseq = c(0,tseq)
  surv_diff = survdiff(Surv(Td,Dd)~A)
  p = pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail=FALSE)
  return(list(time1=tseq,time0=tseq,cif1=cif2,cif0=cif0,se1=se2,se0=se0,
              time=tseq,ate=ate,se=se,p.val=p,cumhaz=cumhaz))
}
